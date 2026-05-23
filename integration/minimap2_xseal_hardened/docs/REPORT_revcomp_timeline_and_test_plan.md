# XSeal × minimap2 集成：双重互补时间线报告与修订测试计划

> 状态：Harvester SIMD 主路径已于 2026-05 修复（去掉 `andnot`、LUT 与 `xseal_syncmer.hpp` 对齐）。  
> 本文档说明 **bug 发生时互补在何时发生**、与扫描阶段的区别，以及修订后的 benchmark 计划。

---

## 1. 集成流水线中「互补」发生的时间点

对 **minimap2_xseal**（minimizer 路径）一次完整建索引，canonical 相关逻辑按时间顺序如下：

| 阶段 | 时刻 | 组件 | 互补次数 | 修复前 | 修复后 |
|------|------|------|----------|--------|--------|
| T0 | 读入参考序列 | `fill_sequences_from_fasta` | 0 | 无 | 无 |
| T1 | 解析 FASTA | `XsealFrontendDistributor::process_file` | 0 | 无 | 无 |
| T2 | 2-bit 编码 | `XSealEncoder::encode_chunk` | 0 | 无 | 无 |
| T3 | **窗口扫描（选命中位置）** | `XSealMinimizer::scan` | **1×（LUT）** | ✅ 正确 | ✅ 正确 |
| T4 | **64-bit 收割（写索引键）** | `xseal_harvest_kmers` SIMD 主循环 | **bug 时 2×** | ❌ 双重互补 + 错误 LUT | ✅ 仅 LUT 1× |
| T4′ | 收割尾循环 | `xseal_harvest_kmers` scalar tail | **1×（`~c`）** | ✅ 正确 | ✅ 正确 |
| T5 | 写入索引 | `mm_idx_add`（键 = Murmur3(canonical)） | 0 | — | — |

**结论（回答「两次互补都发生在什么时候」）：**

- **双重互补只发生在 T4**：修复前的 `xseal_harvest_kmers` **SIMD 主循环**（约 `n - (n%4)` 条 hit，占 hg38 的 ~99%+）。
- **T3 扫描阶段从未双重互补**：`XSealMinimizer` / `XSealSyncmer::extract_8` 仅通过 LUT 做 **一次** 碱基互补（`V_RC_LUT_L`/`V_RC_LUT_H` 为 `0xF0…` / `0x0B…`）。
- **T4′ 尾循环** 始终只有一次 `~c & k_mask`（与 minimap2 语义一致），但体量极小，**无法纠正** T4 的错误。

因此：之前 **distinct 暴跌（~19M）** 与 **singleton ~49%（在 harvester 修复后仍 ~52%）** 分两阶段：

1. **T4 harvester bug**（双重互补 + 错误 LUT）→ distinct ~19M、错误 singleton。
2. **T3 minimizer store bug**（`261ba83`：shuffle 后误用 `temp_p[i]` 而非 `temp_p[0..popcount-1]`）→ 修复后 singleton **~94%**，与 native 一致。详见 `docs/BENCHMARK_hg38_report_2026-05-18.md`。

---

## 2. 漏洞 A：双重互补（仅 Harvester SIMD 主循环，已修复）

### 修复前代码路径（`xseal_kmer_harvester.hpp`，`i + 4 <= n`）

```
vf = extract_kmer(2-bit)           // 前向 k-mer
vr = (~vf) & k_mask                // 第 1 次互补：按位 NOT
vr = xseal_revcomp_avx2_fast(vr)  // 第 2 次互补：LUT 查表交换碱基
canonical = min(vf, vr)
hash = Murmur3(canonical)
```

数学上：第 2 次互补抵消第 1 次 → **仅碱基顺序反转，未做真正反向互补** → canonical 与 Murmur3 键空间塌缩。

### 修复后

```
vf = extract_kmer(2-bit)
vr = xseal_revcomp_avx2_fast(vf)   // 仅 LUT 互补 1 次 + 字节反序 + r_shift
canonical = min(vf, vr)
```

与 **标量尾循环**（`~c & k_mask` + 位置反序）语义一致。

---

## 3. 漏洞 B：错误 LUT（仅 Harvester，已修复）

| 模块 | `V_RC_LUT_L` | `V_RC_LUT_H` | 结果 |
|------|----------------|----------------|------|
| `xseal_syncmer.hpp` / `xseal_minimizer.hpp` | `0xF0, 0xB0, …` | `0x0F, 0x0B, …` | 高低半字节各存互补碱基 ✅ |
| Harvester（修复前） | `15, 11, 7, 3…` | **与 L 相同** | `OR` 后高半字节恒 0 ❌ |

修复后 Harvester 与 Syncmer 使用同一套 LUT 定义。

---

## 4. 与「单例率」指标的关系（修复后应重新对比）

| 指标 | 修复前（harvester bug） | 修复后预期 |
|------|-------------------------|------------|
| `distinct_minimizers` | ~19M（假碰撞） | 应接近 native 量级（~3–4×10⁸，仍可能因 Murmur3≠hash64 有偏差） |
| `singleton_pct` minimizer | ~49%（异常低） | 应与 native ~94% **可比**（需统一索引键为 `hash64` 才严格可比） |
| `singleton_pct` syncmer | ~91% | 修复后仍可能略高于 minimizer（采样算法不同），但不应再因 T4 bug 失真 |

**单例率仍需比较**——此前异常来自 harvester bug，不是「不必比」。修复后应用同一 k/w、同一 hg38，对比 native / xseal-min / xseal-sync 的 `singleton_pct` 与 `distinct_minimizers`。

---

## 5. OOM 现象（多次测试后）— 修订理解

**观测**：第一次测试正常；连续多组后出现 OOM、swap 打满。

**更符合的解释**（非单次建索引必然 OOM）：

1. **单次进程峰值仍高**（~14–20 GiB RSS）：`mi->S` + 索引桶 + 批处理 `m_vec`，但第一次可完成。
2. **连续 6 配置 × 多 trial**：每组全基因组建索引，内核 **page cache + swap** 累积；swap 用满后后续进程即使用户态已释放，仍极慢或 OOM-killer。
3. **`/dev/shm` 拷贝 hg38**（~3 GiB）在 ram 用例中常驻，占用可回收内存池。
4. **索引临时文件**（若未 `rm`）与 **glibc arena** 碎片次要。

**测试规程修订**（不写进核心代码即可缓解）：

- 每组测试 **独立进程、独立 shell 会话**；组间 `sync; sleep 30`；观察 `free -h` / `swapon --show`。
- 连跑不超过 **2–3 组**；或每组后 **重启 benchmark 脚本**（新进程树）。
- 内存紧张时：**不要** `cp` 到 `/dev/shm`；disk 用 `mmap` + 可选 `drop_caches` 即可。
- 记录 `RESULT|*_peak_rss_gib`；若 swap>0 且可用内存<4GiB，**停止后续组**。

---

## 6. 修订后的实施计划（待你确认后编码）

### 阶段 0：验证 harvester 修复（你已具备代码）

- [ ] 跑 `revcomp_verify`：broken vs fixed SIMD 对拍
- [ ] 单组 `xseal_min_disk`：`distinct` / `singleton_pct` 应相对修复前明显改善

### 阶段 A：指标与公平性（小改 `index.cpp`）

- [ ] 拆分 `RESULT`：`io_s` / `parse_s` / `compute_s` / `index_insert_s` / `frontend_s` / `e2e_s`
- [ ] 索引键改用 minimap2 `hash64(kmer, mask)`（与 `sketch.c` 一致），便于 **单例率与 distinct 和 native 对比**
- [ ] 保留 `seed_hits`、`*_singleton_pct`、`*_distinct_minimizers`

### 阶段 B：集成性能（消除双重 FASTA + 锁）

- [ ] 从 `mi->S` 分块编码扫描，取消第二次 `process_file`
- [ ] `mm_idx_add` 按 bucket 分片锁或 thread-local 桶缓冲

### 阶段 C：论文口径复现（可选）

- [ ] `XSEAL_BENCH_PAPER=1`：仅计 fused mmap 扫描（无 fill、无 add），复现 ~2.13s / ~14×

### 阶段 D：Benchmark 脚本

- [ ] `run_hg38_compare.sh`：组间冷却、可选跳过 shm、默认 **2 组/配置**、OOM 时 abort
- [ ] 输出 `summary` + 修复前后对比表（singleton、distinct、compute_s）

---

## 7. 建议的验证顺序（最小风险）

1. **只跑 1 次** `xseal_min_disk`（修复后 + 可选 hash64）→ 确认 singleton/distinct 正常、无 swap。
2. 再跑 **native_disk + xseal_min_disk** → 对比单例率、前端分段计时。
3. 确认无 OOM 后，再跑 syncmer / ram 组；**不要一次脚本跑满 6 组**。

---

## 8. 参考：minimap2 原生互补（对照）

`sketch.c`：`hash64(kmer[z], mask)`，`kmer[z]` 为滚动得到的 **canonical 整数 k-mer**（比较 `kmer[0]` 与 `kmer[1]` 后取较小），**无** Murmur3，**无**双重互补。

集成目标：T4 输出的索引键应与 `sketch.c` 的 `m.x >> 8` 语义一致。
