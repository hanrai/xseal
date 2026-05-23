# hg38 集成 Benchmark 报告（2026-05-18）

## 测试环境

| 项 | 值 |
|----|-----|
| 数据 | `/dev/shm/hg38_xseal_bench.fa`（约 3.1 Gbp，内存驻留 mmap） |
| 参数 | k=21, w=11；syncmer s=11（k−w+1） |
| 线程 | `-t 12`，`XSEAL_WORKERS=11`（11 worker + 1 主线程） |
| 索引键 | minimap2 `hash64(canonical_kmer)`（与 `sketch.c` 一致） |
| 二进制 | `minimap2_xseal` / `minimap2_xseal_syncmer`（commit `aa46679` 及之前修复） |
| 原始日志 | `integration/minimap2_xseal_hardened/bench_hg38_final/*.log` |

### 指标单位与含义（必读）

下文表格中的 **`compute_wall`、`index_insert`、`frontend`、`e2e` 等列，单位均为秒（s）**，由 `std::chrono::steady_clock` 墙钟计时累加后输出为 `RESULT|*_s|mean|…`。  
吞吐列 `*_gbp_s` 单位为 **Gbp/s（十亿碱基对每秒）**，分子为参考序列总碱基数（hg38 约 **3.209×10⁹ bp**），分母为对应阶段的秒数。

#### 四个核心指标（统一墙钟口径）

| 指标名 | 单位 | 时钟类型 | 准确含义 |
|--------|------|----------|----------|
| **`compute_wall_s`** | 秒 | **墙钟** | XSeal：`stop()` 等待全部 worker 结束（扫描+收割+**并行**入桶都在此窗口内）。Native：`mm_sketch` 阶段累计。 |
| **`index_insert_s`** | 秒 | **墙钟** | **与 `index_insert_wall_s` 相同**：从**第一次**到最后一次入桶回调的**时间跨度**（不是各线程相加）。 |
| **`index_insert_accum_s`** | 秒 | **线程时间之和** | 每个 worker 回调 `(建 m_vec + 入桶)` 耗时**累加**；含锁等待，可 ≫ 墙钟（仅诊断锁竞争）。 |
| **`frontend_s`** | 秒 | **墙钟（分项相加）** | 见下方公式；**禁止**再把 `index_insert_accum_s` 加进 `frontend_s`。 |
| **`e2e_s`** | 秒 | **墙钟（分项相加）** | `io_s + frontend_s + backend_s`。 |

**为何会出现「112 秒入桶 vs 9 秒 native」的假象（已修复）**

旧实现把 **每个线程** 回调耗时（含排队等锁）**全部相加**，得到 `index_insert_accum_s ≈ 112s`，同时又把 **`compute_wall_s`（~11s 墙钟）** 和该累加值 **一起加进 `frontend_s`**。  
但入桶发生在 worker 内部，与 `compute_wall_s` **同一墙钟窗口重叠** → **重复计时**，导致 `frontend_s` 虚高到 124s。

**当前代数关系（`index.cpp` → `emit_index_timing`）：**

```
# XSeal（统一口径）
frontend_s = file_map_s + parse_s + compute_wall_s     # 入桶已含在 compute_wall 墙钟内
index_insert_s = index_insert_wall_s                   # 首末回调跨度，仅作对照，不加入 frontend

# Native（pipeline 分阶段，阶段和可大于墙钟，因 kt_pipeline 并行）
frontend_s = compute_wall_s + index_insert_accum_s    # sketch + dispatch

e2e_s = io_s + frontend_s + backend_s
```

**对比 native 9s vs XSeal 入桶（修复后，`bench_quick` /shm）：**

| | `index_insert_wall_s` | `index_insert_accum_s` | 说明 |
|--|----------------------|------------------------|------|
| native | ~9.2s | ~9.2s | 单阶段累加≈墙钟 |
| xseal minimizer（分桶锁） | **~14.0s** | ~134s（诊断） | 墙钟与 `compute_wall` ~13s 同量级 |
| xseal minimizer（**TLS 按桶缓冲**） | **~8.8s** | ~78s（诊断） | 每 task `push`→`flush`；墙钟与 `compute_wall` ~8.0s 同量级 |

#### `compute_wall_s` — 各路径具体测什么

**XSeal（minimizer / syncmer）**

- **起止**：`process_file()` 内，`ingest_buffer()` 结束之后调用 `stop()` → `stop()` 返回为止（`xseal_frontend_distributor.hpp`）。
- **包含**：11 个 worker 线程上的 **2-bit 编码 + 窗口扫描（minimizer/syncmer）+ canonical k-mer 收割**；主线程在 `stop()` 中 **等待所有 worker `join`** 的时间。
- **不包含**：`open`/`mmap`（计入 `file_map_s`）、FASTA 解析与任务切片（`parse_s`）、`fill_sequences_from_fasta`（`io_s`）、`mm_idx_post`（`backend_s`）。
- **注意**：worker 在扫描**之后同一任务内**会调用入桶回调，故 **入桶墙钟与 `compute_wall_s` 重叠**；报告 `index_insert_wall_s` 仅供分拆观察，**不计入 `frontend_s`**。
- **与论文 “sketching ~2.13s” 最接近的口径**：应用 **`compute_wall_s`**（或 distributor 打印的 `End-to-End Time` ≈ `parse_s + compute_wall_s`，仍不含 `index_insert`）。

**Native（`XSEAL_USE_NATIVE=1`）**

- **等于** minimap2 `kt_pipeline` **step 1（`mm_sketch`）** 的墙钟累计（`g_native_sketch_ns`）。
- **不包含**：step 0 读序列（`io_s`）、step 2 把 sketch 结果 `mm_idx_add` 入桶（计入 **`index_insert_s`**，即 `g_native_dispatch_ns`）。
- 对比时：**native `compute_wall_s` ≈ XSeal `compute_wall_s`**（都是“生成 minimizer/syncmer 命中”的核心计算，而非整个 native pipeline）。

#### `index_insert_s` / `index_insert_accum_s`

**XSeal**

- **`index_insert_wall_s`**：`min(回调开始)` → `max(回调结束)` 的墙钟跨度。
- **`index_insert_accum_s`**：每次回调耗时之和（11 个 worker 并行时，总和可达 **10× 墙钟**）。
- **实现**：`xseal_index_insert_tls.hpp` — 每 worker **thread_local** 按桶 `vector<mm128_t>` 缓冲，每个 fused task 结束 `flush()` 时持**桶锁**批量 `mm_idx_add`；去掉全局 `mm_idx_lock`。

**Native**

- `index_insert_accum_s` = pipeline step 2（`mm_idx_add`）各线程耗时之和；无单独 wall 跨度时与 accum 相同。

#### `frontend_s` — 为何故意去掉 `io_s`

集成设计里参考序列被读 **两次**：

1. **`io_s`**：`fill_sequences_from_fasta` → 填充 `mi->seq` / `mi->S`（与 minimap2 原生索引所需序列存储一致）。
2. **第二次**：XSeal `process_file` → `file_map` + `parse` + `compute_wall` + `index_insert`。

论文对比 “minimap2 前端 vs XSeal 前端” 时，通常 **不把第一次读盘算进 XSeal 前端**（native 侧对应部分在 `io_s`）。若要比 “含读盘的完整建索引”，应使用 **`e2e_s`**。

**Native 的 `frontend_s`**：`file_map_s` 与 `parse_s` 恒为 0，故  
`frontend_s = compute_wall_s + index_insert_accum_s`（sketch + dispatch，pipeline 可分阶段并行）。

#### `e2e_s` — 与进程 `Real time` 的区别

| | `e2e_s` | `[M::main] Real time` |
|--|---------|------------------------|
| 定义 | `io_s + frontend_s + backend_s` 分项相加 | 操作系统看到的进程总墙钟 |
| 是否含 CLI/重复 `mm_idx_stat` | 否（仅建索引主路径） | **可能偏多**（日志末尾常再打一次 stat） |
| 用途 | 论文式分层、与 `RESULT` 行一致 | 运维/整机感受 |

#### 相关辅助指标

| 字段 | 单位 | 含义 |
|------|------|------|
| `io_s` | 秒 | 第一次 FASTA → `mi->S`（XSeal：`fill_sequences_from_fasta`；Native：pipeline step 0） |
| `file_map_s` | 秒 | 仅 XSeal：`open` + `mmap` 或 preload 读入 |
| `parse_s` | 秒 | 仅 XSeal：`ingest_buffer`（解析 header/切片、推任务队列） |
| `worker_compute_s` | 秒 | 仅 XSeal：各 worker **任务内**耗时之和（**CPU 时间叠加**，可 ≫ `compute_wall_s`） |
| `backend_s` | 秒 | `mm_idx_post`：按桶排序、建 `idx` 哈希表、统计单例率等 |
| `compute_gbp_s` | Gbp/s | `总碱基数 / compute_wall_s` |
| `frontend_gbp_s` | Gbp/s | `总碱基数 / frontend_s` |
| `e2e_gbp_s` | Gbp/s | `总碱基数 / e2e_s` |

---

## 1. 墙钟时间（秒）

### 1.1 端到端与分层（统一口径，`bench_tls` smoke，/dev/shm，n=1）

| 配置 | io_s | parse_s | compute_wall_s | index_insert_wall_s | index_insert_accum_s | **frontend_s** | backend_s | **e2e_s** |
|------|------|---------|----------------|---------------------|----------------------|----------------|-----------|-----------|
| **native** | 9.87 | 0.00 | 29.95 | 9.13 | 9.13 | **39.08** | 3.04 | **51.99** |
| **xseal minimizer** | 9.27 | 0.86 | 8.02 | 8.81 | 78.07 | **8.88** | 2.76 | **20.90** |

> 分桶锁阶段见 `bench_quick`（frontend ~14s）；上表为 **TLS 按桶缓冲** 后单次试跑。完整三轮可重跑：`OUT_DIR=bench_hg38_unified … ./run_hg38_compare.sh`。

**解读**

- **论文式 fused 对比**：`compute_wall_s` native 30.0s vs XSeal 8.0s → **~3.7×**。
- **集成 frontend（墙钟）**：XSeal **8.9s** vs native **39.1s** → **~4.4×**（`frontend_s` 不再含重复的 `index_insert_accum`）。
- **e2e**：XSeal **20.9s** vs native **52.0s**（仍含双方 `io_s` ~9–10s 读盘；XSeal 第二次 mmap 路径 `file_map_s≈0`）。

### 1.2 吞吐（Gbp/s）

| 配置 | compute_gbp_s | frontend_gbp_s | e2e_gbp_s |
|------|---------------|----------------|-----------|
| native | 0.107 | 0.082 | 0.062 |
| xseal minimizer（TLS） | **0.400** | **0.361** | **0.154** |
| xseal syncmer | — | — | — |

### 1.3 并行度（仅 XSeal）

| 配置 | compute_wall_s | worker_compute_s | 并行效率† |
|------|----------------|------------------|-----------|
| xseal minimizer（TLS） | 8.02 | 94.16 | ~11.7× |
| xseal syncmer | — | — | — |

† `worker_compute_s / compute_wall_s`，接近 11 个 worker。

### 1.4 Fused 流水线内部计时（stderr，不含 index_insert）

| 配置 | 内部 E2E‡ | 内部 MiB/s | seed_hits |
|------|-----------|------------|-----------|
| xseal minimizer | 12.07 s | 241.0 | 512,205,311 |
| xseal syncmer | 13.42 s | 216.7 | 565,698,074 |

‡ `[XSeal HARDENED BENCHMARK] End-to-End Time` ≈ `parse_s + compute_wall_s`（mmap 后 fused 扫描）。

---

## 2. 索引质量（与 native 可比）

| 配置 | distinct | singleton% | avg_occ | avg_spacing | peak_rss |
|------|----------|------------|---------|-------------|----------|
| native | 381,630,162 | 93.71% | 1.335 | 6.297 | 13.9 GiB |
| xseal minimizer | 391,783,782 | **94.26%** | 1.307 | 6.266 | 19.9 GiB |
| xseal syncmer | 412,371,155 | 93.34% | 1.372 | 5.673 | 23.1 GiB |

---

## 3. 与论文数字对照

| 论文/口径 | 本次测量 | 说明 |
|-----------|----------|------|
| XSeal sketch ~2.13 s | `compute_wall_s` ≈ **11.2 s** | 全 hg38 + k21/w11 + 11 worker；非同一 microbench |
| 14.9× frontend | `compute_wall` native/xseal_min ≈ **29.95/8.02 ≈ 3.7×** | 论文多为纯 fused、无二次 `io_s` |
| 集成端到端 | native **52.0 s** vs xseal_min **20.9 s** | 修复计时 + 分桶锁 + TLS 后；`index_insert_wall` 与 compute 重叠，不再拖垮 `frontend_s` |

**论文可辩护口径**：用 `compute_wall_s` / `compute_gbp_s`（≈0.40 Gbp/s vs native 0.11）说明 **fused 扫描** 加速；集成 E2E 需披露 `io_s`（双方均有一次 FASTA→`mi->S`）。

---

## 4. 已修复问题（影响单例率与 distinct）

### 4.1 Harvester 双重互补（T4，`ac71e6d`）

- SIMD 主路径：`andnot` + LUT → 仅 LUT；LUT 与 syncmer 对齐。
- 影响：修复前 distinct ~19M、singleton ~49%（错误键塌缩）。

### 4.2 Minimizer shuffle 写入下标（`261ba83` / `aa46679`）

- **根因**：`shuffle_lut` 将有效 lane 压到 `temp_p[0..popcount-1]`，旧代码 `if (mask>>i) out[k++]=temp_p[i]` 误用 lane 下标 `i`。
- **现象**：同一 hash 多次入索引 → singleton ~52%（**非**窗口未去重）。
- **修复后**：singleton **94.26%**，与 native **93.71%** 一致。

---

## 5. 复现命令

```bash
cd integration/minimap2_xseal_hardened
make -j$(nproc)
OUT_DIR=bench_hg38_final COOLDOWN_SEC=30 TRIALS=1 \
  HG38_FASTA=/dev/shm/hg38_xseal_bench.fa \
  CASES="native_disk xseal_min_disk xseal_syn_disk" \
  ./run_hg38_compare.sh
```

机器可读 `RESULT|...` 行见各 `*.log`；汇总见 `bench_hg38_final/summary.txt`。

---

## 6. 结论摘要

1. **Fused 计算**：minimizer `compute_wall` 约为 native sketch 的 **~3.7×** 快（墙钟，`bench_tls`）。
2. **集成 E2E**：TLS + 统一计时后 XSeal **~21 s** vs native **~52 s**；旧版 ~112 s `index_insert` 为**线程时间累加且重复计入 frontend** 的统计假象。
3. **索引正确性**：修复 harvester + store 后，minimizer **singleton/distinct 与 native 同量级**。
4. **后续**：完整三轮 benchmark、syncmer TLS 验证、可选单次 FASTA / `XSEAL_BENCH_PAPER=1` 纯 fused 口径。
