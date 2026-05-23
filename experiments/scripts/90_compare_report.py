#!/usr/bin/env python3

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

"""Compare experiments/results/metrics_*.json to expected_metrics.yaml → results/report.md."""
from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path
from typing import Any


def load_expected_yaml(path: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    text = path.read_text(encoding="utf-8", errors="replace")
    tol_default: dict[str, float] = {}
    m_warn = re.search(r"warn_rel_pct:\s*([\d.]+)", text)
    m_fail = re.search(r"fail_rel_pct:\s*([\d.]+)", text)
    if m_warn:
        tol_default["warn_rel_pct"] = float(m_warn.group(1))
    if m_fail:
        tol_default["fail_rel_pct"] = float(m_fail.group(1))

    tol_loose = dict(tol_default)
    mloose = re.search(
        r"tolerances_loose:\s*\n\s*warn_rel_pct:\s*([\d.]+)\s*\n\s*fail_rel_pct:\s*([\d.]+)",
        text,
    )
    if mloose:
        tol_loose["warn_rel_pct"] = float(mloose.group(1))
        tol_loose["fail_rel_pct"] = float(mloose.group(2))

    tol_by_group: dict[str, dict[str, float]] = {"default": tol_default, "loose": tol_loose}

    metrics: list[dict[str, Any]] = []
    for block in re.split(r"\n\s*-\s+id:\s*", text):
        if "expected:" not in block:
            continue
        mid = re.match(r"(\S+)\s*\n", block)
        if not mid:
            continue
        mid_s = mid.group(1).strip().strip('"')
        ex = re.search(r"expected:\s*([\d.]+)", block)
        unit_m = re.search(r"unit:\s*(\S+)", block)
        cref_m = re.search(r'cref:\s*"([^"]*)"', block)
        desc_m = re.search(r"description:\s*\"([^\"]*)\"", block)
        tg_m = re.search(r"tolerance_group:\s*(\S+)", block)
        skip_m = re.search(r"skip_compare:\s*true", block)
        reason_m = re.search(r'skip_reason:\s*"([^"]*)"', block)
        if not ex:
            continue
        metrics.append(
            {
                "id": mid_s,
                "expected": float(ex.group(1)),
                "unit": unit_m.group(1) if unit_m else "",
                "cref": cref_m.group(1) if cref_m else "",
                "description": desc_m.group(1) if desc_m else "",
                "tolerance_group": tg_m.group(1) if tg_m else "default",
                "skip_compare": bool(skip_m),
                "skip_reason": reason_m.group(1) if reason_m else "",
            }
        )
    return tol_by_group, metrics


def read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {"error": "invalid_json", "path": str(path)}


def merge_metric_jsons(results: Path) -> dict[str, Any]:
    merged: dict[str, Any] = {}
    for p in sorted(results.glob("metrics_*.json")):
        data = read_json(p)
        for k, v in data.items():
            if k in ("error", "pure_scan_by_threads"):
                continue
            if isinstance(v, dict) and "mean" in v:
                merged[k] = v
    return merged


def pick_measured(all_metrics: dict[str, Any], metric_id: str) -> dict[str, Any] | None:
    if metric_id in all_metrics and isinstance(all_metrics[metric_id], dict):
        return all_metrics[metric_id]  # type: ignore[return-value]
    return None


def rel_pct(expected: float, mean: float) -> float:
    if expected == 0:
        return math.nan
    return abs(mean - expected) / expected * 100.0


def cpu_model() -> str:
    try:
        with open("/proc/cpuinfo", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("model name"):
                    return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return "unknown"


def main() -> int:
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent.parent
    results = root / "results"
    expected_path = root / "expected_metrics.yaml"
    out_md = results / "report.md"

    tol_by_group, expected_list = load_expected_yaml(expected_path)
    default_tol = tol_by_group.get("default", {"warn_rel_pct": 8.0, "fail_rel_pct": 15.0})
    loose_tol = tol_by_group.get("loose", default_tol)

    merged = merge_metric_jsons(results)

    lines: list[str] = []
    lines.append("# Paper reproduction report\n")
    lines.append(f"- Host CPU (from `/proc/cpuinfo`): **{cpu_model()}**")
    lines.append(
        "- 若与论文 Ryzen 5 5600X 不同，下列与 `expected_metrics.yaml` 的相对误差 **仅作回归/对照**，不等同于复现论文绝对数值。\n"
    )
    lines.append(
        "- `skip_compare: true` 的项只展示实测（若有），**不与印刷值算相对误差**（见 `skip_reason`）。\n"
    )
    lines.append("## Per-metric comparison\n")
    lines.append("| id | cref | expected | measured mean ± stdev (n) | rel err % | status |")
    lines.append("|----|------|----------|---------------------------|-----------|--------|")

    any_fail = False
    for m in expected_list:
        mid = m["id"]
        exp = float(m["expected"])
        unit = m.get("unit", "")
        cref = m.get("cref", "")
        tg = m.get("tolerance_group", "default")
        tol = loose_tol if tg == "loose" else default_tol
        warn_pct = float(tol.get("warn_rel_pct", 8.0))
        fail_pct = float(tol.get("fail_rel_pct", 15.0))

        if m.get("skip_compare"):
            obs = pick_measured(merged, mid)
            reason = m.get("skip_reason", "")
            if obs is None or obs.get("mean") is None:
                lines.append(
                    f"| `{mid}` | {cref} | {exp} {unit} | *missing* | — | ⚠ skip ({reason}) |"
                )
            else:
                mean = float(obs["mean"])
                stdev = float(obs.get("stdev", 0.0))
                n = int(obs.get("n", 1))
                lines.append(
                    f"| `{mid}` | {cref} | {exp} {unit} | {mean:.4f} ± {stdev:.4f} ({n}) | — | ℹ️ not compared |"
                )
            continue

        obs = pick_measured(merged, mid)
        if obs is None or obs.get("mean") is None:
            lines.append(f"| `{mid}` | {cref} | {exp} {unit} | *missing* | — | ⚠ skip |")
            continue
        mean = float(obs["mean"])
        stdev = float(obs.get("stdev", 0.0))
        n = int(obs.get("n", 1))
        rp = rel_pct(exp, mean)
        if math.isnan(rp):
            status = "?"
        elif rp > fail_pct:
            status = "🔴 fail"
            any_fail = True
        elif rp > warn_pct:
            status = "🟡 warn"
        else:
            status = "🟢 ok"
        lines.append(
            f"| `{mid}` | {cref} | {exp} {unit} | {mean:.4f} ± {stdev:.4f} ({n}) | {rp:.2f} | {status} |"
        )

    lines.append("\n## Raw artifacts\n")
    for p in sorted(results.iterdir()):
        if p.is_file() and not p.name.startswith("."):
            lines.append(f"- `{p.name}` ({p.stat().st_size} bytes)")

    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {out_md}")
    return 1 if any_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
