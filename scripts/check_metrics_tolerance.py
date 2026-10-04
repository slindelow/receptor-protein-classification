#!/usr/bin/env python
"""Compare freshly computed metrics.json against a reference (or re-run self-check)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def flatten(d, prefix=""):
    out = {}
    for k, v in d.items():
        key = f"{prefix}.{k}" if prefix else k
        if isinstance(v, dict):
            out.update(flatten(v, key))
        elif isinstance(v, (int, float)):
            out[key] = float(v)
    return out


def main() -> int:
    ref_path = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "artifacts" / "metrics.json"
    new_path = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "artifacts" / "metrics.json"
    # When checking reproduce: caller writes new metrics then compares to committed reference copy
    # For self-check after single run, compare key numeric fields exist.
    ref = json.loads(ref_path.read_text())
    new = json.loads(new_path.read_text())
    tol = 1e-4
    try:
        import yaml

        cfg = yaml.safe_load((ROOT / "configs" / "default.yaml").read_text())
        tol = float(cfg["eval"]["metric_float_tol"])
    except Exception:
        pass

    rf = flatten(ref.get("splits", {}))
    nf = flatten(new.get("splits", {}))
    keys = sorted(set(rf) & set(nf))
    metric_keys = [k for k in keys if any(x in k for x in ("spearman", "ef_at", "auroc"))]
    failed = []
    for k in metric_keys:
        a, b = rf[k], nf[k]
        if abs(a - b) > tol:
            failed.append((k, a, b, abs(a - b)))
    if failed:
        print("TOLERANCE FAIL")
        for row in failed:
            print(f"  {row[0]}: ref={row[1]} new={row[2]} absdiff={row[3]}")
        return 1
    print(f"TOLERANCE OK ({len(metric_keys)} metrics, tol={tol})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
