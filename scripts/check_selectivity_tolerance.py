#!/usr/bin/env python
"""Compare v1.4 selectivity metrics (not the splits schema)."""
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
        elif isinstance(v, bool):
            continue
        elif isinstance(v, (int, float)):
            out[key] = float(v)
    return out


def main() -> int:
    ref_path = Path(sys.argv[1])
    new_path = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "artifacts" / "metrics_selectivity.json"
    tol = 1e-4
    if len(sys.argv) > 3:
        tol = float(sys.argv[3])
    ref = json.loads(ref_path.read_text())
    new = json.loads(new_path.read_text())
    # Drop per-ligand table so we compare summary metrics only
    ref.pop("ligand_table_aac_dpc", None)
    new.pop("ligand_table_aac_dpc", None)
    rf, nf = flatten(ref), flatten(new)
    want = ("spearman_delta", "ef_at_", "mean_abs_delta", "window_sequence_identity", "aac_dpc_cosine", "n_selective", "n_dual")
    keys = [k for k in sorted(set(rf) & set(nf)) if any(w in k for w in want)]
    failed = []
    for k in keys:
        a, b = rf[k], nf[k]
        if abs(a - b) > tol and not (a != a and b != b):  # NaN==NaN ok
            failed.append((k, a, b, abs(a - b)))
    if failed:
        print("TOLERANCE FAIL")
        for row in failed:
            print(f"  {row[0]}: ref={row[1]} new={row[2]} absdiff={row[3]}")
        return 1
    print(f"TOLERANCE OK ({len(keys)} metrics, tol={tol})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
