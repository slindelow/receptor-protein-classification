#!/usr/bin/env python
"""Freeze DAVIS (if needed), train, evaluate, export metrics."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from vs.pipeline import run_all  # noqa: E402


def main() -> None:
    results = run_all()
    print("=== metrics summary ===")
    for split, modes in results["splits"].items():
        for mode, m in modes.items():
            print(
                f"{split:14s} {mode:12s}  rho={m['spearman_rho']:.4f}  "
                f"EF1={m['ef_at_1pct']:.4f}  EF5={m['ef_at_5pct']:.4f}  "
                f"AUROC={m['auroc']:.4f}  n_test={int(m['n_test'])}"
            )
    print(f"metrics.md  -> {ROOT / 'artifacts' / 'metrics.md'}")
    print(f"metrics.json-> {ROOT / 'artifacts' / 'metrics.json'}")


if __name__ == "__main__":
    main()
