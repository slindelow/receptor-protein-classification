from __future__ import annotations

import argparse
from pathlib import Path

from .config import ROOT, load_config, resolve
from .pipeline import screen_library


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(
        prog="vs-screen",
        description="Screen a SMILES library against a target sequence with the frozen HistGBM model.",
    )
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--target-seq", type=str, help="Amino-acid sequence of the target protein")
    g.add_argument("--target-seq-file", type=Path, help="File containing the amino-acid sequence")
    p.add_argument("--uniprot", type=str, default=None, help="Optional label only (no network fetch in MVP)")
    p.add_argument("--library", type=Path, required=True, help="CSV with smiles (+ optional mol_id)")
    p.add_argument("--out", type=Path, required=True, help="Output ranked CSV path")
    p.add_argument("--model", type=Path, default=None, help="Override model joblib path")
    p.add_argument("--config", type=Path, default=None)
    args = p.parse_args(argv)

    if args.target_seq_file:
        seq = Path(args.target_seq_file).read_text().strip().replace("\n", "").replace(" ", "")
    else:
        seq = args.target_seq.strip().replace("\n", "").replace(" ", "")

    ranked = screen_library(
        target_seq=seq,
        library_csv=args.library,
        model_path=args.model,
        cfg_path=args.config,
    )
    if args.uniprot:
        ranked.insert(1, "uniprot_label", args.uniprot)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    ranked.to_csv(args.out, index=False)
    print(f"Wrote {len(ranked)} rows -> {args.out}")


if __name__ == "__main__":
    main()
