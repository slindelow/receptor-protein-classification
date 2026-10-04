#!/usr/bin/env python
"""Build a ~400 SMILES demo library (all DAVIS ligands + ESOL pad) and sample ranked.csv."""
from __future__ import annotations

import hashlib
import json
import sys
import urllib.request
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np
import pandas as pd
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from vs.config import load_config, resolve  # noqa: E402
from vs.pipeline import screen_library  # noqa: E402

ESOL_URL = "https://raw.githubusercontent.com/deepchem/deepchem/master/datasets/delaney-processed.csv"
TZ = ZoneInfo("America/Toronto")


def _ok(s: str) -> bool:
    return Chem.MolFromSmiles(str(s)) is not None


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    cfg = load_config()
    raw_dir = resolve(cfg, "raw_dir")
    esol_path = raw_dir / "delaney_esol.csv"
    if not esol_path.exists():
        urllib.request.urlretrieve(ESOL_URL, esol_path)

    davis = pd.read_csv(raw_dir / "davis_pairs.csv")
    davis_lig = (
        davis.drop_duplicates("smiles")[["drug_id", "smiles"]]
        .rename(columns={"drug_id": "mol_id"})
        .assign(
            mol_id=lambda d: d["mol_id"].astype(str).map(lambda x: f"davis_{x}"),
            source="DAVIS",
        )
    )
    davis_lig = davis_lig[davis_lig["smiles"].map(_ok)]

    esol = pd.read_csv(esol_path)[["Compound ID", "smiles"]].dropna().drop_duplicates("smiles")
    esol = esol.rename(columns={"Compound ID": "mol_id"})
    esol["mol_id"] = esol["mol_id"].astype(str).map(lambda x: f"esol_{x}")
    esol["source"] = "ESOL_Delaney"
    esol = esol[esol["smiles"].map(_ok)]

    rng = np.random.default_rng(cfg["seed"])
    n_target = 400
    need = max(0, n_target - len(davis_lig))
    if need < len(esol):
        esol = esol.iloc[rng.choice(len(esol), size=need, replace=False)]
    lib = (
        pd.concat([davis_lig, esol], ignore_index=True)
        .drop_duplicates("smiles")
        .sample(frac=1.0, random_state=cfg["seed"])
        .reset_index(drop=True)
    )
    lib_path = resolve(cfg, "demo_library")
    lib_path.parent.mkdir(parents=True, exist_ok=True)
    lib[["mol_id", "smiles", "source"]].to_csv(lib_path, index=False)

    tid = davis["target_id"].value_counts().index[0]
    seq = davis.loc[davis["target_id"] == tid, "sequence"].iloc[0]
    (lib_path.parent / "demo_target.txt").write_text(f"target_id={tid}\nsequence={seq}\n")
    ranked = screen_library(seq, lib_path)
    out = resolve(cfg, "sample_ranked")
    ranked.to_csv(out, index=False)

    # refresh demo entries in MANIFEST if present
    man_path = ROOT / "data" / "MANIFEST.json"
    if man_path.exists():
        man = json.loads(man_path.read_text())
        man["files"]["delaney_esol.csv"] = {
            "path": "data/raw/delaney_esol.csv",
            "sha256": _sha256(esol_path),
            "bytes": esol_path.stat().st_size,
            "note": "Public Delaney/ESOL SMILES used only to pad demo screening library (not training labels)",
            "source_url": ESOL_URL,
            "download_date_america_toronto": datetime.now(TZ).strftime("%Y-%m-%d"),
        }
        man["files"]["demo_library.csv"] = {
            "path": "data/demo/library.csv",
            "sha256": _sha256(lib_path),
            "bytes": lib_path.stat().st_size,
            "note": f"Demo screen library n={len(lib)} (DAVIS ligands + ESOL pad)",
        }
        man_path.write_text(json.dumps(man, indent=2, sort_keys=True) + "\n")

    print(f"demo library n={len(lib)} -> {lib_path}")
    print(f"sample ranked -> {out} (target={tid})")


if __name__ == "__main__":
    main()
