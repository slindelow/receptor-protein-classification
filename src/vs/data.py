from __future__ import annotations

import hashlib
import json
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np
import pandas as pd

DAVIS_URL = "https://dataverse.harvard.edu/api/access/datafile/5219748"
DAVIS_TDC_FILE_ID = 5219748
TZ = ZoneInfo("America/Toronto")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def kd_nm_to_pkd(kd_nm: np.ndarray | pd.Series) -> np.ndarray:
    """Convert Kd in nM to pKd = 9 - log10(Kd_nM)."""
    kd = np.asarray(kd_nm, dtype=float)
    kd = np.clip(kd, 1e-10, None)
    return 9.0 - np.log10(kd)


def load_davis_tab(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    # TDC DAVIS columns: ID1, X1, ID2, X2, Y  (Y = Kd nM)
    rename = {
        "ID1": "drug_id",
        "X1": "smiles",
        "ID2": "target_id",
        "X2": "sequence",
        "Y": "kd_nm",
    }
    df = df.rename(columns=rename)
    df["smiles"] = df["smiles"].astype(str).str.strip().str.strip('"')
    df["sequence"] = df["sequence"].astype(str).str.strip().str.strip('"')
    df["target_id"] = df["target_id"].astype(str).str.strip().str.strip('"')
    df["drug_id"] = df["drug_id"].astype(str)
    df["kd_nm"] = pd.to_numeric(df["kd_nm"], errors="coerce")
    df = df.dropna(subset=["smiles", "sequence", "kd_nm", "target_id"]).copy()
    df["pkd"] = kd_nm_to_pkd(df["kd_nm"])
    df["pair_id"] = (
        df["drug_id"].astype(str) + "::" + df["target_id"].astype(str) + "::" + df.index.astype(str)
    )
    df = df.reset_index(drop=True)
    return df


def freeze_davis(raw_dir: Path, manifest_path: Path, source_url: str = DAVIS_URL) -> dict:
    raw_dir.mkdir(parents=True, exist_ok=True)
    tab_path = raw_dir / "davis_tdc.tab"
    if not tab_path.exists():
        import urllib.request

        urllib.request.urlretrieve(source_url, tab_path)

    df = load_davis_tab(tab_path)
    frozen_csv = raw_dir / "davis_pairs.csv"
    df.to_csv(frozen_csv, index=False)

    download_date = datetime.now(TZ).strftime("%Y-%m-%d")
    # Prefer mtime of tab if already present from earlier download today
    try:
        mtime = datetime.fromtimestamp(tab_path.stat().st_mtime, TZ).strftime("%Y-%m-%d")
        download_date = mtime
    except Exception:
        pass

    files = {
        "davis_tdc.tab": {
            "path": str(tab_path.relative_to(manifest_path.parent.parent))
            if False
            else "data/raw/davis_tdc.tab",
            "sha256": sha256_file(tab_path),
            "bytes": tab_path.stat().st_size,
            "note": "Raw TDC DAVIS tab from Harvard Dataverse (Kd nM)",
        },
        "davis_pairs.csv": {
            "path": "data/raw/davis_pairs.csv",
            "sha256": sha256_file(frozen_csv),
            "bytes": frozen_csv.stat().st_size,
            "note": "Sanitized pairs with pKd = 9 - log10(Kd_nM)",
        },
    }
    manifest = {
        "freeze_id": "sofia-vs-v1-2026-10",
        "dataset": "DAVIS",
        "source": {
            "name": "Therapeutics Data Commons (TDC) multi_pred.DTI DAVIS",
            "url": source_url,
            "dataverse_file_id": DAVIS_TDC_FILE_ID,
            "download_date_america_toronto": download_date,
            "citation": "Huang et al., Therapeutics Data Commons, NeurIPS 2021; Davis et al., Nat Biotechnol 2011",
        },
        "transform": {
            "pkd": "pkd = 9 - log10(Kd_nM); Kd reported in nM in TDC file",
            "n_pairs": int(len(df)),
            "n_unique_ligands": int(df["smiles"].nunique()),
            "n_unique_targets": int(df["target_id"].nunique()),
        },
        "files": files,
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)
        f.write("\n")
    return manifest
