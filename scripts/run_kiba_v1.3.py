#!/usr/bin/env python
"""v1.3 KIBA secondary freeze: same split protocol as DAVIS MVP, HistGBM ECFP+aac_dpc.

Downloads KIBA via TDC Dataverse if needed, freezes pairs, builds scaffold + cold_protein
splits (new CSVs under data/splits_kiba/), trains ligand-only vs aac_dpc, exports
artifacts/metrics_kiba.md + .json. Does not touch DAVIS splits.
"""
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from vs.config import load_config  # noqa: E402
from vs.data import kd_nm_to_pkd  # noqa: E402
from vs.features import build_pair_matrix, murcko_scaffold  # noqa: E402
from vs.metrics_lib import compute_metrics  # noqa: E402
from vs.model import save_artifact, train_model  # noqa: E402
from vs.splits import cold_protein_split, scaffold_split  # noqa: E402

TZ = ZoneInfo("America/Toronto")

# TDC KIBA file id (Harvard Dataverse) — verified via TDC docs / common mirror
KIBA_URL = "https://dataverse.harvard.edu/api/access/datafile/4156619"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_kiba_tab(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    rename = {"ID1": "drug_id", "X1": "smiles", "ID2": "target_id", "X2": "sequence", "Y": "kd_nm"}
    # KIBA Y is already a continuous score (not necessarily nM Kd). TDC documents Y as KIBA score.
    # We treat Y as affinity score directly as "pkd-like" label for ranking (document honestly).
    if "Y" in df.columns and "X1" in df.columns:
        df = df.rename(columns=rename)
    df["smiles"] = df["smiles"].astype(str).str.strip().str.strip('"')
    df["sequence"] = df["sequence"].astype(str).str.strip().str.strip('"')
    df["target_id"] = df["target_id"].astype(str).str.strip().str.strip('"')
    df["drug_id"] = df["drug_id"].astype(str)
    # KIBA score: higher = stronger. Use as regression target directly; store as pkd column for protocol reuse.
    df["kiba_score"] = pd.to_numeric(df["kd_nm"], errors="coerce")
    df = df.dropna(subset=["smiles", "sequence", "kiba_score", "target_id"]).copy()
    df["pkd"] = df["kiba_score"].astype(float)  # protocol reuse; NOT true pKd
    df["kd_nm"] = np.nan
    df["pair_id"] = (
        df["drug_id"].astype(str) + "::" + df["target_id"].astype(str) + "::" + df.index.astype(str)
    )
    return df.reset_index(drop=True)


def freeze_kiba(raw_dir: Path, manifest_path: Path) -> tuple[pd.DataFrame, dict]:
    raw_dir.mkdir(parents=True, exist_ok=True)
    pairs = raw_dir / "kiba_pairs.csv"
    # Prefer pre-assembled DeepDTA mirror (TDC Dataverse often 403 from this box)
    if pairs.exists() and manifest_path.exists():
        df = pd.read_csv(pairs)
        manifest = json.loads(manifest_path.read_text())
        print(f"Loaded existing KIBA freeze: {pairs} n={len(df)}")
        return df, manifest
    tab = raw_dir / "kiba_tdc.tab"
    if not tab.exists():
        print(f"Downloading KIBA from {KIBA_URL} ...")
        try:
            urllib.request.urlretrieve(KIBA_URL, tab)
        except Exception as e:
            raise SystemExit(
                f"KIBA TDC download failed ({e}). Assemble from DeepDTA mirror into data/raw/kiba_pairs.csv first."
            )
    df = load_kiba_tab(tab)
    df.to_csv(pairs, index=False)
    download_date = datetime.now(TZ).strftime("%Y-%m-%d")
    manifest = {
        "freeze_id": "sofia-vs-kiba-v1.3",
        "dataset": "KIBA",
        "source": {
            "name": "TDC multi_pred.DTI KIBA",
            "url": KIBA_URL,
            "download_date_america_toronto": download_date,
            "citation": "Tang et al., J Chem Inf Model 2014; Huang et al., TDC NeurIPS 2021",
            "label_note": "Y is KIBA score (not Kd nM). Stored as pkd column for protocol reuse only.",
        },
        "transform": {
            "label": "pkd := kiba_score (NOT 9-log10Kd); binder threshold set separately for KIBA",
            "n_pairs": int(len(df)),
            "n_unique_ligands": int(df["smiles"].nunique()),
            "n_unique_targets": int(df["target_id"].nunique()),
        },
        "files": {
            "kiba_tdc.tab": {"path": "data/raw/kiba_tdc.tab", "sha256": sha256_file(tab), "bytes": tab.stat().st_size},
            "kiba_pairs.csv": {"path": "data/raw/kiba_pairs.csv", "sha256": sha256_file(pairs), "bytes": pairs.stat().st_size},
        },
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return df, manifest


def _cached_pair_matrix(df, include_protein, protein_type, radius, n_bits, fp_cache, prot_cache):
    """Build X using caches keyed by smiles / sequence (KIBA has heavy ligand reuse)."""
    from vs.features import smiles_to_ecfp, sequence_to_protein_features, protein_feat_dim
    n = len(df)
    pdim = protein_feat_dim(protein_type) if include_protein else 0
    X = np.zeros((n, n_bits + pdim), dtype=np.float32)
    mask = np.ones(n, dtype=bool)
    smis = df["smiles"].astype(str).tolist()
    seqs = df["sequence"].astype(str).tolist() if include_protein else [None] * n
    for i, (smi, seq) in enumerate(zip(smis, seqs)):
        if smi not in fp_cache:
            fp_cache[smi] = smiles_to_ecfp(smi, radius=radius, n_bits=n_bits)
        fp = fp_cache[smi]
        if fp is None:
            mask[i] = False
            continue
        if include_protein:
            if seq not in prot_cache:
                prot_cache[seq] = sequence_to_protein_features(seq or "", protein_type)
            X[i] = np.concatenate([fp, prot_cache[seq]])
        else:
            X[i] = fp
    return X, mask


def fit_eval(df_split, cfg, include_protein, protein_type="aac_dpc", binder_thr=12.1, fp_cache=None, prot_cache=None):
    radius = cfg["fingerprint"]["radius"]
    n_bits = cfg["fingerprint"]["n_bits"]
    fp_cache = {} if fp_cache is None else fp_cache
    prot_cache = {} if prot_cache is None else prot_cache
    train = df_split[df_split["split"] == "train"]
    test = df_split[df_split["split"] == "test"]
    print(f"  featurize train n={len(train)} cache_fp={len(fp_cache)} ...", flush=True)
    X_tr, m_tr = _cached_pair_matrix(train, include_protein, protein_type, radius, n_bits, fp_cache, prot_cache)
    y_tr = train["pkd"].to_numpy()[m_tr]
    X_tr = X_tr[m_tr]
    print(f"  featurize test n={len(test)} ...", flush=True)
    X_te, m_te = _cached_pair_matrix(test, include_protein, protein_type, radius, n_bits, fp_cache, prot_cache)
    y_te = test["pkd"].to_numpy()[m_te]
    X_te = X_te[m_te]
    print(f"  train HistGBM X={X_tr.shape} ...", flush=True)
    model = train_model(X_tr, y_tr, cfg)
    pred = model.predict(X_te)
    metrics = compute_metrics(y_te, pred, binder_thr, cfg["eval"]["ef_fractions"])
    metrics["n_train"] = float(len(y_tr))
    metrics["n_test"] = float(len(y_te))
    metrics["protein_type"] = protein_type if include_protein else "none"
    return metrics, model


def main():
    cfg = load_config(ROOT / "configs" / "default.yaml")
    # KIBA binder: common literature uses KIBA score >= 12.1
    binder_thr = 12.1
    raw_dir = ROOT / "data" / "raw"
    splits_dir = ROOT / "data" / "splits_kiba"
    models_dir = ROOT / "artifacts" / "models_kiba"
    splits_dir.mkdir(parents=True, exist_ok=True)
    models_dir.mkdir(parents=True, exist_ok=True)

    df, manifest = freeze_kiba(raw_dir, ROOT / "data" / "MANIFEST_kiba.json")
    print(f"KIBA pairs={len(df)} ligands={df['smiles'].nunique()} targets={df['target_id'].nunique()}")

    # Build splits once (same fracs/seeds as DAVIS)
    scaf_path = splits_dir / "scaffold_split.csv"
    cold_path = splits_dir / "cold_protein_split.csv"
    if not scaf_path.exists():
        print("Building scaffold split...")
        scaf = scaffold_split(
            df,
            test_frac=cfg["splits"]["scaffold_test_frac"],
            val_frac=cfg["splits"]["val_frac"],
            seed=cfg["seed"],
        )
        scaf.to_csv(scaf_path, index=False)
    else:
        scaf = pd.read_csv(scaf_path)
        # rejoin labels
        scaf = scaf[["pair_id", "split"]].merge(df, on="pair_id", how="inner")

    if not cold_path.exists():
        print("Building cold_protein split...")
        cold = cold_protein_split(
            df,
            test_frac=cfg["splits"]["cold_protein_test_frac"],
            val_frac=cfg["splits"]["val_frac"],
            seed=cfg["seed"],
        )
        cold.to_csv(cold_path, index=False)
    else:
        cold = pd.read_csv(cold_path)
        cold = cold[["pair_id", "split"]].merge(df, on="pair_id", how="inner")

    # If freshly built, scaffold_split returns full df with split col
    if "pkd" not in scaf.columns:
        scaf = scaf.merge(df, on="pair_id", how="inner")
    if "pkd" not in cold.columns:
        cold = cold.merge(df, on="pair_id", how="inner")

    results = {
        "freeze_id": "sofia-vs-kiba-v1.3",
        "data_freeze_id": "sofia-vs-kiba-v1.3",
        "generated_at_america_toronto": datetime.now(TZ).isoformat(timespec="seconds"),
        "binder_threshold_kiba_score": binder_thr,
        "label_note": "pkd column holds KIBA score (not true pKd)",
        "splits": {},
        "models": {},
        "manifest_sha256": {k: v["sha256"] for k, v in manifest["files"].items()},
    }

    for split_name, dsplit in (("scaffold", scaf), ("cold_protein", cold)):
        results["splits"][split_name] = {}
        fp_cache, prot_cache = {}, {}
        for include, tag, ptype in (
            (False, "ligand_only", "none"),
            (True, "aac_dpc", "aac_dpc"),
        ):
            print(f"=== KIBA {split_name} {tag} ===", flush=True)
            metrics, model = fit_eval(
                dsplit, cfg, include,
                protein_type=ptype if include else "aac_dpc",
                binder_thr=binder_thr,
                fp_cache=fp_cache,
                prot_cache=prot_cache,
            )
            results["splits"][split_name][tag] = metrics
            art = models_dir / f"hgb_kiba_{split_name}_{tag}.joblib"
            save_artifact(art, model, {"split": split_name, "mode": tag, "freeze_id": "sofia-vs-kiba-v1.3", "metrics": metrics})
            results["models"][f"{split_name}_{tag}"] = str(art.relative_to(ROOT))
            print(f"  rho={metrics['spearman_rho']:.4f} EF1={metrics['ef_at_1pct']:.4f}", flush=True)

    js = ROOT / "artifacts" / "metrics_kiba.json"
    md = ROOT / "artifacts" / "metrics_kiba.md"
    js.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")

    def fmt(x):
        return f"{x:.4f}" if isinstance(x, float) else str(x)

    lines = [
        "# Virtual Screening KIBA secondary freeze (v1.3)",
        "",
        f"**Experiment freeze:** `sofia-vs-kiba-v1.3`  ",
        f"**Generated (America/Toronto):** {results['generated_at_america_toronto']}  ",
        f"**Binder threshold:** KIBA score >= {binder_thr}  ",
        "**Label:** KIBA score stored in `pkd` column for protocol reuse (NOT true pKd).  ",
        "**Model:** HistGBM ECFP +/- AAC+DPC (same hyperparams as DAVIS MVP).  ",
        "",
        "| Split | Model | Spearman rho | EF@1% | EF@5% | AUROC | n_test | n_binders |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for split in ("scaffold", "cold_protein"):
        for tag, label in (("ligand_only", "HistGBM ligand-only"), ("aac_dpc", "HistGBM+AAC+DPC")):
            m = results["splits"][split][tag]
            lines.append(
                f"| {split} | {label} | {fmt(m['spearman_rho'])} | {fmt(m['ef_at_1pct'])} | "
                f"{fmt(m['ef_at_5pct'])} | {fmt(m['auroc'])} | {int(m['n_test'])} | {int(m['n_binders'])} |"
            )
    cold = results["splits"]["cold_protein"]
    d = cold["aac_dpc"]["spearman_rho"] - cold["ligand_only"]["spearman_rho"]
    lines += [
        "",
        "## Success (cold_protein protein vs ligand-only)",
        f"- Delta Spearman = {fmt(d)}",
        f"- Pass Delta>=0.05 or EF lift: {'YES' if (d >= 0.05 or (cold['aac_dpc']['ef_at_1pct']-cold['ligand_only']['ef_at_1pct'])>0.5) else 'NO'}",
        "",
        "## Notes",
        "",
        "- DAVIS splits untouched. KIBA splits under `data/splits_kiba/`.",
        "- No wet-lab / NOVA claims.",
        "",
    ]
    md.write_text("\n".join(lines))
    print(f"Wrote {md}")


if __name__ == "__main__":
    main()
