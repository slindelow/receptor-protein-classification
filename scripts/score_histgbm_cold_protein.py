#!/usr/bin/env python
"""Predict frozen v1.1 HistGBM models on the frozen DAVIS cold-protein test.

Loads artifacts/models_v1.1 joblibs. Does not refit. Does not rewrite
metrics_v1.1.json. Writes artifacts/scores_histgbm_cold_protein.csv only
when recomputed Spearman and EF@1% match the published point estimates
within 1e-4. Otherwise writes the mismatch and exits nonzero.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from vs.config import load_config  # noqa: E402
from vs.features import build_pair_matrix  # noqa: E402
from vs.metrics_lib import compute_metrics  # noqa: E402
from vs.model import load_artifact  # noqa: E402

TOL = 1e-4
CFG_PATH = ROOT / "configs" / "v1.1_protein.yaml"
PAIRS = ROOT / "data" / "raw" / "davis_pairs.csv"
SPLIT = ROOT / "data" / "splits" / "cold_protein_split.csv"
MODELS = ROOT / "artifacts" / "models_v1.1"
LIGAND_MODEL = MODELS / "hgb_cold_protein_ligand_only.joblib"
AAC_MODEL = MODELS / "hgb_cold_protein_full_aac_dpc.joblib"
METRICS = ROOT / "artifacts" / "metrics_v1.1.json"
OUT_CSV = ROOT / "artifacts" / "scores_histgbm_cold_protein.csv"
OUT_CHECK = ROOT / "artifacts" / "scores_histgbm_cold_protein_check.json"


def load_cold_test() -> pd.DataFrame:
    pairs = pd.read_csv(PAIRS)
    sid = pd.read_csv(SPLIT)
    keep = ["pair_id", "smiles", "sequence", "pkd", "drug_id", "target_id", "kd_nm"]
    keep = [c for c in keep if c in pairs.columns]
    merged = sid[["pair_id", "split"]].merge(pairs[keep], on="pair_id", how="inner")
    if len(merged) != len(sid):
        raise SystemExit(f"split join mismatch {len(merged)} vs {len(sid)}")
    test = merged[merged["split"] == "test"].reset_index(drop=True)
    if len(test) == 0:
        raise SystemExit("cold-protein test split is empty")
    return test


def predict(model_path: Path, X: np.ndarray) -> np.ndarray:
    art = load_artifact(model_path)
    model = art["model"] if isinstance(art, dict) else art
    pred = np.asarray(model.predict(X), dtype=float)
    if pred.shape != (len(X),):
        raise SystemExit(f"{model_path.name} predict shape {pred.shape} != {(len(X),)}")
    if not np.isfinite(pred).all():
        raise SystemExit(f"{model_path.name} produced non-finite predictions")
    return pred


def main() -> None:
    if not LIGAND_MODEL.exists() or not AAC_MODEL.exists():
        raise SystemExit("v1.1 HistGBM joblibs missing; refusing to refit in this script")
    cfg = load_config(CFG_PATH)
    radius = int(cfg["fingerprint"]["radius"])
    n_bits = int(cfg["fingerprint"]["n_bits"])
    thr = float(cfg["binder"]["pkd_threshold"])
    ef_fracs = list(cfg["eval"]["ef_fractions"])
    published = json.loads(METRICS.read_text())
    cold = published["splits"]["cold_protein"]
    pub_l = cold["ligand_only"]
    pub_p = cold["protein"]["aac_dpc"]

    test = load_cold_test()
    X_l, m_l = build_pair_matrix(
        test["smiles"], None, radius, n_bits, include_protein=False
    )
    X_p, m_p = build_pair_matrix(
        test["smiles"],
        test["sequence"],
        radius,
        n_bits,
        include_protein=True,
        protein_type="aac_dpc",
    )
    if not (m_l.all() and m_p.all()):
        raise SystemExit(
            f"invalid SMILES on cold-protein test: ligand mask {int((~m_l).sum())}, "
            f"protein mask {int((~m_p).sum())}. Published n_test is the full split, so scores were not written."
        )
    y = test["pkd"].to_numpy(dtype=float)
    pred_l = predict(LIGAND_MODEL, X_l)
    pred_p = predict(AAC_MODEL, X_p)
    met_l = compute_metrics(y, pred_l, thr, ef_fracs)
    met_p = compute_metrics(y, pred_p, thr, ef_fracs)

    checks = {
        "ligand_only_spearman": {
            "recomputed": met_l["spearman_rho"],
            "published": pub_l["spearman_rho"],
        },
        "aac_dpc_spearman": {
            "recomputed": met_p["spearman_rho"],
            "published": pub_p["spearman_rho"],
        },
        "ligand_only_ef_at_1pct": {
            "recomputed": met_l["ef_at_1pct"],
            "published": pub_l["ef_at_1pct"],
        },
        "aac_dpc_ef_at_1pct": {
            "recomputed": met_p["ef_at_1pct"],
            "published": pub_p["ef_at_1pct"],
        },
    }
    for name, rec in checks.items():
        rec["abs_diff"] = abs(float(rec["recomputed"]) - float(rec["published"]))
        rec["match"] = bool(rec["abs_diff"] <= TOL)
    matched = all(rec["match"] for rec in checks.values())
    check_doc = {
        "models": "loaded",
        "refit": False,
        "ligand_model": str(LIGAND_MODEL.relative_to(ROOT)),
        "aac_dpc_model": str(AAC_MODEL.relative_to(ROOT)),
        "split": str(SPLIT.relative_to(ROOT)),
        "n_pairs": int(len(test)),
        "n_kinases": int(test["target_id"].nunique()),
        "n_ligands": int(test["drug_id"].nunique()),
        "tolerance": TOL,
        "matched": matched,
        "checks": checks,
        "published_metrics_overwritten": False,
        "score_csv": str(OUT_CSV.relative_to(ROOT)) if matched else None,
    }
    OUT_CHECK.write_text(json.dumps(check_doc, indent=2) + "\n")
    print(json.dumps(check_doc, indent=2))
    if not matched:
        if OUT_CSV.exists():
            OUT_CSV.unlink()
        raise SystemExit(
            "Recomputed Spearman or EF@1% does not match metrics_v1.1.json within 1e-4. "
            "Published metrics were not overwritten. Score CSV was not written."
        )

    out = pd.DataFrame(
        {
            "pair_id": test["pair_id"].astype(str),
            "target_id": test["target_id"].astype(str),
            "drug_id": test["drug_id"].astype(str),
            "smiles": test["smiles"].astype(str),
            "pkd": y,
            "pred_ligand_only": pred_l,
            "pred_aac_dpc": pred_p,
        }
    )
    out.to_csv(OUT_CSV, index=False)
    print(f"wrote {OUT_CSV}")


if __name__ == "__main__":
    main()
