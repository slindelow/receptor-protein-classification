#!/usr/bin/env python
"""One protein-cluster test of the cold-protein pooled Spearman lift.

Does not refit models. Uses saved Chemprop predictions only.
Writes artifacts/metrics_stat_review.json.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[1]
ART = ROOT / "artifacts"
SEED = 42
N_BOOT = 1000


def spearman_rho(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    if y_true.size < 2 or np.allclose(y_true, y_true[0]) or np.allclose(y_pred, y_pred[0]):
        return float("nan")
    ry = rankdata(y_true)
    rp = rankdata(y_pred)
    ry = ry - ry.mean()
    rp = rp - rp.mean()
    denom = float(np.sqrt(np.dot(ry, ry) * np.dot(rp, rp)))
    if denom == 0.0:
        return float("nan")
    return float(np.dot(ry, rp) / denom)


def percentile_ci(samples: np.ndarray) -> dict:
    finite = samples[np.isfinite(samples)]
    if finite.size == 0:
        return {"low": None, "high": None, "n_finite": 0, "includes_0": None}
    low, high = np.percentile(finite, [2.5, 97.5])
    return {
        "low": float(low),
        "high": float(high),
        "n_finite": int(finite.size),
        "includes_0": bool(low <= 0.0 <= high),
    }


def twosided_recentered_p(samples: np.ndarray, observed: float) -> dict:
    finite = samples[np.isfinite(samples)]
    if finite.size == 0 or not np.isfinite(observed):
        return {"p_value": None, "n_finite": int(finite.size), "sidedness": "two-sided"}
    recentered = finite - observed
    count = int(np.sum(np.abs(recentered) >= abs(observed)))
    p = (1 + count) / (finite.size + 1)
    return {
        "p_value": float(p),
        "n_finite": int(finite.size),
        "n_as_extreme": count,
        "sidedness": "two-sided",
        "null": "delta = 0",
        "procedure": (
            "Protein-cluster bootstrap. Recenter each replicate by subtracting the "
            "observed delta, then (1 + count of absolute recentered replicates at "
            "least as large as the absolute observed delta) / (n_finite + 1)."
        ),
    }


def split_census(path: Path) -> dict:
    df = pd.read_csv(path)
    out = {"path": str(path.relative_to(ROOT)), "n_rows": int(len(df)), "by_split": {}}
    for sp, sub in df.groupby("split"):
        rec = {
            "n_rows": int(len(sub)),
            "n_proteins": int(sub["target_id"].nunique()),
            "n_drugs": int(sub["drug_id"].nunique()),
        }
        prot_sizes = sub.groupby("target_id").size()
        drug_sizes = sub.groupby("drug_id").size()
        rec["rows_per_protein_min"] = int(prot_sizes.min())
        rec["rows_per_protein_median"] = float(prot_sizes.median())
        rec["rows_per_protein_max"] = int(prot_sizes.max())
        rec["proteins_per_drug_min"] = int(drug_sizes.min())
        rec["proteins_per_drug_median"] = float(drug_sizes.median())
        rec["proteins_per_drug_max"] = int(drug_sizes.max())
        rec["fraction_drugs_on_more_than_one_protein"] = float((drug_sizes > 1).mean())
        if "scaffold" in sub.columns:
            sc_sizes = sub.groupby("scaffold").size()
            rec["n_scaffolds"] = int(sub["scaffold"].nunique())
            rec["rows_per_scaffold_min"] = int(sc_sizes.min())
            rec["rows_per_scaffold_median"] = float(sc_sizes.median())
            rec["rows_per_scaffold_max"] = int(sc_sizes.max())
        out["by_split"][str(sp)] = rec
    return out


def main() -> None:
    pairs = pd.read_csv(ROOT / "data/raw/davis_pairs.csv")
    sid = pd.read_csv(ROOT / "data/splits/cold_protein_split.csv")
    keep = ["pair_id", "smiles", "sequence", "pkd", "drug_id", "target_id"]
    merged = sid[["pair_id", "split"]].merge(pairs[keep], on="pair_id", how="inner")
    if len(merged) != len(sid):
        raise SystemExit(f"split join mismatch {len(merged)} vs {len(sid)}")
    test = merged[merged["split"] == "test"].reset_index(drop=True)

    zp = np.load(ART / "models_v1.2/preds_cold_protein_protein.npz")
    zl = np.load(ART / "models_v1.2/preds_cold_protein_ligand_only.npz")
    y_true = np.asarray(zp["y_true"], dtype=float)
    y_prot = np.asarray(zp["y_pred"], dtype=float)
    y_lig = np.asarray(zl["y_pred"], dtype=float)
    y_true_l = np.asarray(zl["y_true"], dtype=float)
    if y_true.shape != y_true_l.shape or not np.array_equal(y_true, y_true_l):
        raise SystemExit("protein and ligand-only y_true disagree")
    if len(y_true) != len(test):
        raise SystemExit(f"length mismatch preds {len(y_true)} test {len(test)}")
    max_abs = float(np.max(np.abs(y_true - test["pkd"].to_numpy(dtype=float))))
    if max_abs > 1e-8:
        raise SystemExit(f"y_true does not match cold-protein test pkd in split order (max abs {max_abs})")

    proteins = test["target_id"].astype(str).to_numpy()
    drugs = test["drug_id"].astype(str).to_numpy()
    uniq_proteins = pd.unique(proteins)
    groups = [np.flatnonzero(proteins == p) for p in uniq_proteins]

    # Ligand-only score should be constant across proteins for a drug if the
    # prediction ignores the protein. Record the max within-drug range.
    lig_range = (
        pd.DataFrame({"drug": drugs, "pred": y_lig})
        .groupby("drug")["pred"]
        .agg(lambda s: float(s.max() - s.min()))
    )
    prot_range = (
        pd.DataFrame({"drug": drugs, "pred": y_prot})
        .groupby("drug")["pred"]
        .agg(lambda s: float(s.max() - s.min()))
    )

    rho_p = spearman_rho(y_true, y_prot)
    rho_l = spearman_rho(y_true, y_lig)
    delta = float(rho_p - rho_l)

    within = []
    for p, idx in zip(uniq_proteins, groups):
        rp = spearman_rho(y_true[idx], y_prot[idx])
        rl = spearman_rho(y_true[idx], y_lig[idx])
        within.append(
            {
                "target_id": str(p),
                "n": int(len(idx)),
                "spearman_protein": rp,
                "spearman_ligand_only": rl,
                "delta": float(rp - rl) if np.isfinite(rp) and np.isfinite(rl) else float("nan"),
            }
        )
    deltas_w = np.array([w["delta"] for w in within], dtype=float)
    finite_w = deltas_w[np.isfinite(deltas_w)]

    rng = np.random.default_rng(SEED)
    n_prot = len(groups)
    delta_b = np.empty(N_BOOT, dtype=float)
    rho_p_b = np.empty(N_BOOT, dtype=float)
    rho_l_b = np.empty(N_BOOT, dtype=float)
    for b in range(N_BOOT):
        draw = rng.integers(0, n_prot, size=n_prot)
        idx = np.concatenate([groups[i] for i in draw])
        rho_p_b[b] = spearman_rho(y_true[idx], y_prot[idx])
        rho_l_b[b] = spearman_rho(y_true[idx], y_lig[idx])
        delta_b[b] = rho_p_b[b] - rho_l_b[b]

    # Cited aggregates, not refit.
    v11 = json.loads((ART / "metrics_v1.1.json").read_text())
    cold = v11["splits"]["cold_protein"]
    hgb_delta = float(cold["protein"]["aac_dpc"]["spearman_rho"] - cold["ligand_only"]["spearman_rho"])
    kiba = json.loads((ART / "metrics_kiba.json").read_text())
    kcold = kiba["splits"]["cold_protein"]
    kiba_delta = float(kcold["aac_dpc"]["spearman_rho"] - kcold["ligand_only"]["spearman_rho"])
    unc = json.loads((ART / "metrics_uncertainty.json").read_text())
    screen = json.loads((ART / "metrics_screen_cold.json").read_text())
    cal = json.loads((ART / "metrics_calibration.json").read_text())

    scaffold_census = split_census(ROOT / "data/splits/scaffold_split.csv")
    cold_census = split_census(ROOT / "data/splits/cold_protein_split.csv")
    kiba_cold_census = split_census(ROOT / "data/splits_kiba/cold_protein_split.csv")
    kiba_scaf_census = split_census(ROOT / "data/splits_kiba/scaffold_split.csv")

    # Overlap checks on DAVIS cold split.
    cold_df = pd.read_csv(ROOT / "data/splits/cold_protein_split.csv")
    overlap = {}
    for a, b in (("train", "test"), ("train", "val"), ("val", "test")):
        pa = set(cold_df.loc[cold_df["split"] == a, "target_id"])
        pb = set(cold_df.loc[cold_df["split"] == b, "target_id"])
        da = set(cold_df.loc[cold_df["split"] == a, "drug_id"])
        db = set(cold_df.loc[cold_df["split"] == b, "drug_id"])
        overlap[f"{a}_vs_{b}"] = {
            "n_shared_proteins": int(len(pa & pb)),
            "n_shared_drugs": int(len(da & db)),
        }

    scaf = pd.read_csv(ROOT / "data/splits/scaffold_split.csv")
    scaffold_overlap = {}
    for a, b in (("train", "test"),):
        pa = set(scaf.loc[scaf["split"] == a, "target_id"])
        pb = set(scaf.loc[scaf["split"] == b, "target_id"])
        da = set(scaf.loc[scaf["split"] == a, "drug_id"])
        db = set(scaf.loc[scaf["split"] == b, "drug_id"])
        sa = set(scaf.loc[scaf["split"] == a, "scaffold"])
        sb = set(scaf.loc[scaf["split"] == b, "scaffold"])
        scaffold_overlap[f"{a}_vs_{b}"] = {
            "n_shared_proteins": int(len(pa & pb)),
            "n_shared_drugs": int(len(da & db)),
            "n_shared_scaffolds": int(len(sa & sb)),
        }

    n_cal_cold = int(0.3 * 5168)
    n_cal_scaf = int(0.3 * 5306)

    payload = {
        "seed": SEED,
        "n_bootstrap": N_BOOT,
        "models_refit": False,
        "alignment": {
            "cold_protein_test_rows": int(len(test)),
            "y_true_matches_split_order_pkd": True,
            "max_abs_y_true_minus_pkd": max_abs,
            "n_proteins": int(n_prot),
            "n_drugs": int(pd.unique(drugs).size),
            "ligand_only_score_range_within_drug_max": float(lig_range.max()),
            "protein_score_range_within_drug_max": float(prot_range.max()),
            "note": (
                "Predictions have no identifiers. Alignment is the cold-protein test "
                "frame in split-file order, which is the order run_chemprop_v1.2.py "
                "writes when dropped_smiles_test is 0. HistGBM per-pair scores are not in artifacts."
            ),
        },
        "recommended_test": {
            "name": "DAVIS cold-protein pooled delta Spearman, graph model plus joint composition versus graph ligand-only",
            "status": "computed",
            "unit": "protein",
            "not_computed_for": [
                "DAVIS HistGBM amino-acid plus dipeptide versus ligand-only: per-pair scores are not saved.",
                "KIBA HistGBM: per-pair scores are not saved.",
                "DAVIS scaffold graph comparison: predictions exist, but this review computes only the cold-protein test.",
            ],
            "estimand": (
                "Pooled Spearman rho of the protein-aware score minus pooled Spearman rho "
                "of the ligand-only score, on the cold-protein test rows."
            ),
            "resample": (
                "Draw the 76 test proteins with replacement. Keep every test row of each "
                "drawn protein (a protein drawn twice contributes its rows twice). "
                "Recompute both pooled Spearman correlations on the concatenated rows. "
                "The same draw enters both scores."
            ),
            "why_this_unit": (
                "The cold-protein split holds out whole proteins. Every test ligand is also "
                "in train, so the split does not sample new scaffolds. A pair draw treats "
                "5168 rows as independent. A scaffold draw would resample ligands, which "
                "is not the factor this split randomized."
            ),
            "point": {
                "spearman_protein": rho_p,
                "spearman_ligand_only": rho_l,
                "delta_spearman_protein_minus_ligand": delta,
            },
            "ci95": percentile_ci(delta_b),
            "spearman_protein_ci95": percentile_ci(rho_p_b),
            "spearman_ligand_only_ci95": percentile_ci(rho_l_b),
            "test": twosided_recentered_p(delta_b, delta),
            "descriptive_within_protein_not_a_test": {
                "note": (
                    "Unweighted mean of within-protein delta Spearman. Not the published "
                    "pooled estimand and not given a p-value in this review."
                ),
                "n_proteins_finite": int(finite_w.size),
                "mean_delta": float(np.mean(finite_w)) if finite_w.size else None,
                "median_delta": float(np.median(finite_w)) if finite_w.size else None,
                "std_delta": float(np.std(finite_w, ddof=1)) if finite_w.size > 1 else None,
                "n_positive_delta": int(np.sum(finite_w > 0)),
                "n_negative_delta": int(np.sum(finite_w < 0)),
                "n_zero_delta": int(np.sum(finite_w == 0)),
            },
        },
        "cited_point_estimates": {
            "histgbm_cold_aac_dpc_spearman": cold["protein"]["aac_dpc"]["spearman_rho"],
            "histgbm_cold_ligand_only_spearman": cold["ligand_only"]["spearman_rho"],
            "histgbm_cold_delta_spearman": hgb_delta,
            "histgbm_cold_aac_dpc_ef1": cold["protein"]["aac_dpc"]["ef_at_1pct"],
            "histgbm_cold_ligand_only_ef1": cold["ligand_only"]["ef_at_1pct"],
            "histgbm_cold_aac_dpc_ef5": cold["protein"]["aac_dpc"]["ef_at_5pct"],
            "histgbm_cold_ligand_only_ef5": cold["ligand_only"]["ef_at_5pct"],
            "kiba_label_note": kiba["label_note"],
            "kiba_binder_threshold": kiba["binder_threshold_kiba_score"],
            "kiba_cold_delta_spearman": kiba_delta,
            "kiba_cold_n": kcold["aac_dpc"]["n"],
            "kiba_cold_n_binders": kcold["aac_dpc"]["n_binders"],
            "slk_selection": screen["screen"]["selection"],
            "slk_n_known_binders": screen["recovery_protein"]["n_known_binders_in_library"],
            "slk_n_library": screen["recovery_protein"]["n_library"],
            "chemprop_calibration_n_cal_formula": {
                "cold_protein_n_cal": n_cal_cold,
                "cold_protein_n_eval_exported": next(
                    r["n_eval"] for r in cal["rows"] if r["split"] == "cold_protein" and r["model"] == "chemprop_protein_aac_dpc"
                ),
                "scaffold_n_cal": n_cal_scaf,
                "scaffold_n_eval_exported": next(
                    r["n_eval"] for r in cal["rows"] if r["split"] == "scaffold" and r["model"] == "chemprop_protein_aac_dpc"
                ),
            },
        },
        "dependence": {
            "davis_n_rows": int(len(pairs)),
            "davis_n_proteins": int(pairs["target_id"].nunique()),
            "davis_n_drugs": int(pairs["drug_id"].nunique()),
            "davis_complete_grid": bool(len(pairs) == pairs["target_id"].nunique() * pairs["drug_id"].nunique()),
            "cold_protein": cold_census,
            "scaffold": scaffold_census,
            "cold_protein_overlap": overlap,
            "scaffold_overlap_train_test": scaffold_overlap,
            "kiba_cold": kiba_cold_census,
            "kiba_scaffold": kiba_scaf_census,
            "davis_scaffolds_equal_drugs": bool(scaf["scaffold"].nunique() == scaf["drug_id"].nunique() == 68),
        },
        "existing_pair_bootstrap": {
            "source": "artifacts/metrics_uncertainty.json",
            "note": "Row bootstrap already exported. Not the recommended test. Not rerun.",
            "seed": unc["seed"],
            "n_bootstrap": unc["n_bootstrap"],
            "comparisons": [
                {
                    "name": c["name"],
                    "status": c["status"],
                    "kind": c["kind"],
                    "n_rows": None if c["result"] is None else c["result"].get("n_rows", c["result"].get("n_library")),
                    "delta_point": None
                    if c["result"] is None or c["kind"] != "spearman"
                    else c["result"]["delta_spearman_protein_minus_ligand"]["point"],
                    "delta_p": None
                    if c["result"] is None or c["kind"] != "spearman"
                    else c["result"]["delta_spearman_protein_minus_ligand"]["test"]["p_value"],
                    "delta_includes_0": None
                    if c["result"] is None or c["kind"] != "spearman"
                    else c["result"]["delta_spearman_protein_minus_ligand"]["ci95"]["includes_0"],
                }
                for c in unc["comparisons"]
            ],
        },
    }
    out = ART / "metrics_stat_review.json"
    out.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {out}")
    print("delta", delta)
    print("ci", payload["recommended_test"]["ci95"])
    print("p", payload["recommended_test"]["test"])
    print("within", payload["recommended_test"]["descriptive_within_protein_not_a_test"])
    print("hgb delta", hgb_delta)
    print("lig range max", float(lig_range.max()), "prot range max", float(prot_range.max()))


if __name__ == "__main__":
    main()
