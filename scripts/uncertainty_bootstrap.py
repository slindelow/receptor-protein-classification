#!/usr/bin/env python
"""Bootstrap intervals from saved per-pair predictions. Does not refit models.

Reads npz, csv, and json files under artifacts/. Writes
artifacts/metrics_uncertainty.json and artifacts/metrics_uncertainty.md.
The HistGBM cold-protein contrast is a kinase-cluster bootstrap of
artifacts/scores_histgbm_cold_protein.csv. Previously computed graph
pair-bootstrap and SLK recovery blocks are kept as stored when present.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from vs.metrics_lib import enrichment_factor  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = ROOT / "artifacts"
SEED = 42
N_BOOT = 1000


def spearman_rho(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Average-rank Spearman, matching scipy.stats.spearmanr on finite pairs."""
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
        return {
            "low": None,
            "high": None,
            "n_finite": 0,
            "includes_0": None,
        }
    low, high = np.percentile(finite, [2.5, 97.5])
    return {
        "low": float(low),
        "high": float(high),
        "n_finite": int(finite.size),
        "includes_0": bool(low <= 0.0 <= high),
    }


def twosided_recentered_p(samples: np.ndarray, observed: float) -> dict:
    """Two-sided bootstrap p for a parameter of zero.

    Each replicate is recentered by the observed value. The p-value is
    (1 + count of |recentered| >= |observed|) / (n + 1).
    """
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
            "Recenter each bootstrap replicate by subtracting the observed delta, "
            "then (1 + count of absolute recentered replicates at least as large "
            "as the absolute observed delta) / (n_finite + 1)."
        ),
    }


def bootstrap_spearman_pair(
    y_true: np.ndarray,
    y_protein: np.ndarray,
    y_ligand: np.ndarray,
    n_boot: int,
    seed: int,
) -> dict:
    y_true = np.asarray(y_true, dtype=float)
    y_protein = np.asarray(y_protein, dtype=float)
    y_ligand = np.asarray(y_ligand, dtype=float)
    if not (len(y_true) == len(y_protein) == len(y_ligand)):
        raise ValueError("prediction vectors differ in length")
    if not np.array_equal(np.asarray(y_true), np.asarray(y_true)):
        raise ValueError("non-finite check failed")
    n = len(y_true)
    rho_p = spearman_rho(y_true, y_protein)
    rho_l = spearman_rho(y_true, y_ligand)
    delta = float(rho_p - rho_l)
    rng = np.random.default_rng(seed)
    rho_p_b = np.empty(n_boot, dtype=float)
    rho_l_b = np.empty(n_boot, dtype=float)
    delta_b = np.empty(n_boot, dtype=float)
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        rho_p_b[b] = spearman_rho(y_true[idx], y_protein[idx])
        rho_l_b[b] = spearman_rho(y_true[idx], y_ligand[idx])
        delta_b[b] = rho_p_b[b] - rho_l_b[b]
    return {
        "n_rows": int(n),
        "n_bootstrap": int(n_boot),
        "seed": int(seed),
        "resample": "pairs with replacement, shared row indices for both scores",
        "spearman_protein": {
            "point": float(rho_p),
            "ci95": percentile_ci(rho_p_b),
        },
        "spearman_ligand_only": {
            "point": float(rho_l),
            "ci95": percentile_ci(rho_l_b),
        },
        "delta_spearman_protein_minus_ligand": {
            "point": delta,
            "ci95": percentile_ci(delta_b),
            "test": twosided_recentered_p(delta_b, delta),
        },
    }


def ranks_high_is_best(scores: np.ndarray, tie_codes: np.ndarray) -> np.ndarray:
    """Rank 1 is the highest finite score. Ties break by ascending tie code. Non-finite last."""
    scores = np.asarray(scores, dtype=float)
    nonfinite = ~np.isfinite(scores)
    sort_score = np.where(nonfinite, 0.0, -scores)
    order = np.lexsort((tie_codes, sort_score, nonfinite.astype(np.int8)))
    ranks = np.empty(len(scores), dtype=np.int32)
    ranks[order] = np.arange(1, len(scores) + 1, dtype=np.int32)
    return ranks


def recovery_from_ranks(ranks: np.ndarray, is_binder: np.ndarray, n: int) -> dict:
    k1 = max(1, int(np.ceil(0.01 * n)))
    binder_ranks = ranks[is_binder]
    n_known = int(is_binder.sum())
    if n_known == 0:
        return {
            "n_known": 0,
            "fraction_top_1pct": float("nan"),
            "fraction_top_500": float("nan"),
            "mean_rank": float("nan"),
            "top1_k": k1,
        }
    return {
        "n_known": n_known,
        "fraction_top_1pct": float(np.sum(binder_ranks <= k1) / n_known),
        "fraction_top_500": float(np.sum(binder_ranks <= 500) / n_known),
        "mean_rank": float(np.mean(binder_ranks)),
        "top1_k": k1,
    }


def bootstrap_screen_recovery(csv_path: Path, n_boot: int, seed: int) -> dict:
    df = pd.read_csv(csv_path)
    required = {"pred_pkd", "ligand_only_pred_pkd", "is_known_binder", "smiles", "rank", "ligand_only_rank"}
    missing = required - set(df.columns)
    if missing:
        return {
            "status": "skipped",
            "reason": f"{csv_path.name} is missing columns {sorted(missing)}. Recovery difference was not imputed.",
            "path": str(csv_path.relative_to(ROOT)),
        }
    smiles = df["smiles"].astype(str).to_numpy()
    _, tie_codes = np.unique(smiles, return_inverse=True)
    pred_p = df["pred_pkd"].to_numpy(dtype=float)
    pred_l = df["ligand_only_pred_pkd"].to_numpy(dtype=float)
    is_binder = df["is_known_binder"].astype(str).str.lower().isin(["true", "1", "yes"]).to_numpy()
    n = len(df)
    ranks_p = ranks_high_is_best(pred_p, tie_codes)
    ranks_l = ranks_high_is_best(pred_l, tie_codes)
    # The saved rank columns are the reference. Do not continue if the tie break disagrees.
    if not np.array_equal(ranks_p, df["rank"].to_numpy(dtype=np.int32)):
        n_mismatch = int(np.sum(ranks_p != df["rank"].to_numpy(dtype=np.int32)))
        return {
            "status": "skipped",
            "reason": (
                f"Recomputed protein ranks disagree with the saved rank column on {n_mismatch} rows. "
                "Recovery was not recomputed under a different tie break."
            ),
            "path": str(csv_path.relative_to(ROOT)),
        }
    if not np.array_equal(ranks_l, df["ligand_only_rank"].to_numpy(dtype=np.int32)):
        n_mismatch = int(np.sum(ranks_l != df["ligand_only_rank"].to_numpy(dtype=np.int32)))
        return {
            "status": "skipped",
            "reason": (
                f"Recomputed ligand-only ranks disagree with the saved ligand_only_rank column on {n_mismatch} rows. "
                "Recovery was not recomputed under a different tie break."
            ),
            "path": str(csv_path.relative_to(ROOT)),
        }
    base_p = recovery_from_ranks(ranks_p, is_binder, n)
    base_l = recovery_from_ranks(ranks_l, is_binder, n)
    point = {
        "fraction_top_1pct": base_p["fraction_top_1pct"] - base_l["fraction_top_1pct"],
        "fraction_top_500": base_p["fraction_top_500"] - base_l["fraction_top_500"],
        "mean_rank_protein_minus_ligand": base_p["mean_rank"] - base_l["mean_rank"],
        "mean_rank_improvement_ligand_minus_protein": base_l["mean_rank"] - base_p["mean_rank"],
    }
    rng = np.random.default_rng(seed)
    keys = (
        "fraction_top_1pct",
        "fraction_top_500",
        "mean_rank_protein_minus_ligand",
        "mean_rank_improvement_ligand_minus_protein",
    )
    boots = {k: np.empty(n_boot, dtype=float) for k in keys}
    n_empty = 0
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        rp = ranks_high_is_best(pred_p[idx], tie_codes[idx])
        rl = ranks_high_is_best(pred_l[idx], tie_codes[idx])
        rec_p = recovery_from_ranks(rp, is_binder[idx], n)
        rec_l = recovery_from_ranks(rl, is_binder[idx], n)
        if rec_p["n_known"] == 0 or rec_l["n_known"] == 0:
            n_empty += 1
            for k in keys:
                boots[k][b] = np.nan
            continue
        boots["fraction_top_1pct"][b] = rec_p["fraction_top_1pct"] - rec_l["fraction_top_1pct"]
        boots["fraction_top_500"][b] = rec_p["fraction_top_500"] - rec_l["fraction_top_500"]
        boots["mean_rank_protein_minus_ligand"][b] = rec_p["mean_rank"] - rec_l["mean_rank"]
        boots["mean_rank_improvement_ligand_minus_protein"][b] = rec_l["mean_rank"] - rec_p["mean_rank"]
    deltas = {}
    for k in keys:
        deltas[k] = {
            "point": float(point[k]),
            "ci95": percentile_ci(boots[k]),
            "test": twosided_recentered_p(boots[k], float(point[k])),
            "direction": (
                "protein minus ligand-only; positive favors the protein score"
                if k != "mean_rank_protein_minus_ligand"
                else "protein mean rank minus ligand-only mean rank; negative favors the protein score"
            ),
        }
    return {
        "status": "computed",
        "path": str(csv_path.relative_to(ROOT)),
        "comparison": "SLK library recovery, HistGBM joint composition minus ligand-only",
        "n_library": int(n),
        "n_known_binders": int(is_binder.sum()),
        "n_bootstrap": int(n_boot),
        "seed": int(seed),
        "n_replicates_with_no_known_binder": int(n_empty),
        "resample": (
            "library rows with replacement; ranks recomputed inside each replicate; "
            "top 1% cutoff is ceil(0.01 * n) and the top-500 cutoff stays 500"
        ),
        "point_recovery_protein": base_p,
        "point_recovery_ligand_only": base_l,
        "deltas": deltas,
        "note": (
            "Continuous per-pair pKd is not stored for the library, so delta Spearman was not computed. "
            "The interval is on the recovery difference only."
        ),
    }


def load_aligned_npz(protein_path: Path, ligand_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    zp = np.load(protein_path)
    zl = np.load(ligand_path)
    if "y_true" not in zp.files or "y_pred" not in zp.files:
        raise ValueError(f"{protein_path} missing y_true/y_pred")
    if "y_true" not in zl.files or "y_pred" not in zl.files:
        raise ValueError(f"{ligand_path} missing y_true/y_pred")
    yt_p = np.asarray(zp["y_true"], dtype=float)
    yt_l = np.asarray(zl["y_true"], dtype=float)
    if yt_p.shape != yt_l.shape or not np.array_equal(yt_p, yt_l):
        raise ValueError(f"y_true mismatch between {protein_path.name} and {ligand_path.name}")
    return yt_p, np.asarray(zp["y_pred"], dtype=float), np.asarray(zl["y_pred"], dtype=float)



def bootstrap_kinase_cluster(
    y_true: np.ndarray,
    y_protein: np.ndarray,
    y_ligand: np.ndarray,
    protein_ids: np.ndarray,
    n_boot: int,
    seed: int,
    binder_threshold: float = 7.0,
    ef_fraction: float = 0.01,
) -> dict:
    """Resample held-out kinases, not pairs. Same recentering as the pair bootstrap.

    A kinase drawn twice contributes its pairs twice. Both scores use that draw.
    Spearman and EF@1% are recomputed on the concatenated pairs.
    """
    y_true = np.asarray(y_true, dtype=float)
    y_protein = np.asarray(y_protein, dtype=float)
    y_ligand = np.asarray(y_ligand, dtype=float)
    protein_ids = np.asarray(protein_ids)
    if not (len(y_true) == len(y_protein) == len(y_ligand) == len(protein_ids)):
        raise ValueError("kinase-cluster inputs differ in length")
    # Appearance order, matching scripts/stat_review_cold_protein.py. Do not sort.
    uniq = pd.unique(protein_ids)
    groups = [np.flatnonzero(protein_ids == p) for p in uniq]
    y_bin = (y_true >= binder_threshold).astype(int)
    rho_p = spearman_rho(y_true, y_protein)
    rho_l = spearman_rho(y_true, y_ligand)
    delta = float(rho_p - rho_l)
    ef_p = float(enrichment_factor(y_bin, y_protein, ef_fraction))
    ef_l = float(enrichment_factor(y_bin, y_ligand, ef_fraction))
    delta_ef = float(ef_p - ef_l)
    rng = np.random.default_rng(seed)
    n_prot = len(groups)
    rho_p_b = np.empty(n_boot, dtype=float)
    rho_l_b = np.empty(n_boot, dtype=float)
    delta_b = np.empty(n_boot, dtype=float)
    ef_p_b = np.empty(n_boot, dtype=float)
    ef_l_b = np.empty(n_boot, dtype=float)
    delta_ef_b = np.empty(n_boot, dtype=float)
    for b in range(n_boot):
        draw = rng.integers(0, n_prot, size=n_prot)
        idx = np.concatenate([groups[int(i)] for i in draw])
        rho_p_b[b] = spearman_rho(y_true[idx], y_protein[idx])
        rho_l_b[b] = spearman_rho(y_true[idx], y_ligand[idx])
        delta_b[b] = rho_p_b[b] - rho_l_b[b]
        ef_p_b[b] = enrichment_factor(y_bin[idx], y_protein[idx], ef_fraction)
        ef_l_b[b] = enrichment_factor(y_bin[idx], y_ligand[idx], ef_fraction)
        delta_ef_b[b] = ef_p_b[b] - ef_l_b[b]
    sizes = np.array([len(g) for g in groups], dtype=int)
    return {
        "n_rows": int(len(y_true)),
        "n_kinases": int(n_prot),
        "rows_per_kinase_min": int(sizes.min()) if sizes.size else 0,
        "rows_per_kinase_max": int(sizes.max()) if sizes.size else 0,
        "n_bootstrap": int(n_boot),
        "seed": int(seed),
        "resample_unit": "held-out kinase",
        "resample": (
            "Draw the held-out kinases with replacement. Keep every test pair of each "
            "drawn kinase (a kinase drawn twice contributes its pairs twice). Recompute "
            "pooled Spearman and EF@1% on the concatenated pairs. The same draw enters both scores."
        ),
        "binder_threshold": float(binder_threshold),
        "ef_fraction": float(ef_fraction),
        "spearman_protein": {"point": float(rho_p), "ci95": percentile_ci(rho_p_b)},
        "spearman_ligand_only": {"point": float(rho_l), "ci95": percentile_ci(rho_l_b)},
        "delta_spearman_protein_minus_ligand": {
            "point": delta,
            "ci95": percentile_ci(delta_b),
            "test": twosided_recentered_p(delta_b, delta),
        },
        "ef_at_1pct_protein": {"point": ef_p, "ci95": percentile_ci(ef_p_b)},
        "ef_at_1pct_ligand_only": {"point": ef_l, "ci95": percentile_ci(ef_l_b)},
        "delta_ef_at_1pct_protein_minus_ligand": {
            "point": delta_ef,
            "ci95": percentile_ci(delta_ef_b),
            "test": twosided_recentered_p(delta_ef_b, delta_ef),
        },
    }


HISTGBM_COLD_NAME = "DAVIS cold-protein, HistGBM amino-acid plus dipeptide versus ligand-only"
HISTGBM_SCORE_CSV = ARTIFACTS / "scores_histgbm_cold_protein.csv"
HISTGBM_CHECK_JSON = ARTIFACTS / "scores_histgbm_cold_protein_check.json"
PUBLISHED_V11 = ARTIFACTS / "metrics_v1.1.json"
MATCH_TOL = 1e-4


def _published_cold_histgbm() -> dict:
    published = json.loads(PUBLISHED_V11.read_text())
    cold = published["splits"]["cold_protein"]
    return {
        "ligand_only_spearman": float(cold["ligand_only"]["spearman_rho"]),
        "aac_dpc_spearman": float(cold["protein"]["aac_dpc"]["spearman_rho"]),
        "ligand_only_ef_at_1pct": float(cold["ligand_only"]["ef_at_1pct"]),
        "aac_dpc_ef_at_1pct": float(cold["protein"]["aac_dpc"]["ef_at_1pct"]),
        "mode_spearman": {
            "ligand_only": float(cold["ligand_only"]["spearman_rho"]),
            "aac": float(cold["protein"]["aac"]["spearman_rho"]),
            "dpc": float(cold["protein"]["dpc"]["spearman_rho"]),
            "aac_dpc": float(cold["protein"]["aac_dpc"]["spearman_rho"]),
            "esm2": float(cold["protein"]["esm2"]["spearman_rho"]),
        },
    }


def histgbm_kinase_comparison() -> dict:
    """Kinase-cluster interval for HistGBM AAC+DPC minus ligand-only, or a recorded mismatch."""
    name = HISTGBM_COLD_NAME
    if not HISTGBM_SCORE_CSV.exists():
        return {
            "name": name,
            "status": "skipped",
            "reason": (
                "Skipped. Per-pair scores for this comparison are not in artifacts/. "
                "Aggregate Spearman rho values in metrics_v1.1.json were not turned into an interval."
            ),
            "kind": None,
            "sources": None,
            "result": None,
        }
    df = pd.read_csv(HISTGBM_SCORE_CSV)
    required = {"target_id", "pkd", "pred_ligand_only", "pred_aac_dpc"}
    missing = required - set(df.columns)
    if missing:
        return {
            "name": name,
            "status": "mismatch",
            "reason": f"{HISTGBM_SCORE_CSV.name} is missing columns {sorted(missing)}. No interval was computed.",
            "kind": None,
            "sources": str(HISTGBM_SCORE_CSV.relative_to(ROOT)),
            "result": None,
        }
    y = df["pkd"].to_numpy(dtype=float)
    y_p = df["pred_aac_dpc"].to_numpy(dtype=float)
    y_l = df["pred_ligand_only"].to_numpy(dtype=float)
    pub = _published_cold_histgbm()
    recomputed = {
        "ligand_only_spearman": spearman_rho(y, y_l),
        "aac_dpc_spearman": spearman_rho(y, y_p),
        "ligand_only_ef_at_1pct": float(enrichment_factor((y >= 7.0).astype(int), y_l, 0.01)),
        "aac_dpc_ef_at_1pct": float(enrichment_factor((y >= 7.0).astype(int), y_p, 0.01)),
    }
    match_detail = {}
    matched = True
    for key, got in recomputed.items():
        diff = abs(float(got) - pub[key])
        ok = bool(np.isfinite(got) and diff <= MATCH_TOL)
        match_detail[key] = {"recomputed": float(got), "published": pub[key], "abs_diff": float(diff), "match": ok}
        matched = matched and ok
    check = {}
    if HISTGBM_CHECK_JSON.exists():
        check = json.loads(HISTGBM_CHECK_JSON.read_text())
    models_source = "not recorded"
    if check.get("refit") is False and check.get("models") == "loaded":
        models_source = "loaded from artifacts/models_v1.1 joblibs, not refit"
    elif check.get("refit") is True:
        models_source = "refit"
    if not matched:
        return {
            "name": name,
            "status": "mismatch",
            "reason": (
                "Recomputed full-test Spearman or EF@1% does not match metrics_v1.1.json within 1e-4. "
                "No bootstrap interval was computed. Published metrics were not overwritten."
            ),
            "kind": None,
            "sources": str(HISTGBM_SCORE_CSV.relative_to(ROOT)),
            "result": {
                "point_estimates_matched": False,
                "match_tolerance": MATCH_TOL,
                "match_detail": match_detail,
                "models_source": models_source,
                "n_rows": int(len(df)),
                "n_kinases": int(df["target_id"].nunique()),
            },
        }
    result = bootstrap_kinase_cluster(
        y, y_p, y_l, df["target_id"].astype(str).to_numpy(), N_BOOT, SEED
    )
    result["models_source"] = models_source
    result["point_estimates_matched"] = True
    result["match_tolerance"] = MATCH_TOL
    result["match_detail"] = match_detail
    result["published_mode_spearman"] = pub["mode_spearman"]
    result["pre_specified_contrast"] = "AAC+DPC versus ligand-only"
    modes = pub["mode_spearman"]
    result["selection_note"] = (
        "AAC+DPC was the best cold-protein Spearman among four modes scored on this same split "
        f"(ligand-only {modes['ligand_only']:.6f}, AAC {modes['aac']:.6f}, DPC {modes['dpc']:.6f}, "
        f"AAC+DPC {modes['aac_dpc']:.6f}). ESM was a fifth comparison and was essentially tied with AAC+DPC "
        f"(ESM Spearman {modes['esm2']:.6f}). The contrast tested here is the pre-specified AAC+DPC versus "
        "ligand-only contrast. The interval is not adjusted for the other modes."
    )
    sp_zero = result["delta_spearman_protein_minus_ligand"]["ci95"]["includes_0"]
    ef_zero = result["delta_ef_at_1pct_protein_minus_ligand"]["ci95"]["includes_0"]
    def _zero_sentence(includes, label):
        if includes:
            return f"The {label} interval includes zero."
        return (
            f"The {label} interval excludes zero. That exclusion is not multiplicity-adjusted: "
            "AAC+DPC was chosen as the best of four modes on this split, with ESM a fifth comparison "
            "essentially tied with it."
        )
    result["interval_note"] = (
        _zero_sentence(sp_zero, "delta Spearman")
        + " "
        + _zero_sentence(ef_zero, "delta EF@1%")
        + " This file does not call either contrast significant."
    )
    return {
        "name": name,
        "status": "computed",
        "reason": None,
        "kind": "kinase_spearman",
        "sources": str(HISTGBM_SCORE_CSV.relative_to(ROOT)),
        "result": result,
    }


def preserved_computed(previous: dict | None, name: str) -> dict | None:
    if not previous:
        return None
    for comp in previous.get("comparisons", []):
        if comp.get("name") == name and comp.get("status") == "computed":
            return comp
    return None


def discover_prediction_files() -> list[str]:
    found = []
    for path in sorted(ARTIFACTS.rglob("*")):
        if path.suffix.lower() in {".npz", ".csv", ".json"} and path.is_file():
            found.append(str(path.relative_to(ROOT)))
    return found


def fmt(x, digits=4) -> str:
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    return f"{float(x):.{digits}f}"


def fmt_p(x) -> str:
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    return f"{float(x):.4f}"


def yes_no(flag) -> str:
    if flag is None:
        return "NA"
    return "yes" if flag else "no"


def markdown_report(payload: dict) -> str:
    lines = [
        "# Uncertainty from saved predictions",
        "",
        "Numbers in this file were printed by `scripts/uncertainty_bootstrap.py`. Models were not refit.",
        "",
        (
            f"Seed {payload['seed']}, {payload['n_bootstrap']} resamples, percentile 95% interval. "
            "The two-sided p-value recenters each replicate by the observed delta and uses a plus-one correction. "
            "Graph Spearman blocks and the SLK recovery block resample rows. The HistGBM cold-protein block, "
            "when computed, resamples held-out kinases. Delta is the protein score minus the ligand-only score "
            "unless a row says otherwise."
        ),
        "",
        "## What was saved",
        "",
    ]
    for item in payload["inventory_notes"]:
        lines.append(f"- {item}")
    lines.append("")
    lines.append("## Computed comparisons")
    lines.append("")
    for comp in payload["comparisons"]:
        lines.append(f"### {comp['name']}")
        lines.append("")
        lines.append(f"Status: {comp['status']}.")
        if comp["status"] != "computed":
            lines.append("")
            lines.append(comp["reason"])
            lines.append("")
            continue
        lines.append("")
        lines.append(f"Source: `{comp['sources']}`.")
        lines.append("")
        if comp["kind"] == "kinase_spearman":
            res = comp["result"]
            d = res["delta_spearman_protein_minus_ligand"]
            e = res["delta_ef_at_1pct_protein_minus_ligand"]
            sp = res["spearman_protein"]
            sl = res["spearman_ligand_only"]
            lines.append(
                f"Resample unit: {res['resample_unit']}. n kinases: {res['n_kinases']}. "
                f"n pairs: {res['n_rows']}. Seed: {res['seed']}. n resamples: {res['n_bootstrap']}. "
                f"Models: {res['models_source']}. "
                f"Point estimates matched published metrics_v1.1.json within {res['match_tolerance']}: "
                f"{'yes' if res['point_estimates_matched'] else 'no'}."
            )
            lines.append("")
            lines.append(res["resample"])
            lines.append("")
            lines.append("| Quantity | Point | 95% CI low | 95% CI high | Includes 0 | Two-sided p |")
            lines.append("|---|---:|---:|---:|---|---:|")
            lines.append(
                f"| Spearman rho, AAC+DPC | {fmt(sp['point'])} | {fmt(sp['ci95']['low'])} | {fmt(sp['ci95']['high'])} | {yes_no(sp['ci95']['includes_0'])} |  |"
            )
            lines.append(
                f"| Spearman rho, ligand-only | {fmt(sl['point'])} | {fmt(sl['ci95']['low'])} | {fmt(sl['ci95']['high'])} | {yes_no(sl['ci95']['includes_0'])} |  |"
            )
            lines.append(
                f"| Delta Spearman (AAC+DPC minus ligand-only) | {fmt(d['point'])} | {fmt(d['ci95']['low'])} | {fmt(d['ci95']['high'])} | {yes_no(d['ci95']['includes_0'])} | {fmt_p(d['test']['p_value'])} |"
            )
            lines.append(
                f"| Delta EF@1% (AAC+DPC minus ligand-only) | {fmt(e['point'])} | {fmt(e['ci95']['low'])} | {fmt(e['ci95']['high'])} | {yes_no(e['ci95']['includes_0'])} | {fmt_p(e['test']['p_value'])} |"
            )
            lines.append("")
            lines.append(
                f"Unrounded delta Spearman: {d['point']}. Unrounded delta EF@1%: {e['point']}. "
                f"Finite delta-Spearman replicates: {d['ci95']['n_finite']}."
            )
            lines.append("")
            lines.append(res["selection_note"])
            lines.append("")
            lines.append(res["interval_note"])
            lines.append("")
        elif comp["kind"] == "spearman":
            sp = comp["result"]["spearman_protein"]
            sl = comp["result"]["spearman_ligand_only"]
            d = comp["result"]["delta_spearman_protein_minus_ligand"]
            lines.append("| Quantity | Point | 95% CI low | 95% CI high | Includes 0 | Two-sided p |")
            lines.append("|---|---:|---:|---:|---|---:|")
            lines.append(
                f"| Spearman rho, protein | {fmt(sp['point'])} | {fmt(sp['ci95']['low'])} | {fmt(sp['ci95']['high'])} | {yes_no(sp['ci95']['includes_0'])} |  |"
            )
            lines.append(
                f"| Spearman rho, ligand-only | {fmt(sl['point'])} | {fmt(sl['ci95']['low'])} | {fmt(sl['ci95']['high'])} | {yes_no(sl['ci95']['includes_0'])} |  |"
            )
            lines.append(
                f"| Delta Spearman (protein minus ligand-only) | {fmt(d['point'])} | {fmt(d['ci95']['low'])} | {fmt(d['ci95']['high'])} | {yes_no(d['ci95']['includes_0'])} | {fmt_p(d['test']['p_value'])} |"
            )
            lines.append("")
            lines.append(
                f"Rows: {comp['result']['n_rows']}. Finite delta replicates: {d['ci95']['n_finite']}."
            )
            lines.append("")
        elif comp["kind"] == "recovery":
            res = comp["result"]
            lines.append("| Quantity | Point | 95% CI low | 95% CI high | Includes 0 | Two-sided p |")
            lines.append("|---|---:|---:|---:|---|---:|")
            for key, label in (
                ("fraction_top_1pct", "Delta fraction in top 1% (protein minus ligand-only)"),
                ("fraction_top_500", "Delta fraction in top 500 (protein minus ligand-only)"),
                (
                    "mean_rank_improvement_ligand_minus_protein",
                    "Mean-rank improvement (ligand-only minus protein; positive favors protein)",
                ),
            ):
                d = res["deltas"][key]
                lines.append(
                    f"| {label} | {fmt(d['point'])} | {fmt(d['ci95']['low'])} | {fmt(d['ci95']['high'])} | {yes_no(d['ci95']['includes_0'])} | {fmt_p(d['test']['p_value'])} |"
                )
            lines.append("")
            lines.append(res["note"])
            lines.append("")
            lines.append(
                f"Library rows: {res['n_library']}. Known binders: {res['n_known_binders']}. "
                f"Replicates with no known binder: {res['n_replicates_with_no_known_binder']}."
            )
            lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def update_stat_review(comp: dict) -> None:
    """Append or replace section 6 from the computed comparison. No hand-typed statistics."""
    path = ARTIFACTS / "STAT_REVIEW.md"
    heading = "## 6. HistGBM cold-protein kinase-cluster bootstrap"
    lines = [heading, ""]
    lines.append(
        "This section supersedes the statement in section 3 that the histogram lift cannot be tested. "
        "It was written by `scripts/uncertainty_bootstrap.py` from the comparison block in "
        "`artifacts/metrics_uncertainty.json`. The graph protein-cluster result already in "
        "`artifacts/metrics_stat_review.json` was not recomputed."
    )
    lines.append("")
    if comp["status"] == "skipped":
        lines.append(comp["reason"])
        lines.append("")
    elif comp["status"] == "mismatch":
        lines.append(comp["reason"])
        lines.append("")
        if comp.get("result"):
            lines.append("Match detail:")
            lines.append("")
            lines.append("```")
            lines.append(json.dumps(comp["result"].get("match_detail"), indent=2))
            lines.append("```")
            lines.append("")
    elif comp["status"] == "computed":
        res = comp["result"]
        d = res["delta_spearman_protein_minus_ligand"]
        e = res["delta_ef_at_1pct_protein_minus_ligand"]
        lines.append(
            f"Score file: `{comp['sources']}`. Models: {res['models_source']}. "
            f"Point estimates matched `metrics_v1.1.json` within {res['match_tolerance']}: "
            f"{'yes' if res['point_estimates_matched'] else 'no'}."
        )
        lines.append("")
        lines.append(
            f"Resample unit: {res['resample_unit']}. n kinases: {res['n_kinases']}. "
            f"n pairs: {res['n_rows']}. Seed: {res['seed']}. n resamples: {res['n_bootstrap']}."
        )
        lines.append("")
        lines.append(res["resample"])
        lines.append("")
        lines.append(
            f"Delta Spearman (AAC+DPC minus ligand-only): point {d['point']}, "
            f"95% percentile interval {d['ci95']['low']} to {d['ci95']['high']}, "
            f"includes 0: {d['ci95']['includes_0']}, "
            f"two-sided recentered p {d['test']['p_value']} "
            f"({d['test']['n_as_extreme']} of {d['test']['n_finite']} recentered replicates at least as large, plus one)."
        )
        lines.append("")
        lines.append(
            f"Component Spearman on the original test: AAC+DPC {res['spearman_protein']['point']}, "
            f"ligand-only {res['spearman_ligand_only']['point']}."
        )
        lines.append("")
        lines.append(
            f"Delta EF@1% (AAC+DPC minus ligand-only): point {e['point']}, "
            f"95% percentile interval {e['ci95']['low']} to {e['ci95']['high']}, "
            f"includes 0: {e['ci95']['includes_0']}, "
            f"two-sided recentered p {e['test']['p_value']}."
        )
        lines.append("")
        lines.append(res["selection_note"])
        lines.append("")
        lines.append(res["interval_note"])
        lines.append("")
        lines.append(
            "Rerun from the repository root: "
            "`.venv/bin/python scripts/score_histgbm_cold_protein.py && "
            ".venv/bin/python scripts/uncertainty_bootstrap.py`."
        )
        lines.append("")
    else:
        lines.append(f"Status: {comp['status']}.")
        lines.append("")
    existing = path.read_text() if path.exists() else ""
    if heading in existing:
        existing = existing[: existing.index(heading)].rstrip() + "\n\n"
    else:
        existing = existing.rstrip() + "\n\n"
    path.write_text(existing + "\n".join(lines).rstrip() + "\n")


def main() -> None:
    previous = None
    previous_path = ARTIFACTS / "metrics_uncertainty.json"
    if previous_path.exists():
        previous = json.loads(previous_path.read_text())
    notes = []
    files = discover_prediction_files()
    npz_files = [f for f in files if f.endswith(".npz")]
    notes.append(
        "NPZ prediction files found: " + (", ".join(f"`{f}`" for f in npz_files) if npz_files else "none") + "."
    )
    if HISTGBM_SCORE_CSV.exists():
        notes.append(
            "Per-pair HistGBM AAC+DPC and ligand-only scores for the DAVIS cold-protein test are in "
            f"`{HISTGBM_SCORE_CSV.relative_to(ROOT)}`. Models were loaded from artifacts/models_v1.1, not refit. "
            "Scaffold HistGBM scores and KIBA per-pair scores are still absent. "
            "The warm LCK screen was not bootstrapped."
        )
    else:
        notes.append(
            "No npz, csv, or json under artifacts/ stores per-pair HistGBM amino-acid-plus-dipeptide "
            "and ligand-only scores on the DAVIS cold-protein test or the DAVIS scaffold test. "
            "metrics_v1.1.json stores aggregate metrics only. Those models were not refit."
        )
    notes.append(
        "No per-pair KIBA scores were stored. metrics_kiba.json stores aggregate metrics only. "
        "KIBA cold-protein delta Spearman was not resampled."
    )
    notes.append(
        "artifacts/models_v1.2/preds_*.npz are the Chemprop runs named in metrics_v1.2.json "
        "(graph model, joint composition versus ligand-only), not the histogram model."
    )

    comparisons = []

    def skipped(name: str, reason: str) -> None:
        comparisons.append({"name": name, "status": "skipped", "reason": reason, "kind": None, "sources": None, "result": None})

    comparisons.append(histgbm_kinase_comparison())
    skipped(
        "DAVIS scaffold, HistGBM amino-acid plus dipeptide versus ligand-only",
        "Skipped. Per-pair scores for this comparison are not in artifacts/.",
    )
    skipped(
        "KIBA cold-protein, HistGBM joint composition versus ligand-only",
        "Skipped. Per-pair KIBA scores are not in artifacts/. The recovery-style or delta-Spearman interval was not imputed from aggregate metrics.",
    )

    chemprop_specs = [
        (
            "DAVIS cold-protein, graph model plus joint composition versus graph ligand-only",
            ARTIFACTS / "models_v1.2" / "preds_cold_protein_protein.npz",
            ARTIFACTS / "models_v1.2" / "preds_cold_protein_ligand_only.npz",
        ),
        (
            "DAVIS scaffold, graph model plus joint composition versus graph ligand-only",
            ARTIFACTS / "models_v1.2" / "preds_scaffold_protein.npz",
            ARTIFACTS / "models_v1.2" / "preds_scaffold_ligand_only.npz",
        ),
    ]
    for name, prot, lig in chemprop_specs:
        kept = preserved_computed(previous, name)
        if kept is not None:
            comparisons.append(kept)
            continue
        y_true, y_p, y_l = load_aligned_npz(prot, lig)
        result = bootstrap_spearman_pair(y_true, y_p, y_l, N_BOOT, SEED)
        comparisons.append(
            {
                "name": name,
                "status": "computed",
                "reason": None,
                "kind": "spearman",
                "sources": f"{prot.relative_to(ROOT)} and {lig.relative_to(ROOT)}",
                "result": result,
            }
        )

    screen_path = ARTIFACTS / "screen_ranked_cold.csv"
    screen_name = "SLK library recovery, HistGBM joint composition versus ligand-only"
    kept_screen = preserved_computed(previous, screen_name)
    screen = None if kept_screen is not None else bootstrap_screen_recovery(screen_path, N_BOOT, SEED)
    if kept_screen is not None:
        comparisons.append(kept_screen)
    elif screen["status"] != "computed":
        comparisons.append(
            {
                "name": "SLK library recovery, HistGBM joint composition versus ligand-only",
                "status": "skipped",
                "reason": screen["reason"],
                "kind": None,
                "sources": screen.get("path"),
                "result": None,
            }
        )
    else:
        comparisons.append(
            {
                "name": "SLK library recovery, HistGBM joint composition versus ligand-only",
                "status": "computed",
                "reason": None,
                "kind": "recovery",
                "sources": screen["path"],
                "result": screen,
            }
        )

    payload = {
        "seed": SEED,
        "n_bootstrap": N_BOOT,
        "ci": "percentile 2.5 to 97.5",
        "models_refit": False,
        "test_sentence": (
            "Graph comparisons in this file remain a pairs bootstrap: the same resample enters both scores, "
            "the statistic is Spearman rho of the protein score minus Spearman rho of the ligand-only score, "
            "and the two-sided p-value for a difference of zero recenters each replicate by the observed difference. "
            "The HistGBM cold-protein comparison uses the same recentering, but the resample unit is the held-out kinase."
        ),
        "inventory_notes": notes,
        "prediction_files": files,
        "comparisons": comparisons,
    }
    json_path = ARTIFACTS / "metrics_uncertainty.json"
    md_path = ARTIFACTS / "metrics_uncertainty.md"
    json_path.write_text(json.dumps(payload, indent=2) + "\n")
    md_path.write_text(markdown_report(payload))
    hist = next(c for c in comparisons if c["name"] == HISTGBM_COLD_NAME)
    update_stat_review(hist)
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")
    print(f"updated {ARTIFACTS / 'STAT_REVIEW.md'}")
    for comp in comparisons:
        if comp["status"] != "computed":
            print(f"SKIP {comp['name']}")
            continue
        if comp["kind"] == "kinase_spearman":
            d = comp["result"]["delta_spearman_protein_minus_ligand"]
            e = comp["result"]["delta_ef_at_1pct_protein_minus_ligand"]
            print(
                f"KINASE DELTA {comp['name']}: {d['point']:.6f} "
                f"CI [{d['ci95']['low']:.6f}, {d['ci95']['high']:.6f}] "
                f"includes0={d['ci95']['includes_0']} p={d['test']['p_value']:.6f}"
            )
            print(
                f"KINASE EF1 {comp['name']}: {e['point']:.6f} "
                f"CI [{e['ci95']['low']:.6f}, {e['ci95']['high']:.6f}] "
                f"includes0={e['ci95']['includes_0']} p={e['test']['p_value']:.6f}"
            )
        elif comp["kind"] == "spearman":
            d = comp["result"]["delta_spearman_protein_minus_ligand"]
            print(
                f"DELTA {comp['name']}: {d['point']:.6f} "
                f"CI [{d['ci95']['low']:.6f}, {d['ci95']['high']:.6f}] "
                f"includes0={d['ci95']['includes_0']} p={d['test']['p_value']:.6f}"
            )
        else:
            for key, d in comp["result"]["deltas"].items():
                print(
                    f"SCREEN {key}: {d['point']:.6f} "
                    f"CI [{d['ci95']['low']:.6f}, {d['ci95']['high']:.6f}] "
                    f"includes0={d['ci95']['includes_0']} p={d['test']['p_value']:.6f}"
                )


if __name__ == "__main__":
    main()
