#!/usr/bin/env python
"""ECE / reliability calibration on best VS models (v1.1 HistGBM + v1.2 Chemprop).

Uses frozen DAVIS splits. Binary labels: pKd >= binder threshold.
Exports artifacts/metrics_calibration.md (+ .json). Does not retrain models unless
predictions are missing (then loads saved preds or re-infers from checkpoints/joblibs).
"""
from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from vs.config import load_config, resolve  # noqa: E402
from vs.features import build_pair_matrix, try_init_esm2  # noqa: E402
from vs.metrics_lib import compute_metrics  # noqa: E402
from vs.model import load_artifact  # noqa: E402

TZ = ZoneInfo("America/Toronto")


def expected_calibration_error(y_true: np.ndarray, y_prob: np.ndarray, n_bins: int = 10) -> float:
    y_true = np.asarray(y_true).astype(int)
    y_prob = np.clip(np.asarray(y_prob, dtype=float), 0.0, 1.0)
    bins = np.linspace(0.0, 1.0, n_bins + 1)
    ece = 0.0
    n = len(y_true)
    if n == 0:
        return float("nan")
    for i in range(n_bins):
        lo, hi = bins[i], bins[i + 1]
        if i == n_bins - 1:
            mask = (y_prob >= lo) & (y_prob <= hi)
        else:
            mask = (y_prob >= lo) & (y_prob < hi)
        if not np.any(mask):
            continue
        conf = float(y_prob[mask].mean())
        acc = float(y_true[mask].mean())
        ece += (mask.sum() / n) * abs(acc - conf)
    return float(ece)


def brier_score(y_true: np.ndarray, y_prob: np.ndarray) -> float:
    y_true = np.asarray(y_true, dtype=float)
    y_prob = np.asarray(y_prob, dtype=float)
    return float(np.mean((y_prob - y_true) ** 2))


def scores_to_probs_val_calibrated(
    y_val: np.ndarray,
    s_val: np.ndarray,
    s_test: np.ndarray,
    method: str = "isotonic",
) -> tuple[np.ndarray, np.ndarray]:
    """Map continuous scores to [0,1] probs via val-set calibration."""
    thr_scores = None  # unused; we calibrate score->P(binder) using val labels
    y_bin_val = y_val  # already binary when passed
    if method == "platt":
        lr = LogisticRegression(max_iter=1000)
        lr.fit(s_val.reshape(-1, 1), y_bin_val)
        p_val = lr.predict_proba(s_val.reshape(-1, 1))[:, 1]
        p_test = lr.predict_proba(s_test.reshape(-1, 1))[:, 1]
        return p_val, p_test
    # isotonic
    iso = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0)
    iso.fit(s_val, y_bin_val)
    return iso.predict(s_val), iso.predict(s_test)


def load_split(pairs: pd.DataFrame, split_csv: Path) -> pd.DataFrame:
    sid = pd.read_csv(split_csv)
    keep = [c for c in ["pair_id", "smiles", "sequence", "pkd"] if c in pairs.columns]
    return sid[["pair_id", "split"]].merge(pairs[keep], on="pair_id", how="inner")


def histgbm_predict(model_path: Path, df: pd.DataFrame, protein_type: str | None, cfg: dict, esm=None, esm_cache=None):
    art = load_artifact(model_path)
    model = art["model"] if isinstance(art, dict) and "model" in art else art
    include = protein_type is not None and protein_type != "none"
    X, mask = build_pair_matrix(
        df["smiles"],
        df["sequence"] if include else None,
        cfg["fingerprint"]["radius"],
        cfg["fingerprint"]["n_bits"],
        include,
        protein_type=protein_type or "aac",
        esm_embedder=esm,
        esm_cache=esm_cache,
    )
    pred = np.full(len(df), np.nan)
    pred[mask] = model.predict(X[mask])
    return pred, mask


def main() -> None:
    cfg = load_config(ROOT / "configs" / "v1.1_protein.yaml")
    thr = float(cfg["binder"]["pkd_threshold"])
    pairs = pd.read_csv(ROOT / "data" / "raw" / "davis_pairs.csv")
    splits_dir = ROOT / "data" / "splits"

    # Models to calibrate: v1.1 winners + chemprop if preds exist
    targets = []
    v11 = ROOT / "artifacts" / "models_v1.1"
    for split in ("scaffold", "cold_protein"):
        targets.append({"split": split, "name": "hgb_ligand_only", "kind": "hgb", "path": v11 / f"hgb_{split}_ligand_only.joblib", "protein": None})
        targets.append({"split": split, "name": "hgb_aac_dpc", "kind": "hgb", "path": v11 / f"hgb_{split}_full_aac_dpc.joblib", "protein": "aac_dpc"})
        targets.append({"split": split, "name": "hgb_esm2", "kind": "hgb", "path": v11 / f"hgb_{split}_full_esm2.joblib", "protein": "esm2"})

    v12 = ROOT / "artifacts" / "models_v1.2"
    for split in ("scaffold", "cold_protein"):
        for tag, pname in (("ligand_only", None), ("protein", "aac_dpc")):
            pred_path = v12 / f"preds_{split}_{tag}.npz"
            targets.append({
                "split": split,
                "name": f"chemprop_{tag}" + (f"_{pname}" if pname else ""),
                "kind": "chemprop_preds",
                "path": pred_path,
                "protein": pname,
            })

    esm_embedder, esm_cache = None, None
    if any(t.get("protein") == "esm2" for t in targets):
        emb, cut = try_init_esm2(ROOT / "data" / "cache" / "esm2_t6")
        if emb is not None:
            esm_embedder = emb
            # cache will be filled on demand

    rows = []
    for t in targets:
        split_csv = splits_dir / f"{t['split']}_split.csv"
        if not t["path"].exists():
            rows.append({"split": t["split"], "model": t["name"], "status": f"MISSING {t['path'].name}"})
            continue
        df = load_split(pairs, split_csv)
        val = df[df["split"] == "val"]
        test = df[df["split"] == "test"]
        y_val = (val["pkd"].to_numpy() >= thr).astype(int)
        y_test = (test["pkd"].to_numpy() >= thr).astype(int)

        if t["kind"] == "chemprop_preds":
            z = np.load(t["path"])
            # preds were saved for test only; re-calibrate needs val scores.
            # For chemprop we only have test preds in npz; use rank-based sigmoid on test as uncalibrated,
            # and fit isotonic on a holdout 20% of test for honesty when val preds absent.
            s_test = z["y_pred"]
            y_true_test = z["y_true"]
            y_bin = (y_true_test >= thr).astype(int)
            # Split test into cal/eval 30/70 with fixed seed for reporting when val scores unavailable
            rng = np.random.RandomState(42)
            idx = rng.permutation(len(s_test))
            n_cal = max(50, int(0.3 * len(s_test)))
            cal_i, ev_i = idx[:n_cal], idx[n_cal:]
            p_cal, p_ev = scores_to_probs_val_calibrated(
                y_bin[cal_i], s_test[cal_i], s_test[ev_i], method="isotonic"
            )
            # Uncalibrated: min-max to [0,1] on eval portion using cal stats
            lo, hi = float(s_test[cal_i].min()), float(s_test[cal_i].max())
            p_unc = np.clip((s_test[ev_i] - lo) / (hi - lo + 1e-12), 0, 1)
            y_ev = y_bin[ev_i]
            ece_unc = expected_calibration_error(y_ev, p_unc)
            ece_cal = expected_calibration_error(y_ev, p_ev)
            brier_unc = brier_score(y_ev, p_unc)
            brier_cal = brier_score(y_ev, p_ev)
            ranking = compute_metrics(y_true_test, s_test, thr, [0.01, 0.05])
            rows.append({
                "split": t["split"], "model": t["name"], "status": "ok",
                "cal_protocol": "isotonic_on_30pct_test_holdout",
                "n_eval": int(len(y_ev)), "n_binders_eval": int(y_ev.sum()),
                "ece_uncalibrated": ece_unc, "ece_isotonic": ece_cal,
                "brier_uncalibrated": brier_unc, "brier_isotonic": brier_cal,
                "spearman_rho": ranking["spearman_rho"],
                "ef_at_1pct": ranking["ef_at_1pct"], "auroc": ranking["auroc"],
            })
            continue

        # HistGBM: predict val+test
        s_val, m_val = histgbm_predict(t["path"], val, t["protein"], cfg, esm_embedder, esm_cache)
        s_test, m_test = histgbm_predict(t["path"], test, t["protein"], cfg, esm_embedder, esm_cache)
        # drop invalid
        s_val, y_val = s_val[m_val], y_val[m_val]
        s_test, y_test_f = s_test[m_test], y_test[m_test]
        y_pkd_test = test["pkd"].to_numpy()[m_test]

        p_val, p_test = scores_to_probs_val_calibrated(y_val, s_val, s_test, method="isotonic")
        lo, hi = float(s_val.min()), float(s_val.max())
        p_unc = np.clip((s_test - lo) / (hi - lo + 1e-12), 0, 1)
        ece_unc = expected_calibration_error(y_test_f, p_unc)
        ece_cal = expected_calibration_error(y_test_f, p_test)
        ranking = compute_metrics(y_pkd_test, s_test, thr, [0.01, 0.05])
        rows.append({
            "split": t["split"], "model": t["name"], "status": "ok",
            "cal_protocol": "isotonic_fit_on_val",
            "n_eval": int(len(y_test_f)), "n_binders_eval": int(y_test_f.sum()),
            "ece_uncalibrated": ece_unc, "ece_isotonic": ece_cal,
            "brier_uncalibrated": brier_score(y_test_f, p_unc),
            "brier_isotonic": brier_score(y_test_f, p_test),
            "spearman_rho": ranking["spearman_rho"],
            "ef_at_1pct": ranking["ef_at_1pct"], "auroc": ranking["auroc"],
        })

    out = {
        "freeze_id": "sofia-vs-calibration",
        "data_freeze_id": "sofia-vs-v1-2026-10",
        "generated_at_america_toronto": datetime.now(TZ).isoformat(timespec="seconds"),
        "binder_pkd_threshold": thr,
        "n_bins_ece": 10,
        "rows": rows,
    }
    md_path = ROOT / "artifacts" / "metrics_calibration.md"
    js_path = ROOT / "artifacts" / "metrics_calibration.json"
    with open(js_path, "w") as f:
        json.dump(out, f, indent=2, sort_keys=True)
        f.write("\n")

    def fmt(x):
        if isinstance(x, float):
            return f"{x:.4f}"
        return str(x)

    lines = [
        "# Calibration / ECE (Sofia VS)",
        "",
        f"**Generated (America/Toronto):** {out['generated_at_america_toronto']}  ",
        f"**Data freeze:** `{out['data_freeze_id']}`  ",
        f"**Binder threshold:** pKd >= {thr}  ",
        "**ECE bins:** 10  ",
        "",
        "Continuous affinity scores mapped to P(binder) via isotonic regression.",
        "HistGBM: fit on val split. Chemprop: val scores not saved; isotonic on 30% test holdout (documented).",
        "",
        "| Split | Model | Protocol | ECE unc | ECE iso | Brier unc | Brier iso | AUROC | n_eval |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        if r.get("status") != "ok":
            lines.append(f"| {r['split']} | {r['model']} | {r.get('status','?')} |  |  |  |  |  |  |")
            continue
        lines.append(
            "| {split} | {model} | {proto} | {eu} | {ei} | {bu} | {bi} | {auc} | {n} |".format(
                split=r["split"], model=r["model"], proto=r["cal_protocol"],
                eu=fmt(r["ece_uncalibrated"]), ei=fmt(r["ece_isotonic"]),
                bu=fmt(r["brier_uncalibrated"]), bi=fmt(r["brier_isotonic"]),
                auc=fmt(r["auroc"]), n=r["n_eval"],
            )
        )
    lines += ["", "## Notes", "",
              "- Ranking metrics (Spearman/EF) unchanged by calibration; ECE/Brier measure probability quality.",
              "- Honest negative: if ECE stays high after isotonic, scores are poorly ordered for probability.",
              ""]
    md_path.write_text("\n".join(lines))
    print(f"Wrote {md_path}")
    print(f"Wrote {js_path}")
    for r in rows:
        if r.get("status") == "ok":
            print(f"{r['split']:14s} {r['model']:28s} ECE {r['ece_uncalibrated']:.4f}->{r['ece_isotonic']:.4f}")
        else:
            print(f"{r['split']:14s} {r['model']:28s} {r['status']}")


if __name__ == "__main__":
    main()
