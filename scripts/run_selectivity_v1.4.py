#!/usr/bin/env python
"""v1.4 selectivity sketch on one DAVIS kinase pair (cold-pair holdout).

For ligands with measured pKd on BOTH targets:
  Δtrue = pKd_A - pKd_B
  Δpred = score(lig, A) - score(lig, B)   # HistGBM+aac_dpc vs ligand-only
Report Spearman(Δpred, Δtrue) and top-k enrichment for |Δtrue| ≥ threshold.
"""
from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from vs.config import load_config, resolve  # noqa: E402
from vs.features import build_pair_matrix, sequence_to_aac_dpc, smiles_to_ecfp  # noqa: E402
from vs.metrics_lib import enrichment_factor  # noqa: E402
from vs.model import save_artifact, train_model  # noqa: E402

TZ = ZoneInfo("America/Toronto")
AA = set("ACDEFGHIKLMNPQRSTVWY")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_data_freeze(cfg: dict) -> dict:
    manifest_path = ROOT / "data" / "MANIFEST.json"
    manifest = json.loads(manifest_path.read_text())
    expected = cfg.get("data_freeze_id", "sofia-vs-v1-2026-10")
    if manifest.get("freeze_id") != expected:
        raise SystemExit(f"MANIFEST freeze_id mismatch: {manifest.get('freeze_id')!r} vs {expected!r}")
    for key in ("davis_tdc.tab", "davis_pairs.csv"):
        meta = manifest["files"][key]
        path = ROOT / meta["path"]
        if not path.exists():
            path = resolve(cfg, "raw_dir") / Path(meta["path"]).name
        got = sha256_file(path)
        if got != meta["sha256"]:
            raise SystemExit(f"Hash mismatch for {key}")
    return manifest


def clean_seq(s: str) -> str:
    return "".join(c for c in (s or "").upper() if c in AA)


def best_window_identity(a: str, b: str) -> float:
    a, b = clean_seq(a), clean_seq(b)
    if not a or not b:
        return 0.0
    if len(a) > len(b):
        a, b = b, a
    n, m = len(a), len(b)
    best = 0
    for start in range(m - n + 1):
        matches = sum(1 for i in range(n) if a[i] == b[start + i])
        best = max(best, matches)
    return best / n


def aac_dpc_cosine(a: str, b: str) -> float:
    va, vb = sequence_to_aac_dpc(a), sequence_to_aac_dpc(b)
    return float(np.dot(va, vb) / (np.linalg.norm(va) * np.linalg.norm(vb) + 1e-12))


def fit_model(train: pd.DataFrame, cfg: dict, include_protein: bool, protein_type: str):
    radius = cfg["fingerprint"]["radius"]
    n_bits = cfg["fingerprint"]["n_bits"]
    X, m = build_pair_matrix(
        train["smiles"],
        train["sequence"] if include_protein else None,
        radius,
        n_bits,
        include_protein=include_protein,
        protein_type=protein_type,
    )
    y = train["pkd"].to_numpy()[m]
    X = X[m]
    print(f"  train n={len(y)} include_protein={include_protein} dim={X.shape[1]}", flush=True)
    model = train_model(X, y, cfg)
    return model


def predict_pairs(model, smiles: list[str], sequences: list[str], cfg: dict, include_protein: bool, protein_type: str):
    radius = cfg["fingerprint"]["radius"]
    n_bits = cfg["fingerprint"]["n_bits"]
    X, m = build_pair_matrix(
        smiles,
        sequences if include_protein else None,
        radius,
        n_bits,
        include_protein=include_protein,
        protein_type=protein_type,
    )
    preds = np.full(len(smiles), np.nan, dtype=float)
    if m.any():
        preds[m] = model.predict(X[m])
    return preds


def main():
    cfg_path = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "configs" / "v1.4_selectivity.yaml"
    cfg = load_config(cfg_path)
    manifest = verify_data_freeze(cfg)
    sel = cfg["selectivity"]
    ta, tb = sel["target_a"], sel["target_b"]
    dthr = float(sel["delta_pkd_threshold"])
    protein_type = sel.get("protein_type", "aac_dpc")
    ef_fracs = cfg["eval"]["ef_fractions"]

    pairs_path = resolve(cfg, "raw_dir") / "davis_pairs.csv"
    df = pd.read_csv(pairs_path)
    print(f"Loaded DAVIS n={len(df)} targets={df['target_id'].nunique()} ligands={df['smiles'].nunique()}")

    if ta not in set(df["target_id"]) or tb not in set(df["target_id"]):
        raise SystemExit(f"Targets {ta}/{tb} not in DAVIS")

    seq_a = df.loc[df["target_id"] == ta, "sequence"].iloc[0]
    seq_b = df.loc[df["target_id"] == tb, "sequence"].iloc[0]
    win_id = best_window_identity(seq_a, seq_b)
    cos = aac_dpc_cosine(seq_a, seq_b)

    # Dual-measured ligands
    piv = df.pivot_table(index="smiles", columns="target_id", values="pkd", aggfunc="mean")
    both = piv[[ta, tb]].dropna().copy()
    both["delta_true"] = both[ta] - both[tb]
    both["abs_delta_true"] = both["delta_true"].abs()
    both["selective"] = (both["abs_delta_true"] >= dthr).astype(int)
    n_lig = len(both)
    n_sel = int(both["selective"].sum())
    print(f"Pair {ta}/{tb}: n_dual={n_lig} n_selective(|Δ|≥{dthr})={n_sel}")
    print(f"  window_seq_id={win_id:.4f} aac_dpc_cos={cos:.4f}")

    # Cold-pair holdout: train on all pairs except target A or B
    train = df[~df["target_id"].isin([ta, tb])].copy()
    print(f"Cold-pair train n={len(train)} (excluded {ta},{tb})")

    models_dir = resolve(cfg, "models_dir")
    models_dir.mkdir(parents=True, exist_ok=True)

    print("Training HistGBM + aac_dpc ...", flush=True)
    model_full = fit_model(train, cfg, include_protein=True, protein_type=protein_type)
    save_artifact(
        models_dir / f"hgb_selectivity_{ta}_{tb}_aac_dpc.joblib",
        model_full,
        {"include_protein": True, "protein_type": protein_type, "pair": [ta, tb], "protocol": "cold_pair_holdout"},
    )

    print("Training HistGBM ligand-only ...", flush=True)
    model_lig = fit_model(train, cfg, include_protein=False, protein_type=protein_type)
    save_artifact(
        models_dir / f"hgb_selectivity_{ta}_{tb}_ligand_only.joblib",
        model_lig,
        {"include_protein": False, "pair": [ta, tb], "protocol": "cold_pair_holdout"},
    )

    smiles = both.index.astype(str).tolist()
    seqs_a = [seq_a] * len(smiles)
    seqs_b = [seq_b] * len(smiles)

    def eval_model(model, include_protein: bool, name: str) -> dict:
        pred_a = predict_pairs(model, smiles, seqs_a, cfg, include_protein, protein_type)
        pred_b = predict_pairs(model, smiles, seqs_b, cfg, include_protein, protein_type)
        dpred = pred_a - pred_b
        dtrue = both["delta_true"].to_numpy()
        abs_dpred = np.abs(dpred)
        y_sel = both["selective"].to_numpy().astype(int)

        if np.allclose(dpred, dpred[0], equal_nan=False) or np.isnan(dpred).all():
            rho = float("nan")
            rho_note = "constant_or_nan_delta_pred"
        else:
            rho, _ = spearmanr(dtrue, dpred, nan_policy="omit")
            rho = float(rho)
            rho_note = "ok"

        # Absolute selectivity enrichment: rank by |Δpred|
        efs = {}
        for frac in ef_fracs:
            efs[f"ef_at_{int(frac * 100)}pct"] = enrichment_factor(y_sel, abs_dpred, frac)

        # Directed: SRC-preferring (Δtrue ≥ +thr) ranked by Δpred
        y_a_pref = (dtrue >= dthr).astype(int)
        efs_dir = {}
        for frac in ef_fracs:
            efs_dir[f"ef_A_pref_at_{int(frac * 100)}pct"] = enrichment_factor(y_a_pref, dpred, frac)

        out = {
            "model": name,
            "n_ligands": int(len(smiles)),
            "n_selective_abs": int(y_sel.sum()),
            "n_A_preferring": int(y_a_pref.sum()),
            "spearman_delta": rho,
            "spearman_note": rho_note,
            "mean_abs_delta_pred": float(np.nanmean(abs_dpred)),
            "mean_abs_delta_true": float(np.mean(np.abs(dtrue))),
            "std_delta_pred": float(np.nanstd(dpred)),
            "enrichment_abs_selective": efs,
            "enrichment_A_preferring": efs_dir,
            "delta_pred_constant": bool(np.allclose(dpred, dpred[0], equal_nan=False)),
        }
        # stash arrays for export sample
        out["_arrays"] = {
            "smiles": smiles,
            "pkd_a": both[ta].to_numpy().tolist(),
            "pkd_b": both[tb].to_numpy().tolist(),
            "delta_true": dtrue.tolist(),
            "pred_a": pred_a.tolist(),
            "pred_b": pred_b.tolist(),
            "delta_pred": dpred.tolist(),
            "selective": y_sel.tolist(),
        }
        return out

    res_full = eval_model(model_full, True, "HistGBM+aac_dpc")
    res_lig = eval_model(model_lig, False, "HistGBM ligand-only")

    now = datetime.now(TZ).strftime("%Y-%m-%d %H:%M %Z")
    primary_pass = (
        not np.isnan(res_full["spearman_delta"])
        and res_full["spearman_delta"] > 0.1
        and (
            res_full["enrichment_abs_selective"].get("ef_at_20pct", 0) > 1.2
            or res_full["enrichment_abs_selective"].get("ef_at_10pct", 0) > 1.2
        )
    )
    # Soft sketch gate: positive Spearman OR any EF>1.2; honest null OK
    soft_signal = (
        (not np.isnan(res_full["spearman_delta"]) and res_full["spearman_delta"] > 0.05)
        or any(v > 1.1 for v in res_full["enrichment_abs_selective"].values() if not np.isnan(v))
    )

    payload = {
        "freeze_id": cfg["freeze_id"],
        "data_freeze_id": cfg.get("data_freeze_id"),
        "timestamp_america_toronto": now,
        "pair": {
            "target_a": ta,
            "target_b": tb,
            "rationale": sel.get("rationale"),
            "window_sequence_identity": win_id,
            "aac_dpc_cosine": cos,
            "len_a": len(clean_seq(seq_a)),
            "len_b": len(clean_seq(seq_b)),
            "n_dual_ligands": n_lig,
            "delta_pkd_threshold": dthr,
            "n_selective_abs": n_sel,
        },
        "protocol": {
            "name": "cold_pair_holdout",
            "train_n": int(len(train)),
            "protein_type": protein_type,
            "scorer": "HistGradientBoostingRegressor + ECFP(r=2,2048)",
            "delta_definition": f"Δ = score({ta}) - score({tb})",
            "selective_definition": f"|ΔpKd| >= {dthr}",
        },
        "results": {
            "aac_dpc": {k: v for k, v in res_full.items() if not k.startswith("_")},
            "ligand_only": {k: v for k, v in res_lig.items() if not k.startswith("_")},
        },
        "success": {
            "primary_sketch_gate": bool(primary_pass),
            "soft_signal": bool(soft_signal),
            "note": "Sketch gate: Spearman(Δ)>0.1 AND EF@10/20%>1.2. Soft: Spearman>0.05 OR any EF>1.1. Honest null OK.",
        },
        "data_manifest_freeze_id": manifest.get("freeze_id"),
        "ligand_table_aac_dpc": res_full["_arrays"],
    }

    metrics_json = resolve(cfg, "metrics_json")
    metrics_md = resolve(cfg, "metrics_md")
    metrics_json.parent.mkdir(parents=True, exist_ok=True)
    # JSON without huge duplication risk — keep ligand table
    metrics_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    def fmt_ef(d):
        return ", ".join(f"{k}={v:.4f}" if not (isinstance(v, float) and np.isnan(v)) else f"{k}=nan" for k, v in d.items())

    rho_f = res_full["spearman_delta"]
    rho_l = res_lig["spearman_delta"]
    md = f"""# Selectivity metrics — sofia-vs-v1.4-selectivity

**Timestamp (America/Toronto):** {now}  
**Data freeze:** `{cfg.get("data_freeze_id")}` (MANIFEST verified)  
**Experiment freeze:** `{cfg["freeze_id"]}`

## Kinase pair

| Field | Value |
|---|---|
| Pair | **{ta} / {tb}** |
| Why | {sel.get("rationale")} |
| Best-window sequence identity | {win_id:.4f} |
| AAC+DPC cosine | {cos:.4f} |
| Dual-measured ligands | {n_lig} |
| Selective (\\|ΔpKd\\| ≥ {dthr}) | {n_sel} |
| Protocol | cold-pair holdout (train excludes {ta} and {tb}) |
| Scorer | HistGBM + ECFP + **aac_dpc** (v1.1 cold winner) |

Δ := score({ta}) − score({tb})  (same for true pKd).

## Results

| Model | Spearman(Δpred, Δtrue) | mean\\|Δpred\\| | EF abs-selective | EF {ta}-pref |
|---|---:|---:|---|---|
| HistGBM+aac_dpc | {rho_f if not (isinstance(rho_f, float) and np.isnan(rho_f)) else "nan"} | {res_full["mean_abs_delta_pred"]:.4f} | {fmt_ef(res_full["enrichment_abs_selective"])} | {fmt_ef(res_full["enrichment_A_preferring"])} |
| HistGBM ligand-only | {rho_l if not (isinstance(rho_l, float) and np.isnan(rho_l)) else "nan"} | {res_lig["mean_abs_delta_pred"]:.4f} | {fmt_ef(res_lig["enrichment_abs_selective"])} | {fmt_ef(res_lig["enrichment_A_preferring"])} |

Ligand-only Δpred is (near-)constant because scores ignore protein → cannot rank selectivity (honest null).

## Success (sketch)

- Primary gate (ρ>0.1 and EF@10/20%>1.2): **{"YES" if primary_pass else "NO"}**
- Soft signal: **{"YES" if soft_signal else "NO"}**
- Honest null is an acceptable outcome for this P2 sketch.

## Takeaway

See `artifacts/DECISIONS.md` / `LAB_NOTEBOOK.md` (v1.4 section).
"""
    metrics_md.write_text(md)
    print(md)
    print(f"Wrote {metrics_md} and {metrics_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
