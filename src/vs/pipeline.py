from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np
import pandas as pd

from .config import ROOT, load_config, resolve
from .data import freeze_davis, load_davis_tab
from .features import build_pair_matrix
from .metrics_lib import compute_metrics
from .model import load_artifact, save_artifact, train_model
from .splits import cold_protein_split, scaffold_split, write_split_ids

TZ = ZoneInfo("America/Toronto")


def _fit_eval_one(
    df_split: pd.DataFrame,
    cfg: dict,
    include_protein: bool,
    protein_type: str | None = None,
) -> tuple[dict, object]:
    radius = cfg["fingerprint"]["radius"]
    n_bits = cfg["fingerprint"]["n_bits"]
    thr = cfg["binder"]["pkd_threshold"]
    ef_fracs = cfg["eval"]["ef_fractions"]
    ptype = protein_type or cfg.get("protein_features", {}).get("type", "aac")

    train = df_split[df_split["split"] == "train"]
    test = df_split[df_split["split"] == "test"]
    # Keep train only for fit; HGB does internal validation_fraction split.

    X_tr, m_tr = build_pair_matrix(
        train["smiles"],
        train["sequence"] if include_protein else None,
        radius,
        n_bits,
        include_protein,
        protein_type=ptype,
    )
    y_tr = train["pkd"].to_numpy()[m_tr]
    X_tr = X_tr[m_tr]

    X_te, m_te = build_pair_matrix(
        test["smiles"],
        test["sequence"] if include_protein else None,
        radius,
        n_bits,
        include_protein,
        protein_type=ptype,
    )
    y_te = test["pkd"].to_numpy()[m_te]
    X_te = X_te[m_te]

    model = train_model(X_tr, y_tr, cfg)
    pred = model.predict(X_te)
    metrics = compute_metrics(y_te, pred, thr, ef_fracs)
    metrics["n_train"] = float(len(y_tr))
    metrics["n_test"] = float(len(y_te))
    return metrics, model


def run_all(cfg_path: Path | None = None) -> dict:
    cfg = load_config(cfg_path)
    raw_dir = resolve(cfg, "raw_dir")
    splits_dir = resolve(cfg, "splits_dir")
    models_dir = resolve(cfg, "models_dir")
    metrics_md = resolve(cfg, "metrics_md")
    metrics_json = resolve(cfg, "metrics_json")

    manifest_path = ROOT / "data" / "MANIFEST.json"
    freeze_davis(raw_dir, manifest_path)

    df = pd.read_csv(raw_dir / "davis_pairs.csv")
    seed = cfg["seed"]
    test_frac = cfg["splits"]["scaffold_test_frac"]
    val_frac = cfg["splits"]["val_frac"]
    cold_test = cfg["splits"]["cold_protein_test_frac"]

    scaf = scaffold_split(df, test_frac, val_frac, seed)
    cold = cold_protein_split(df, cold_test, val_frac, seed + 1)

    write_split_ids(scaf, splits_dir / "scaffold_split.csv", "murcko_scaffold")
    write_split_ids(cold, splits_dir / "cold_protein_split.csv", "cold_protein")

    results: dict = {
        "freeze_id": cfg["freeze_id"],
        "generated_at_america_toronto": datetime.now(TZ).isoformat(timespec="seconds"),
        "config": {
            "seed": seed,
            "fingerprint": cfg["fingerprint"],
            "binder_pkd_threshold": cfg["binder"]["pkd_threshold"],
            "model": cfg["model"]["name"],
        },
        "splits": {},
        "models": {},
    }

    # Full model (ECFP + AAC) on both splits; save scaffold full model as screening weights
    for split_name, dsplit in (("scaffold", scaf), ("cold_protein", cold)):
        results["splits"][split_name] = {}
        for mode, include_protein in (("full", True), ("ligand_only", False)):
            metrics, model = _fit_eval_one(dsplit, cfg, include_protein=include_protein)
            results["splits"][split_name][mode] = metrics
            tag = f"{split_name}_{mode}"
            art_path = models_dir / f"hgb_{tag}.joblib"
            meta = {
                "split": split_name,
                "mode": mode,
                "include_protein": include_protein,
                "fingerprint": cfg["fingerprint"],
                "seed": seed,
                "metrics": metrics,
            }
            save_artifact(art_path, model, meta)
            results["models"][tag] = str(art_path.relative_to(ROOT))

    # Canonical screening model: full features trained on all non-test scaffold rows
    # Retrain on train+val of scaffold for deployment
    deploy = scaf[scaf["split"] != "test"].copy()
    radius = cfg["fingerprint"]["radius"]
    n_bits = cfg["fingerprint"]["n_bits"]
    ptype = cfg.get("protein_features", {}).get("type", "aac")
    X_d, m_d = build_pair_matrix(
        deploy["smiles"], deploy["sequence"], radius, n_bits, True, protein_type=ptype
    )
    y_d = deploy["pkd"].to_numpy()[m_d]
    X_d = X_d[m_d]
    deploy_model = train_model(X_d, y_d, cfg)
    deploy_path = models_dir / "hgb_screen.joblib"
    save_artifact(
        deploy_path,
        deploy_model,
        {
            "split": "scaffold_train_val",
            "mode": "full",
            "include_protein": True,
            "fingerprint": cfg["fingerprint"],
            "seed": seed,
            "n_train": int(len(y_d)),
            "note": "Screening weights: ECFP+AAC HistGBM fit on scaffold train+val",
        },
    )
    results["models"]["screen"] = str(deploy_path.relative_to(ROOT))

    # Export metrics JSON + MD (no hand-edited numbers)
    metrics_json.parent.mkdir(parents=True, exist_ok=True)
    with open(metrics_json, "w") as f:
        json.dump(results, f, indent=2, sort_keys=True)
        f.write("\n")

    md = _format_metrics_md(results, cfg)
    with open(metrics_md, "w") as f:
        f.write(md)

    return results


def _fmt(x) -> str:
    if x is None or (isinstance(x, float) and (np.isnan(x) or np.isinf(x))):
        return "nan"
    if isinstance(x, float):
        return f"{x:.4f}"
    return str(x)


def _format_metrics_md(results: dict, cfg: dict) -> str:
    lines = [
        "# Virtual Screening MVP — Metrics",
        "",
        f"**Freeze:** `{results['freeze_id']}`  ",
        f"**Generated (America/Toronto):** {results['generated_at_america_toronto']}  ",
        f"**Binder threshold:** pKd ≥ {cfg['binder']['pkd_threshold']}  ",
        f"**Fingerprint:** ECFP radius={cfg['fingerprint']['radius']}, bits={cfg['fingerprint']['n_bits']}  ",
        f"**Model:** {cfg['model']['name']}  ",
        f"**Seed:** {cfg['seed']}  ",
        "",
        "All numbers exported by `scripts/run_train_eval.py` / `reproduce.sh`. Do not hand-edit.",
        "",
        "## Results",
        "",
        "| Split | Model | Spearman ρ | EF@1% | EF@5% | AUROC | n_test | n_binders |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for split_name in ("scaffold", "cold_protein"):
        for mode in ("full", "ligand_only"):
            m = results["splits"][split_name][mode]
            label = "ECFP+AAC" if mode == "full" else "ligand-only (ECFP)"
            lines.append(
                "| {split} | {label} | {rho} | {ef1} | {ef5} | {auc} | {n} | {nb} |".format(
                    split=split_name,
                    label=label,
                    rho=_fmt(m.get("spearman_rho")),
                    ef1=_fmt(m.get("ef_at_1pct")),
                    ef5=_fmt(m.get("ef_at_5pct")),
                    auc=_fmt(m.get("auroc")),
                    n=int(m.get("n_test", m.get("n", 0))),
                    nb=int(m.get("n_binders", 0)),
                )
            )
    lines += [
        "",
        "## Notes",
        "",
        "- **Scaffold split:** Murcko scaffolds held out (test ≈ 20%, val ≈ 10%).",
        "- **Cold-protein split:** entire target proteins held out (test ≈ 20%, val ≈ 10%).",
        "- **Ligand-only ablation:** same splits, protein AAC features removed — control for ligand memorization.",
        f"- Reproduce float tolerance: `{cfg['eval']['metric_float_tol']}` (absolute) on Spearman/AUROC/EF.",
        "",
        "## Artifacts",
        "",
    ]
    for k, v in sorted(results["models"].items()):
        lines.append(f"- `{k}`: `{v}`")
    lines.append("")
    return "\n".join(lines)


def screen_library(
    target_seq: str,
    library_csv: Path,
    model_path: Path | None = None,
    cfg_path: Path | None = None,
    smiles_col: str = "smiles",
    id_col: str | None = "mol_id",
) -> pd.DataFrame:
    cfg = load_config(cfg_path)
    if model_path is None:
        model_path = resolve(cfg, "models_dir") / "hgb_screen.joblib"
    art = load_artifact(model_path)
    model = art["model"]
    meta = art["meta"]
    radius = meta["fingerprint"]["radius"]
    n_bits = meta["fingerprint"]["n_bits"]
    include_protein = bool(meta.get("include_protein", True))

    lib = pd.read_csv(library_csv)
    if smiles_col not in lib.columns:
        raise ValueError(f"Library missing column '{smiles_col}'")
    if id_col is None or id_col not in lib.columns:
        lib = lib.copy()
        lib["mol_id"] = [f"mol_{i}" for i in range(len(lib))]
        id_col = "mol_id"

    seqs = [target_seq] * len(lib) if include_protein else None
    ptype = meta.get("protein_type") or cfg.get("protein_features", {}).get("type", "aac")
    X, mask = build_pair_matrix(
        lib[smiles_col], seqs, radius, n_bits, include_protein, protein_type=ptype
    )
    scores = np.full(len(lib), np.nan, dtype=float)
    if mask.any():
        scores[mask] = model.predict(X[mask])
    out = pd.DataFrame(
        {
            "mol_id": lib[id_col],
            "smiles": lib[smiles_col],
            "score_pkd": scores,
            "valid_smiles": mask,
        }
    )
    out = out.sort_values("score_pkd", ascending=False, na_position="last").reset_index(drop=True)
    out.insert(0, "rank", np.arange(1, len(out) + 1))
    return out
