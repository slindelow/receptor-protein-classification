#!/usr/bin/env python
"""v1.1 protein-feature ablation: load existing splits, train/eval, export metrics_v1.1.*"""
from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from vs.config import load_config, resolve  # noqa: E402
from vs.features import (  # noqa: E402
    ESM2_EMBED_DIM,
    ESM2_MODEL_NAME,
    build_pair_matrix,
    protein_feat_dim,
    try_init_esm2,
)
from vs.metrics_lib import compute_metrics  # noqa: E402
from vs.model import save_artifact, train_model  # noqa: E402

TZ = ZoneInfo("America/Toronto")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_data_freeze(cfg: dict) -> dict:
    """Confirm raw DAVIS hashes match committed MANIFEST (do not rewrite MANIFEST)."""
    manifest_path = ROOT / "data" / "MANIFEST.json"
    manifest = json.loads(manifest_path.read_text())
    expected_data_freeze = cfg.get("data_freeze_id", "sofia-vs-v1-2026-10")
    if manifest.get("freeze_id") != expected_data_freeze:
        raise SystemExit(
            f"MANIFEST freeze_id={manifest.get('freeze_id')!r} != data_freeze_id={expected_data_freeze!r}"
        )
    raw_dir = resolve(cfg, "raw_dir")
    for key in ("davis_tdc.tab", "davis_pairs.csv"):
        meta = manifest["files"][key]
        path = ROOT / meta["path"]
        if not path.exists():
            # fall back to raw_dir layout
            path = raw_dir / Path(meta["path"]).name
        got = sha256_file(path)
        if got != meta["sha256"]:
            raise SystemExit(f"Hash mismatch for {key}: got {got} expected {meta['sha256']}")
    return manifest


def load_existing_split(pairs: pd.DataFrame, split_csv: Path) -> pd.DataFrame:
    """Join frozen pair labels onto existing split ID CSV. Never regenerates splits."""
    if not split_csv.exists():
        raise SystemExit(f"Missing split CSV (do not regenerate): {split_csv}")
    sid = pd.read_csv(split_csv)
    if "pair_id" not in sid.columns or "split" not in sid.columns:
        raise SystemExit(f"Split CSV missing pair_id/split: {split_csv}")
    keep = ["pair_id", "smiles", "sequence", "pkd", "drug_id", "target_id", "kd_nm"]
    keep = [c for c in keep if c in pairs.columns]
    merged = sid[["pair_id", "split"]].merge(pairs[keep], on="pair_id", how="inner")
    if len(merged) != len(sid):
        raise SystemExit(
            f"Split join size mismatch for {split_csv.name}: sid={len(sid)} merged={len(merged)}"
        )
    if set(merged["split"].unique()) - {"train", "val", "test"}:
        raise SystemExit(f"Unexpected split labels in {split_csv}")
    return merged


def fit_eval(
    df_split: pd.DataFrame,
    cfg: dict,
    include_protein: bool,
    protein_type: str = "aac",
    esm_embedder=None,
    esm_cache=None,
) -> tuple[dict, object]:
    radius = cfg["fingerprint"]["radius"]
    n_bits = cfg["fingerprint"]["n_bits"]
    thr = cfg["binder"]["pkd_threshold"]
    ef_fracs = cfg["eval"]["ef_fractions"]

    train = df_split[df_split["split"] == "train"]
    test = df_split[df_split["split"] == "test"]

    X_tr, m_tr = build_pair_matrix(
        train["smiles"],
        train["sequence"] if include_protein else None,
        radius,
        n_bits,
        include_protein,
        protein_type=protein_type,
        esm_embedder=esm_embedder,
        esm_cache=esm_cache,
    )
    y_tr = train["pkd"].to_numpy()[m_tr]
    X_tr = X_tr[m_tr]

    X_te, m_te = build_pair_matrix(
        test["smiles"],
        test["sequence"] if include_protein else None,
        radius,
        n_bits,
        include_protein,
        protein_type=protein_type,
        esm_embedder=esm_embedder,
        esm_cache=esm_cache,
    )
    y_te = test["pkd"].to_numpy()[m_te]
    X_te = X_te[m_te]

    model = train_model(X_tr, y_tr, cfg)
    pred = model.predict(X_te)
    metrics = compute_metrics(y_te, pred, thr, ef_fracs)
    metrics["n_train"] = float(len(y_tr))
    metrics["n_test"] = float(len(y_te))
    metrics["protein_type"] = protein_type if include_protein else "none"
    metrics["protein_dim"] = float(protein_feat_dim(protein_type) if include_protein else 0)
    metrics["feat_dim"] = float(X_tr.shape[1])
    return metrics, model


def _fmt(x) -> str:
    if x is None or (isinstance(x, float) and (np.isnan(x) or np.isinf(x))):
        return "nan"
    if isinstance(x, float):
        return f"{x:.4f}"
    return str(x)


def format_metrics_md(results: dict, cfg: dict) -> str:
    lines = [
        "# Virtual Screening v1.1 protein feature ablation",
        "",
        f"**Experiment freeze:** `{results['freeze_id']}`  ",
        f"**Data freeze:** `{results['data_freeze_id']}`  ",
        f"**Generated (America/Toronto):** {results['generated_at_america_toronto']}  ",
        f"**Binder threshold:** pKd ≥ {cfg['binder']['pkd_threshold']}  ",
        f"**Fingerprint:** ECFP radius={cfg['fingerprint']['radius']}, bits={cfg['fingerprint']['n_bits']}  ",
        f"**Model:** {cfg['model']['name']}  ",
        f"**Seed:** {cfg['seed']}  ",
        f"**ESM-2:** {results.get('esm2_status', 'n/a')}  ",
        "",
        "All numbers exported by `scripts/run_protein_ablation.py` / `reproduce_v1.1.sh`. Do not hand-edit.",
        "Splits loaded from existing `data/splits/*.csv` (not regenerated).",
        "",
        "## Results",
        "",
        "| Split | Protein mode | Model | Spearman ρ | EF@1% | EF@5% | AUROC | n_test | n_binders | feat_dim |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for split_name in ("scaffold", "cold_protein"):
        block = results["splits"][split_name]
        lo = block["ligand_only"]
        lines.append(
            "| {split} | n/a | ligand-only (ECFP) | {rho} | {ef1} | {ef5} | {auc} | {n} | {nb} | {fd} |".format(
                split=split_name,
                rho=_fmt(lo.get("spearman_rho")),
                ef1=_fmt(lo.get("ef_at_1pct")),
                ef5=_fmt(lo.get("ef_at_5pct")),
                auc=_fmt(lo.get("auroc")),
                n=int(lo.get("n_test", lo.get("n", 0))),
                nb=int(lo.get("n_binders", 0)),
                fd=int(lo.get("feat_dim", cfg["fingerprint"]["n_bits"])),
            )
        )
        for mode in results["protein_modes_ran"]:
            m = block["protein"][mode]
            label = f"ECFP+{mode.upper() if mode != 'aac_dpc' else 'AAC+DPC'}"
            if mode == "aac":
                label = "ECFP+AAC"
            elif mode == "dpc":
                label = "ECFP+DPC"
            elif mode == "aac_dpc":
                label = "ECFP+AAC+DPC"
            elif mode == "esm2":
                label = "ECFP+ESM2"
            lines.append(
                "| {split} | {mode} | {label} | {rho} | {ef1} | {ef5} | {auc} | {n} | {nb} | {fd} |".format(
                    split=split_name,
                    mode=mode,
                    label=label,
                    rho=_fmt(m.get("spearman_rho")),
                    ef1=_fmt(m.get("ef_at_1pct")),
                    ef5=_fmt(m.get("ef_at_5pct")),
                    auc=_fmt(m.get("auroc")),
                    n=int(m.get("n_test", m.get("n", 0))),
                    nb=int(m.get("n_binders", 0)),
                    fd=int(m.get("feat_dim", 0)),
                )
            )

    # Success criteria block (computed, not hand-edited)
    lines += ["", "## Success criteria (auto)", ""]
    cold = results["splits"]["cold_protein"]
    scaf = results["splits"]["scaffold"]
    lo_rho = cold["ligand_only"]["spearman_rho"]
    lo_ef1 = cold["ligand_only"]["ef_at_1pct"]
    lo_ef5 = cold["ligand_only"]["ef_at_5pct"]
    best_mode = None
    best_rho = -999.0
    for mode in results["protein_modes_ran"]:
        rho = cold["protein"][mode]["spearman_rho"]
        if rho > best_rho:
            best_rho = rho
            best_mode = mode
    delta = best_rho - lo_rho if best_mode is not None else float("nan")
    best_ef1 = cold["protein"][best_mode]["ef_at_1pct"] if best_mode else float("nan")
    best_ef5 = cold["protein"][best_mode]["ef_at_5pct"] if best_mode else float("nan")
    primary_rho = bool(delta >= 0.05) if best_mode is not None else False
    primary_ef = bool(
        (best_ef1 - lo_ef1) > 0.5 or (best_ef5 - lo_ef5) > 0.5
    ) if best_mode is not None else False
    lines.append(
        f"- Primary (cold_protein): best protein mode=`{best_mode}` ρ={_fmt(best_rho)} "
        f"vs ligand-only ρ={_fmt(lo_rho)} (Δ={_fmt(delta)}). "
        f"EF@1% {_fmt(best_ef1)} vs {_fmt(lo_ef1)}; EF@5% {_fmt(best_ef5)} vs {_fmt(lo_ef5)}."
    )
    lines.append(
        f"- Primary pass (ΔSpearman≥0.05 OR clear EF lift): "
        f"{'YES' if (primary_rho or primary_ef) else 'NO'} "
        f"(rho_gate={primary_rho}, ef_gate={primary_ef})."
    )
    # secondary: each full mode vs ligand-only on scaffold
    for mode in results["protein_modes_ran"]:
        s_rho = scaf["protein"][mode]["spearman_rho"]
        lo_s = scaf["ligand_only"]["spearman_rho"]
        ok = s_rho >= lo_s
        lines.append(
            f"- Secondary scaffold `{mode}` ≥ ligand-only: "
            f"{'YES' if ok else 'NO'} (ρ={_fmt(s_rho)} vs {_fmt(lo_s)})."
        )

    lines += [
        "",
        "## Notes",
        "",
        "- Same raw DAVIS file + MANIFEST hashes as MVP data freeze.",
        "- Same split CSVs under `data/splits/` (loaded, not rewritten).",
        "- Ligand-only control shared across protein modes (trained once per split).",
        "- HistGBM hyperparams and seeds unchanged from MVP.",
        f"- Reproduce float tolerance: `{cfg['eval']['metric_float_tol']}` (absolute).",
        "",
        "## Artifacts",
        "",
    ]
    for k, v in sorted(results.get("models", {}).items()):
        lines.append(f"- `{k}`: `{v}`")
    lines.append("")
    return "\n".join(lines)


def run(cfg_path: Path | None = None) -> dict:
    cfg_path = cfg_path or (ROOT / "configs" / "v1.1_protein.yaml")
    cfg = load_config(cfg_path)
    manifest = verify_data_freeze(cfg)

    raw_dir = resolve(cfg, "raw_dir")
    splits_dir = resolve(cfg, "splits_dir")
    models_dir = resolve(cfg, "models_dir")
    metrics_md = resolve(cfg, "metrics_md")
    metrics_json = resolve(cfg, "metrics_json")
    models_dir.mkdir(parents=True, exist_ok=True)

    pairs = pd.read_csv(raw_dir / "davis_pairs.csv")
    scaf = load_existing_split(pairs, splits_dir / "scaffold_split.csv")
    cold = load_existing_split(pairs, splits_dir / "cold_protein_split.csv")

    modes_wanted = list(cfg.get("protein_ablation", {}).get("modes", ["aac", "dpc", "aac_dpc", "esm2"]))
    esm_cfg = cfg.get("protein_features", {}).get("esm2", {})
    cache_dir = ROOT / esm_cfg.get("cache_dir", "data/cache/esm2_t6")

    esm_embedder = None
    esm_cache = None
    esm2_status = "not_requested"
    modes_ran = []
    cut_reasons = []

    if "esm2" in modes_wanted:
        print("Attempting ESM-2 init (CPU)...")
        esm_embedder, cut = try_init_esm2(cache_dir)
        if esm_embedder is None:
            esm2_status = f"CUT: {cut}"
            cut_reasons.append(esm2_status)
            print(esm2_status)
            modes_wanted = [m for m in modes_wanted if m != "esm2"]
        else:
            esm2_status = f"ran ({ESM2_MODEL_NAME}, dim={ESM2_EMBED_DIM}, cache={cache_dir})"
            print(esm2_status)
            # Precompute unique sequence embeddings once
            uniq_seqs = pd.concat([scaf["sequence"], cold["sequence"]]).unique().tolist()
            print(f"Caching ESM-2 embeddings for {len(uniq_seqs)} unique sequences...")
            esm_cache = esm_embedder.embed_many(uniq_seqs)
            print(f"ESM-2 cache entries: {len(list(cache_dir.glob('*.npy')))}")

    results: dict = {
        "freeze_id": cfg["freeze_id"],
        "data_freeze_id": cfg.get("data_freeze_id", manifest.get("freeze_id")),
        "generated_at_america_toronto": datetime.now(TZ).isoformat(timespec="seconds"),
        "config": {
            "seed": cfg["seed"],
            "fingerprint": cfg["fingerprint"],
            "binder_pkd_threshold": cfg["binder"]["pkd_threshold"],
            "model": cfg["model"]["name"],
            "histgbm": {k: cfg["model"][k] for k in cfg["model"] if k != "name"},
        },
        "esm2_status": esm2_status,
        "protein_modes_requested": list(cfg.get("protein_ablation", {}).get("modes", [])),
        "protein_modes_ran": [],
        "cuts": cut_reasons,
        "splits": {},
        "models": {},
        "manifest_sha256": {
            "davis_tdc.tab": manifest["files"]["davis_tdc.tab"]["sha256"],
            "davis_pairs.csv": manifest["files"]["davis_pairs.csv"]["sha256"],
        },
    }

    split_frames = (("scaffold", scaf), ("cold_protein", cold))

    for split_name, dsplit in split_frames:
        results["splits"][split_name] = {"protein": {}}
        print(f"=== {split_name}: ligand_only ===")
        metrics, model = fit_eval(dsplit, cfg, include_protein=False)
        results["splits"][split_name]["ligand_only"] = metrics
        tag = f"{split_name}_ligand_only"
        art_path = models_dir / f"hgb_{tag}.joblib"
        save_artifact(
            art_path,
            model,
            {
                "split": split_name,
                "mode": "ligand_only",
                "include_protein": False,
                "protein_type": None,
                "fingerprint": cfg["fingerprint"],
                "seed": cfg["seed"],
                "metrics": metrics,
                "freeze_id": cfg["freeze_id"],
            },
        )
        results["models"][tag] = str(art_path.relative_to(ROOT))

        for mode in modes_wanted:
            print(f"=== {split_name}: full + {mode} ===")
            metrics, model = fit_eval(
                dsplit,
                cfg,
                include_protein=True,
                protein_type=mode,
                esm_embedder=esm_embedder,
                esm_cache=esm_cache,
            )
            results["splits"][split_name]["protein"][mode] = metrics
            tag = f"{split_name}_full_{mode}"
            art_path = models_dir / f"hgb_{tag}.joblib"
            save_artifact(
                art_path,
                model,
                {
                    "split": split_name,
                    "mode": "full",
                    "include_protein": True,
                    "protein_type": mode,
                    "fingerprint": cfg["fingerprint"],
                    "seed": cfg["seed"],
                    "metrics": metrics,
                    "freeze_id": cfg["freeze_id"],
                },
            )
            results["models"][tag] = str(art_path.relative_to(ROOT))

    results["protein_modes_ran"] = list(modes_wanted)
    modes_ran = modes_wanted

    metrics_json.parent.mkdir(parents=True, exist_ok=True)
    with open(metrics_json, "w") as f:
        json.dump(results, f, indent=2, sort_keys=True)
        f.write("\n")

    md = format_metrics_md(results, cfg)
    with open(metrics_md, "w") as f:
        f.write(md)

    return results


def main() -> None:
    cfg_arg = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    results = run(cfg_arg)
    print("=== v1.1 metrics summary ===")
    for split, block in results["splits"].items():
        lo = block["ligand_only"]
        print(
            f"{split:14s} ligand_only   rho={lo['spearman_rho']:.4f}  "
            f"EF1={lo['ef_at_1pct']:.4f}  EF5={lo['ef_at_5pct']:.4f}"
        )
        for mode, m in block["protein"].items():
            print(
                f"{split:14s} full/{mode:8s} rho={m['spearman_rho']:.4f}  "
                f"EF1={m['ef_at_1pct']:.4f}  EF5={m['ef_at_5pct']:.4f}  "
                f"dim={int(m['feat_dim'])}"
            )
    print(f"metrics_v1.1.md  -> {ROOT / 'artifacts' / 'metrics_v1.1.md'}")
    print(f"metrics_v1.1.json-> {ROOT / 'artifacts' / 'metrics_v1.1.json'}")
    print(f"esm2_status: {results['esm2_status']}")


if __name__ == "__main__":
    main()
