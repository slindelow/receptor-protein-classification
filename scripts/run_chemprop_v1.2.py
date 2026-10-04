#!/usr/bin/env python
"""v1.2 Chemprop D-MPNN: ligand-only vs Chemprop+aac_dpc on frozen DAVIS splits."""
from __future__ import annotations

import hashlib
import json
import sys
import warnings
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np
import pandas as pd
import torch
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from vs.config import load_config, resolve  # noqa: E402
from vs.features import sequence_to_aac_dpc  # noqa: E402
from vs.metrics_lib import compute_metrics  # noqa: E402

TZ = ZoneInfo("America/Toronto")
warnings.filterwarnings("ignore", category=UserWarning)


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
        raise SystemExit(
            f"MANIFEST freeze_id={manifest.get('freeze_id')!r} != data_freeze_id={expected!r}"
        )
    for key in ("davis_tdc.tab", "davis_pairs.csv"):
        meta = manifest["files"][key]
        path = ROOT / meta["path"]
        if not path.exists():
            path = resolve(cfg, "raw_dir") / Path(meta["path"]).name
        got = sha256_file(path)
        if got != meta["sha256"]:
            raise SystemExit(f"Hash mismatch for {key}: got {got} expected {meta['sha256']}")
    return manifest


def load_existing_split(pairs: pd.DataFrame, split_csv: Path) -> pd.DataFrame:
    if not split_csv.exists():
        raise SystemExit(f"Missing split CSV (do not regenerate): {split_csv}")
    sid = pd.read_csv(split_csv)
    keep = ["pair_id", "smiles", "sequence", "pkd", "drug_id", "target_id", "kd_nm"]
    keep = [c for c in keep if c in pairs.columns]
    merged = sid[["pair_id", "split"]].merge(pairs[keep], on="pair_id", how="inner")
    if len(merged) != len(sid):
        raise SystemExit(
            f"Split join size mismatch for {split_csv.name}: sid={len(sid)} merged={len(merged)}"
        )
    return merged


def set_seed(seed: int) -> None:
    np.random.seed(seed)
    torch.manual_seed(seed)
    try:
        import lightning.pytorch as pl

        pl.seed_everything(seed, workers=True)
    except Exception:
        pass


def _valid_smiles_mask(smiles: list[str]) -> np.ndarray:
    from rdkit import Chem

    mask = np.ones(len(smiles), dtype=bool)
    for i, s in enumerate(smiles):
        if Chem.MolFromSmiles(str(s)) is None:
            mask[i] = False
    return mask


def build_datapoints(df: pd.DataFrame, include_protein: bool):
    """Build Chemprop MoleculeDatapoint list; drop invalid SMILES."""
    from chemprop import data as cp_data

    smis = df["smiles"].astype(str).tolist()
    ys = df["pkd"].to_numpy(dtype=float).reshape(-1, 1)
    mask = _valid_smiles_mask(smis)
    smis_v = [s for s, m in zip(smis, mask) if m]
    ys_v = ys[mask]
    x_d = None
    if include_protein:
        feats = np.stack(
            [sequence_to_aac_dpc(seq) for seq, m in zip(df["sequence"].tolist(), mask) if m]
        ).astype(np.float32)
        x_d = feats

    dps = []
    for i, (smi, y) in enumerate(zip(smis_v, ys_v)):
        kwargs = {"y": y}
        if x_d is not None:
            kwargs["x_d"] = x_d[i]
        dp = cp_data.MoleculeDatapoint.from_smi(smi, **kwargs)
        if dp.mol is None:
            continue
        dps.append(dp)
    return dps, int(mask.sum()), int((~mask).sum())


def train_and_predict(
    train_df: pd.DataFrame,
    val_df: pd.DataFrame,
    test_df: pd.DataFrame,
    cfg: dict,
    include_protein: bool,
    ckpt_dir: Path,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Train Chemprop MPNN; return (y_true_test, y_pred_test, train_info)."""
    from chemprop import data as cp_data
    from chemprop import featurizers, models, nn
    from lightning import pytorch as pl
    from lightning.pytorch.callbacks import EarlyStopping, ModelCheckpoint

    cp = cfg["chemprop"]
    train_dps, n_tr, drop_tr = build_datapoints(train_df, include_protein)
    val_dps, n_va, drop_va = build_datapoints(val_df, include_protein)
    test_dps, n_te, drop_te = build_datapoints(test_df, include_protein)

    featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
    train_dset = cp_data.MoleculeDataset(train_dps, featurizer)
    val_dset = cp_data.MoleculeDataset(val_dps, featurizer)
    test_dset = cp_data.MoleculeDataset(test_dps, featurizer)

    scaler = train_dset.normalize_targets()
    val_dset.normalize_targets(scaler)

    x_d_transform = None
    n_xd = 0
    if include_protein:
        xd_scaler = train_dset.normalize_inputs("X_d")
        val_dset.normalize_inputs("X_d", xd_scaler)
        # Do NOT normalize test_dset: model X_d_transform applies at predict time
        x_d_transform = nn.ScaleTransform.from_standard_scaler(xd_scaler)
        n_xd = int(train_dps[0].x_d.shape[0]) if train_dps[0].x_d is not None else 420

    train_dset.cache = True
    val_dset.cache = True
    test_dset.cache = True

    bs = int(cp["batch_size"])
    nw = int(cp.get("num_workers", 0))
    train_loader = cp_data.build_dataloader(train_dset, batch_size=bs, num_workers=nw, shuffle=True)
    val_loader = cp_data.build_dataloader(val_dset, batch_size=bs, num_workers=nw, shuffle=False)
    test_loader = cp_data.build_dataloader(test_dset, batch_size=bs, num_workers=nw, shuffle=False)

    mp = nn.BondMessagePassing(
        d_h=int(cp.get("hidden_dim", 300)),
        depth=int(cp.get("depth", 3)),
        dropout=float(cp.get("dropout", 0.0)),
    )
    agg = nn.MeanAggregation()
    output_transform = nn.UnscaleTransform.from_standard_scaler(scaler)
    ffn_in = mp.output_dim + n_xd
    ffn = nn.RegressionFFN(
        input_dim=ffn_in,
        hidden_dim=int(cp.get("hidden_dim", 300)),
        n_layers=1,
        dropout=float(cp.get("dropout", 0.0)),
        output_transform=output_transform,
    )
    mpnn = models.MPNN(
        mp,
        agg,
        ffn,
        batch_norm=True,
        metrics=[nn.metrics.RMSE(), nn.metrics.MAE()],
        warmup_epochs=int(cp.get("warmup_epochs", 2)),
        init_lr=float(cp.get("init_lr", 1e-4)),
        max_lr=float(cp.get("max_lr", 1e-3)),
        final_lr=float(cp.get("final_lr", 1e-4)),
        X_d_transform=x_d_transform if include_protein else None,
    )

    ckpt_dir.mkdir(parents=True, exist_ok=True)
    checkpointing = ModelCheckpoint(
        dirpath=str(ckpt_dir),
        filename="best-{epoch}-{val_loss:.4f}",
        monitor="val_loss",
        mode="min",
        save_last=True,
        save_top_k=1,
    )
    early = EarlyStopping(
        monitor="val_loss",
        patience=int(cp.get("early_stopping_patience", 5)),
        mode="min",
    )
    trainer = pl.Trainer(
        logger=False,
        enable_checkpointing=True,
        enable_progress_bar=True,
        accelerator=cp.get("accelerator", "cpu"),
        devices=1,
        max_epochs=int(cp["max_epochs"]),
        callbacks=[checkpointing, early],
        deterministic=True,
    )
    trainer.fit(mpnn, train_loader, val_loader)

    # Predict with best checkpoint (torch>=2.6 defaults weights_only=True; Chemprop ckpts need False)
    best_path = checkpointing.best_model_path
    if best_path:
        try:
            model = models.MPNN.load_from_checkpoint(best_path, weights_only=False)
        except TypeError:
            # older lightning without weights_only kw
            model = models.MPNN.load_from_checkpoint(best_path)
        preds = trainer.predict(model, test_loader)
    else:
        model = mpnn
        preds = trainer.predict(model, test_loader)
    # preds: list of tensors (batch, n_tasks)
    y_pred = torch.cat(preds, dim=0).cpu().numpy().reshape(-1)
    y_true = np.array([float(dp.y[0]) for dp in test_dps], dtype=float)

    info = {
        "n_train": n_tr,
        "n_val": n_va,
        "n_test": n_te,
        "dropped_smiles_train": drop_tr,
        "dropped_smiles_val": drop_va,
        "dropped_smiles_test": drop_te,
        "n_xd": n_xd,
        "best_ckpt": str(best_path) if best_path else None,
        "best_val_loss": float(checkpointing.best_model_score)
        if checkpointing.best_model_score is not None
        else None,
        "epochs_run": int(trainer.current_epoch) + 1 if trainer.current_epoch is not None else None,
    }
    return y_true, y_pred, info, model, best_path


def _fmt(x) -> str:
    if x is None or (isinstance(x, float) and (np.isnan(x) or np.isinf(x))):
        return "nan"
    if isinstance(x, float):
        return f"{x:.4f}"
    return str(x)


def format_metrics_md(results: dict, cfg: dict) -> str:
    lines = [
        "# Virtual Screening v1.2 Chemprop (D-MPNN)",
        "",
        f"**Experiment freeze:** `{results['freeze_id']}`  ",
        f"**Data freeze:** `{results['data_freeze_id']}`  ",
        f"**Generated (America/Toronto):** {results['generated_at_america_toronto']}  ",
        f"**Binder threshold:** pKd >= {cfg['binder']['pkd_threshold']}  ",
        f"**Ligand encoder:** Chemprop BondMessagePassing (D-MPNN)  ",
        f"**Protein extras:** `{results.get('protein_mode', 'aac_dpc')}` as Chemprop x_d concat  ",
        f"**Seed:** {cfg['seed']}  ",
        f"**Max epochs:** {cfg['chemprop']['max_epochs']}  ",
        f"**Install:** chemprop=={results.get('chemprop_version', '?')}  ",
        "",
        "All numbers exported by `scripts/run_chemprop_v1.2.py` / `reproduce_v1.2.sh`. Do not hand-edit.",
        "Splits loaded from existing `data/splits/*.csv` (not regenerated).",
        "",
        "## Results (Chemprop)",
        "",
        "| Split | Model | Spearman rho | EF@1% | EF@5% | AUROC | n_test | n_binders |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for split_name in ("scaffold", "cold_protein"):
        block = results["splits"][split_name]
        for key, label in (
            ("ligand_only", "Chemprop ligand-only"),
            ("protein", f"Chemprop+{results.get('protein_mode', 'aac_dpc')}"),
        ):
            m = block[key]
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

    # Cite v1.1 HistGBM baselines
    lines += [
        "",
        "## v1.1 HistGBM baselines (cited; not retrained)",
        "",
        "| Split | Model | Spearman rho | EF@1% | EF@5% | AUROC |",
        "|---|---|---:|---:|---:|---:|",
    ]
    v11 = results.get("v1_1_baselines", {})
    for split_name in ("scaffold", "cold_protein"):
        b = v11.get(split_name, {})
        for key, label in (
            ("ligand_only", "HistGBM ligand-only (ECFP)"),
            ("aac_dpc", "HistGBM+AAC+DPC"),
        ):
            m = b.get(key, {})
            if not m:
                continue
            lines.append(
                "| {split} | {label} | {rho} | {ef1} | {ef5} | {auc} |".format(
                    split=split_name,
                    label=label,
                    rho=_fmt(m.get("spearman_rho")),
                    ef1=_fmt(m.get("ef_at_1pct")),
                    ef5=_fmt(m.get("ef_at_5pct")),
                    auc=_fmt(m.get("auroc")),
                )
            )

    # Success criteria
    cold = results["splits"]["cold_protein"]
    lo = cold["ligand_only"]
    pr = cold["protein"]
    d_rho = pr["spearman_rho"] - lo["spearman_rho"]
    d_ef1 = pr["ef_at_1pct"] - lo["ef_at_1pct"]
    d_ef5 = pr["ef_at_5pct"] - lo["ef_at_5pct"]
    primary = bool(d_rho >= 0.05 or d_ef1 > 0.5 or d_ef5 > 0.5)

    v11_cold_aac = v11.get("cold_protein", {}).get("aac_dpc", {})
    sec_rho = pr["spearman_rho"] - v11_cold_aac.get("spearman_rho", float("nan"))

    lines += [
        "",
        "## Success criteria (auto)",
        "",
        f"- Primary (cold_protein Chemprop+protein vs Chemprop ligand-only): "
        f"rho {_fmt(pr['spearman_rho'])} vs {_fmt(lo['spearman_rho'])} (Delta={_fmt(d_rho)}); "
        f"EF@1% {_fmt(pr['ef_at_1pct'])} vs {_fmt(lo['ef_at_1pct'])} (Delta={_fmt(d_ef1)}); "
        f"EF@5% {_fmt(pr['ef_at_5pct'])} vs {_fmt(lo['ef_at_5pct'])} (Delta={_fmt(d_ef5)}).",
        f"- Primary pass (DeltaSpearman>=0.05 OR clear EF lift): {'YES' if primary else 'NO'}.",
        f"- Secondary vs HistGBM+aac_dpc cold_protein: Chemprop+prot rho={_fmt(pr['spearman_rho'])} "
        f"vs HistGBM rho={_fmt(v11_cold_aac.get('spearman_rho'))} "
        f"(Delta={_fmt(sec_rho)}).",
        "",
        "## Notes",
        "",
        "- Same raw DAVIS + MANIFEST hashes; same split CSVs (loaded, not rewritten).",
        f"- Protein mode for x_d: `{results.get('protein_mode', 'aac_dpc')}` (best cheap v1.1 mode).",
        "- Chemprop/Lightning CPU runs are not guaranteed bit-stable; seed fixed; see reproduce_v1.2.sh.",
        f"- Reproduce float tolerance (strict): `{cfg['eval']['metric_float_tol']}`; "
        f"nondeterministic fallback: `{cfg['eval'].get('reproduce_tol_nondeterministic', 0.05)}`.",
        "",
        "## Artifacts",
        "",
    ]
    for k, v in sorted(results.get("models", {}).items()):
        lines.append(f"- `{k}`: `{v}`")
    lines.append("")
    return "\n".join(lines)


def load_v11_baselines(path: Path) -> dict:
    if not path.exists():
        return {}
    raw = json.loads(path.read_text())
    out = {}
    for split in ("scaffold", "cold_protein"):
        block = raw.get("splits", {}).get(split, {})
        out[split] = {
            "ligand_only": block.get("ligand_only", {}),
            "aac_dpc": block.get("protein", {}).get("aac_dpc", {}),
        }
    return out


def run(cfg_path: Path | None = None) -> dict:
    cfg_path = cfg_path or (ROOT / "configs" / "v1.2_chemprop.yaml")
    cfg = load_config(cfg_path)
    manifest = verify_data_freeze(cfg)
    set_seed(int(cfg["seed"]))

    import chemprop

    raw_dir = resolve(cfg, "raw_dir")
    splits_dir = resolve(cfg, "splits_dir")
    models_dir = resolve(cfg, "models_dir")
    metrics_md = resolve(cfg, "metrics_md")
    metrics_json = resolve(cfg, "metrics_json")
    models_dir.mkdir(parents=True, exist_ok=True)

    protein_mode = cfg.get("protein_features", {}).get("type", "aac_dpc")
    pairs = pd.read_csv(raw_dir / "davis_pairs.csv")
    scaf = load_existing_split(pairs, splits_dir / "scaffold_split.csv")
    cold = load_existing_split(pairs, splits_dir / "cold_protein_split.csv")

    v11_path = ROOT / cfg["paths"].get("v1_1_metrics_json", "artifacts/metrics_v1.1.json")
    v11 = load_v11_baselines(v11_path)

    results: dict = {
        "freeze_id": cfg["freeze_id"],
        "data_freeze_id": cfg.get("data_freeze_id", manifest.get("freeze_id")),
        "generated_at_america_toronto": datetime.now(TZ).isoformat(timespec="seconds"),
        "chemprop_version": chemprop.__version__,
        "protein_mode": protein_mode,
        "config": {
            "seed": cfg["seed"],
            "binder_pkd_threshold": cfg["binder"]["pkd_threshold"],
            "chemprop": cfg["chemprop"],
            "protein_features": cfg.get("protein_features"),
        },
        "splits": {},
        "models": {},
        "train_info": {},
        "v1_1_baselines": v11,
        "manifest_sha256": {
            "davis_tdc.tab": manifest["files"]["davis_tdc.tab"]["sha256"],
            "davis_pairs.csv": manifest["files"]["davis_pairs.csv"]["sha256"],
        },
        "nondeterminism_note": (
            "Chemprop+Lightning on CPU is not bit-stable across machines/runs; "
            "metrics exported from one fixed-seed run; reproduce_v1.2.sh documents tolerance."
        ),
    }

    thr = cfg["binder"]["pkd_threshold"]
    ef_fracs = cfg["eval"]["ef_fractions"]

    for split_name, dsplit in (("scaffold", scaf), ("cold_protein", cold)):
        train = dsplit[dsplit["split"] == "train"]
        val = dsplit[dsplit["split"] == "val"]
        test = dsplit[dsplit["split"] == "test"]
        results["splits"][split_name] = {}

        for include_protein, tag in ((False, "ligand_only"), (True, "protein")):
            label = "ligand_only" if not include_protein else f"protein_{protein_mode}"
            print(f"=== {split_name}: Chemprop {label} ===", flush=True)
            ckpt_dir = models_dir / f"ckpt_{split_name}_{tag}"
            y_true, y_pred, info, model, best_path = train_and_predict(
                train, val, test, cfg, include_protein, ckpt_dir
            )
            metrics = compute_metrics(y_true, y_pred, thr, ef_fracs)
            metrics["n_train"] = float(info["n_train"])
            metrics["n_test"] = float(info["n_test"])
            metrics["protein_type"] = protein_mode if include_protein else "none"
            metrics["protein_dim"] = float(info["n_xd"] if include_protein else 0)
            results["splits"][split_name][tag] = metrics
            results["train_info"][f"{split_name}_{tag}"] = info

            # Save model checkpoint path reference + predictions
            art_tag = f"{split_name}_{tag}"
            if best_path:
                dest = models_dir / f"chemprop_{art_tag}.ckpt"
                import shutil

                shutil.copy2(best_path, dest)
                results["models"][art_tag] = str(dest.relative_to(ROOT))
            pred_path = models_dir / f"preds_{art_tag}.npz"
            np.savez_compressed(pred_path, y_true=y_true, y_pred=y_pred)
            results["models"][f"{art_tag}_preds"] = str(pred_path.relative_to(ROOT))

            print(
                f"  rho={metrics['spearman_rho']:.4f} EF1={metrics['ef_at_1pct']:.4f} "
                f"EF5={metrics['ef_at_5pct']:.4f} AUROC={metrics['auroc']:.4f}",
                flush=True,
            )

    metrics_json.parent.mkdir(parents=True, exist_ok=True)
    with open(metrics_json, "w") as f:
        json.dump(results, f, indent=2, sort_keys=True)
        f.write("\n")
    with open(metrics_md, "w") as f:
        f.write(format_metrics_md(results, cfg))
    return results


def main() -> None:
    cfg_arg = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    results = run(cfg_arg)
    print("=== v1.2 Chemprop summary ===")
    for split, block in results["splits"].items():
        for tag, m in block.items():
            print(
                f"{split:14s} {tag:12s} rho={m['spearman_rho']:.4f}  "
                f"EF1={m['ef_at_1pct']:.4f}  EF5={m['ef_at_5pct']:.4f}"
            )
    print(f"metrics -> {ROOT / 'artifacts' / 'metrics_v1.2.md'}")


if __name__ == "__main__":
    main()
