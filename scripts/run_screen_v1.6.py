#!/usr/bin/env python
"""v1.6 cold-target screen: hold out one DAVIS kinase, train HistGBM+ECFP+aac_dpc, rank the v1.5 library.

Model scores only. Not wet-lab hits.
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from vs.config import load_config, resolve  # noqa: E402
from vs.features import build_pair_matrix, sequence_to_aac_dpc  # noqa: E402
from vs.model import save_artifact, train_model  # noqa: E402
from rdkit import Chem  # noqa: E402

TZ = ZoneInfo("America/Toronto")


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
    raw_dir = resolve(cfg, "raw_dir")
    for key in ("davis_tdc.tab", "davis_pairs.csv"):
        meta = manifest["files"][key]
        path = ROOT / meta["path"]
        if not path.exists():
            path = raw_dir / Path(meta["path"]).name
        got = sha256_file(path)
        if got != meta["sha256"]:
            raise SystemExit(f"Hash mismatch for {key}: got {got} expected {meta['sha256']}")
    return manifest


def canonical(s: str) -> str | None:
    mol = Chem.MolFromSmiles(str(s))
    if mol is None:
        return None
    try:
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        return None


def binder_counts(davis: pd.DataFrame, threshold: float) -> pd.Series:
    return davis.groupby("target_id")["pkd"].apply(lambda s: int((s >= threshold).sum()))


def select_cold_target(davis: pd.DataFrame, cfg: dict) -> tuple[str, dict]:
    """Locked rule (applied before any model score).

    1. Drop LCK (v1.5 demo; it sat inside scaffold train).
    2. Count binders at pKd >= threshold.
    3. Drop target_ids that contain '(' (mutant or domain constructs).
    4. Keep targets with at least min_binders.
    5. Pick the binder count closest to LCK. Ties: higher binder count, then target_id.
    """
    screen = cfg["screen"]
    threshold = float(cfg["binder"]["pkd_threshold"])
    exclude = str(screen.get("exclude_target_id", "LCK"))
    min_binders = int(screen.get("min_binders", 10))
    counts = binder_counts(davis, threshold)
    if exclude not in counts.index:
        raise SystemExit(f"exclude target {exclude!r} not in DAVIS")
    ref = int(counts.loc[exclude])
    plain = counts.drop(labels=[exclude])
    plain = plain[~plain.index.to_series().str.contains(r"\(", regex=True)]
    eligible = plain[plain >= min_binders]
    if eligible.empty:
        raise SystemExit("No eligible cold target under the selection rule")
    # distance to ref, then -count, then name
    ranked = sorted(eligible.items(), key=lambda kv: (abs(int(kv[1]) - ref), -int(kv[1]), str(kv[0])))
    chosen, n_bind = ranked[0]
    info = {
        "rule": (
            "Exclude LCK. Drop target_ids containing '('. "
            f"Require n_binders >= {min_binders} at pKd >= {threshold}. "
            f"Pick binder count closest to LCK ({ref}). "
            "Tie break: more binders, then target_id."
        ),
        "excluded": exclude,
        "lck_n_binders": ref,
        "min_binders": min_binders,
        "n_eligible": int(len(eligible)),
        "chosen": str(chosen),
        "chosen_n_binders": int(n_bind),
        "runners_up": [
            {"target_id": str(t), "n_binders": int(c)} for t, c in ranked[1:6]
        ],
    }
    locked = screen.get("target_id")
    if locked and str(locked) != str(chosen):
        raise SystemExit(
            f"Config target_id={locked!r} does not match selection rule result {chosen!r}"
        )
    return str(chosen), info


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    na = float(np.linalg.norm(a))
    nb = float(np.linalg.norm(b))
    if na == 0 or nb == 0:
        return float("nan")
    return float(np.dot(a, b) / (na * nb))


def nearest_aac_dpc(davis: pd.DataFrame, target_id: str) -> dict:
    seqs = davis.groupby("target_id")["sequence"].first()
    query = sequence_to_aac_dpc(seqs.loc[target_id])
    best_id = None
    best = -1.0
    for tid, seq in seqs.items():
        if tid == target_id:
            continue
        c = cosine(query, sequence_to_aac_dpc(seq))
        if c > best:
            best = c
            best_id = str(tid)
    return {"nearest_target_id": best_id, "aac_dpc_cosine": best}


def load_existing_library(cfg: dict) -> tuple[pd.DataFrame, dict]:
    lib_path = resolve(cfg, "library_csv")
    man_path = resolve(cfg, "library_manifest")
    if not lib_path.exists():
        raise SystemExit(
            f"Screen library missing at {lib_path}. Refusing to redownload in v1.6; restore v1.5 library."
        )
    lib = pd.read_csv(lib_path)
    need = {"mol_id", "smiles", "source"}
    if not need.issubset(lib.columns):
        raise SystemExit(f"Library columns {list(lib.columns)} missing {need}")
    if lib["smiles"].duplicated().any():
        raise SystemExit("Library has duplicate SMILES; v1.5 file expected unique")
    manifest = json.loads(man_path.read_text()) if man_path.exists() else {}
    got = sha256_file(lib_path)
    expected = (manifest.get("library") or {}).get("library_sha256")
    if expected and got != expected:
        raise SystemExit(
            f"Library hash mismatch (will not rebuild): got {got} expected {expected}"
        )
    print(f"Reused library n={len(lib)} path={lib_path} sha256={got}")
    return lib, manifest


def train_cold_models(train_df: pd.DataFrame, cfg: dict, target_id: str):
    radius = cfg["fingerprint"]["radius"]
    n_bits = cfg["fingerprint"]["n_bits"]
    ptype = cfg["protein_features"]["type"]
    print(f"Featurizing n={len(train_df)} train pairs (ECFP + {ptype}) ...")
    X_full, m_full = build_pair_matrix(
        train_df["smiles"], train_df["sequence"], radius, n_bits, True, protein_type=ptype
    )
    y = train_df["pkd"].to_numpy()[m_full]
    X_full = X_full[m_full]
    print(f"Training protein model on {len(y)} rows, dim={X_full.shape[1]} ...")
    model_full = train_model(X_full, y, cfg)

    print("Featurizing ligand-only control ...")
    X_lig, m_lig = build_pair_matrix(
        train_df["smiles"], None, radius, n_bits, False, protein_type=ptype
    )
    y_lig = train_df["pkd"].to_numpy()[m_lig]
    X_lig = X_lig[m_lig]
    print(f"Training ligand-only on {len(y_lig)} rows, dim={X_lig.shape[1]} ...")
    model_lig = train_model(X_lig, y_lig, cfg)

    meta_full = {
        "split": "cold_target_holdout",
        "mode": "full",
        "include_protein": True,
        "protein_type": ptype,
        "fingerprint": cfg["fingerprint"],
        "seed": cfg["seed"],
        "n_train": int(len(y)),
        "held_out_target_id": target_id,
        "freeze_id": cfg["freeze_id"],
        "note": "HistGBM ECFP+aac_dpc fit on all DAVIS pairs except the held-out target_id",
    }
    meta_lig = {
        "split": "cold_target_holdout",
        "mode": "ligand_only",
        "include_protein": False,
        "protein_type": None,
        "fingerprint": cfg["fingerprint"],
        "seed": cfg["seed"],
        "n_train": int(len(y_lig)),
        "held_out_target_id": target_id,
        "freeze_id": cfg["freeze_id"],
        "note": "Ligand-only control. Same ligands may still appear on other kinases.",
    }
    models_dir = resolve(cfg, "models_dir")
    save_artifact(models_dir / "hgb_cold_screen_aac_dpc.joblib", model_full, meta_full)
    save_artifact(models_dir / "hgb_cold_screen_ligand_only.joblib", model_lig, meta_lig)
    return model_full, model_lig, meta_full, meta_lig


def score_library(model, include_protein, protein_type, lib, target_seq, cfg):
    radius = cfg["fingerprint"]["radius"]
    n_bits = cfg["fingerprint"]["n_bits"]
    seqs = [target_seq] * len(lib) if include_protein else None
    ptype = protein_type or "aac"
    X, mask = build_pair_matrix(lib["smiles"], seqs, radius, n_bits, include_protein, protein_type=ptype)
    scores = np.full(len(lib), np.nan, dtype=float)
    if mask.any():
        scores[mask] = model.predict(X[mask])
    return scores, mask


def assign_ranks(scores: np.ndarray, smiles: pd.Series) -> np.ndarray:
    """Rank 1 is best (highest score). NaN scores rank last. Ties broken by SMILES."""
    n = len(scores)
    order = sorted(
        range(n),
        key=lambda i: (
            1 if not np.isfinite(scores[i]) else 0,
            -scores[i] if np.isfinite(scores[i]) else 0.0,
            str(smiles.iloc[i]),
        ),
    )
    ranks = np.empty(n, dtype=int)
    for rank, i in enumerate(order, start=1):
        ranks[i] = rank
    return ranks


def recovery_block(ranks: np.ndarray, is_binder: np.ndarray, n: int) -> dict:
    k1 = max(1, int(np.ceil(0.01 * n)))
    binder_ranks = ranks[is_binder].astype(int)
    n_known = int(is_binder.sum())
    n_top1 = int(np.sum(binder_ranks <= k1)) if n_known else 0
    n_top500 = int(np.sum(binder_ranks <= 500)) if n_known else 0
    mean_rank = float(np.mean(binder_ranks)) if n_known else None
    median_rank = float(np.median(binder_ranks)) if n_known else None
    random_mean = float((n + 1) / 2.0)
    return {
        "n_library": n,
        "top1_k": k1,
        "top1_definition": "k = max(1, ceil(0.01 * n)); same cutoff as EF@1%",
        "n_known_binders_in_library": n_known,
        "n_in_top_1pct": n_top1,
        "fraction_top_1pct": float(n_top1 / n_known) if n_known else None,
        "random_expected_in_top_1pct": float(k1 * n_known / n) if n else None,
        "random_expected_fraction_top_1pct": float(k1 / n) if n else None,
        "n_in_top_500": n_top500,
        "fraction_top_500": float(n_top500 / n_known) if n_known else None,
        "random_expected_in_top_500": float(500.0 * n_known / n) if n else None,
        "random_expected_fraction_top_500": float(500.0 / n) if n else None,
        "mean_rank": mean_rank,
        "median_rank": median_rank,
        "random_expected_mean_rank": random_mean,
        "mean_rank_vs_random": (
            float(random_mean / mean_rank) if mean_rank and mean_rank > 0 else None
        ),
    }


def beats(protein: dict, ligand: dict) -> dict:
    """Protein beats ligand-only only if it is better or tied on all three recovery metrics and strict on at least one.

    Better: higher top-1% fraction, higher top-500 fraction, lower mean rank.
    """
    checks = {
        "top_1pct_fraction": protein["fraction_top_1pct"] - ligand["fraction_top_1pct"],
        "top_500_fraction": protein["fraction_top_500"] - ligand["fraction_top_500"],
        "mean_rank": ligand["mean_rank"] - protein["mean_rank"],
    }
    # positive delta means protein is better
    better = {k: v > 1e-12 for k, v in checks.items()}
    worse = {k: v < -1e-12 for k, v in checks.items()}
    if any(worse.values()) and any(better.values()):
        verdict = "mixed"
    elif any(better.values()) and not any(worse.values()):
        verdict = "yes"
    elif any(worse.values()) and not any(better.values()):
        verdict = "no"
    else:
        verdict = "tie"
    return {
        "protein_beats_ligand_only": verdict,
        "deltas_protein_minus_ligand_higher_is_better": {
            "top_1pct_fraction": float(checks["top_1pct_fraction"]),
            "top_500_fraction": float(checks["top_500_fraction"]),
            "mean_rank_improvement": float(checks["mean_rank"]),
        },
        "rule": (
            "yes if protein is strictly better on at least one of "
            "{top 1% fraction, top 500 fraction, mean rank} and worse on none; "
            "no if ligand-only is strictly better on at least one and protein on none; "
            "mixed if each wins at least one; tie if all equal. "
            "Lower mean rank is better."
        ),
    }


def fmt(x) -> str:
    if x is None:
        return "n/a"
    if isinstance(x, float):
        if np.isnan(x) or np.isinf(x):
            return "nan"
        return f"{x:.4f}"
    return str(x)


def write_metrics_md(results: dict, cfg: dict, path: Path) -> None:
    sc = results["screen"]
    prot = results["recovery_protein"]
    lig = results["recovery_ligand_only"]
    cmp_ = results["comparison"]
    lines = [
        "# Virtual Screening v1.6 - cold-target screen",
        "",
        f"**Experiment freeze:** `{results['freeze_id']}`  ",
        f"**Data freeze:** `{results['data_freeze_id']}`  ",
        f"**Generated (America/Toronto):** {results['generated_at_america_toronto']}  ",
        f"**Model:** HistGBM + ECFP(r={cfg['fingerprint']['radius']},{cfg['fingerprint']['n_bits']}) + aac_dpc  ",
        f"**Weights:** `{results['models']['full']}` (cold-target train, n_train={results['n_train']})  ",
        f"**Ligand-only control:** `{results['models']['ligand_only']}` (n_train={results['n_train_ligand_only']})  ",
        "",
        "Honest scope: these are **model scores** (predicted pKd), not wet-lab hits. No prospective claim.",
        "",
        "## Why this is not the v1.5 LCK screen",
        "",
        "v1.5 fit HistGBM on scaffold train+val, which still contained LCK pairs, then ranked a library for LCK. Recovery there is not a cold-target test.",
        "v1.6 holds out every pair for one other kinase and trains only on the remaining DAVIS pairs.",
        "",
        "## Held-out target",
        "",
        f"- **target_id:** `{sc['target_id']}`",
        f"- **selection rule:** {sc['selection']['rule']}",
        f"- **eligible kinases:** {sc['selection']['n_eligible']}; chosen binders={sc['selection']['chosen_n_binders']} (LCK had {sc['selection']['lck_n_binders']})",
        f"- **runners-up:** "
        + ", ".join(f"{r['target_id']} ({r['n_binders']})" for r in sc["selection"]["runners_up"]),
        f"- **sequence length:** {sc['seq_len']}",
        f"- **nearest other kinase by AAC+DPC cosine:** `{sc['nearest_target_id']}` ({fmt(sc['aac_dpc_cosine'])})",
        f"- **DAVIS labeled pairs held out:** {sc['n_held_out_pairs']}",
        f"- **DAVIS binders (pKd >= {cfg['binder']['pkd_threshold']}):** {sc['n_davis_binders']}",
        f"- **train pairs (all other target_ids):** {results['n_train']}",
        "",
        "Cold-target, not cold-ligand: the same 68 ligands still have labels on other kinases, so a ligand-only model can rank chemotypes that were potent elsewhere. The protein model never sees this target_id.",
        "",
        "## Screening library (reused, not rebuilt)",
        "",
        f"- **n molecules:** {sc['n_library']}",
        f"- **library CSV:** `{sc['library_csv']}`",
        f"- **library SHA256:** `{sc['library_sha256']}`",
        f"- **source note:** v1.5 MoleculeNet HIV subsample plus DAVIS ligand spike-in. HIV labels unused. Not redownloaded.",
        "",
        "## Recovery of known binders",
        "",
        f"Top 1% cutoff: rank <= {prot['top1_k']} ({prot['top1_definition']}). Top 500: rank <= 500.",
        "Fractions use known binders present in the library as the denominator.",
        "",
        "| model | n binders in library | fraction top 1% | n top 1% | fraction top 500 | n top 500 | mean rank | random mean rank |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
        (
            f"| protein (ECFP+aac_dpc) | {prot['n_known_binders_in_library']} | "
            f"{fmt(prot['fraction_top_1pct'])} | {prot['n_in_top_1pct']} | "
            f"{fmt(prot['fraction_top_500'])} | {prot['n_in_top_500']} | "
            f"{fmt(prot['mean_rank'])} | {fmt(prot['random_expected_mean_rank'])} |"
        ),
        (
            f"| ligand-only | {lig['n_known_binders_in_library']} | "
            f"{fmt(lig['fraction_top_1pct'])} | {lig['n_in_top_1pct']} | "
            f"{fmt(lig['fraction_top_500'])} | {lig['n_in_top_500']} | "
            f"{fmt(lig['mean_rank'])} | {fmt(lig['random_expected_mean_rank'])} |"
        ),
        "",
        (
            f"Random expected counts: top 1% ~ {fmt(prot['random_expected_in_top_1pct'])} "
            f"(fraction {fmt(prot['random_expected_fraction_top_1pct'])}); "
            f"top 500 ~ {fmt(prot['random_expected_in_top_500'])} "
            f"(fraction {fmt(prot['random_expected_fraction_top_500'])})."
        ),
        (
            f"Mean rank vs random (random/model, higher is better than chance): "
            f"protein {fmt(prot['mean_rank_vs_random'])}, ligand-only {fmt(lig['mean_rank_vs_random'])}."
        ),
        "",
        f"**Does the protein model beat ligand-only on recovery?** {cmp_['protein_beats_ligand_only']}",
        "",
        cmp_["rule"],
        "",
        (
            f"Deltas (positive means protein better): "
            f"top 1% fraction {fmt(cmp_['deltas_protein_minus_ligand_higher_is_better']['top_1pct_fraction'])}, "
            f"top 500 fraction {fmt(cmp_['deltas_protein_minus_ligand_higher_is_better']['top_500_fraction'])}, "
            f"mean-rank improvement {fmt(cmp_['deltas_protein_minus_ligand_higher_is_better']['mean_rank_improvement'])}."
        ),
        "",
        "## Known binder detail (sorted by protein rank)",
        "",
        "| protein_rank | ligand_only_rank | true_pkd | pred_pkd | ligand_only_pred_pkd | drug_id | smiles |",
        "|---:|---:|---:|---:|---:|---|---|",
    ]
    for d in results["known_binder_detail"]:
        lines.append(
            f"| {d['rank']} | {d['ligand_only_rank']} | {fmt(d['true_pkd'])} | "
            f"{fmt(d['pred_pkd'])} | {fmt(d['ligand_only_pred_pkd'])} | {d['drug_id']} | `{d['smiles']}` |"
        )
    lines += [
        "",
        "## Top-10 by protein score",
        "",
        "| rank | pred_pkd | ligand_only_rank | ligand_only_pred_pkd | source | known_binder | smiles |",
        "|---:|---:|---:|---:|---|---|---|",
    ]
    for row in results["top10"]:
        lines.append(
            f"| {row['rank']} | {fmt(row['pred_pkd'])} | {row['ligand_only_rank']} | "
            f"{fmt(row['ligand_only_pred_pkd'])} | {row['source']} | {row['is_known_binder']} | `{row['smiles']}` |"
        )
    lines += [
        "",
        "## Artifacts",
        "",
        f"- Ranked CSV: `{results['ranked_csv']}`",
        f"- Models dir: `{results['models_dir']}`",
        "",
        "Scores are predicted pKd from a model. They are not confirmed binders and not wet-lab hits.",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines))


def main() -> None:
    cfg_path = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "configs" / "v1.6_cold_screen.yaml"
    cfg = load_config(cfg_path)
    print(f"Config: {cfg_path}")
    verify_data_freeze(cfg)

    davis = pd.read_csv(resolve(cfg, "raw_dir") / "davis_pairs.csv")
    target_id, selection = select_cold_target(davis, cfg)
    print(f"Cold target: {target_id} ({selection})")

    held = davis[davis["target_id"] == target_id].copy()
    train_df = davis[davis["target_id"] != target_id].copy()
    if held.empty or (train_df["target_id"] == target_id).any():
        raise SystemExit("Holdout split failed")
    if len(train_df) + len(held) != len(davis):
        raise SystemExit("Train+holdout does not cover DAVIS")
    print(f"Holdout pairs={len(held)} train pairs={len(train_df)}")

    threshold = float(cfg["binder"]["pkd_threshold"])
    binders = held[held["pkd"] >= threshold].copy()
    binder_can = {}
    for _, row in binders.iterrows():
        can = canonical(row["smiles"])
        if can:
            binder_can[can] = row

    lib, lib_manifest = load_existing_library(cfg)
    near = nearest_aac_dpc(davis, target_id)
    seq = held["sequence"].iloc[0]

    model_full, model_lig, meta_full, meta_lig = train_cold_models(train_df, cfg, target_id)

    print(f"Scoring library vs {target_id} (seq_len={len(seq)}) ...")
    pred_full, mask_full = score_library(model_full, True, "aac_dpc", lib, seq, cfg)
    pred_lig, mask_lig = score_library(model_lig, False, None, lib, seq, cfg)

    smiles = lib["smiles"]
    rank_full = assign_ranks(pred_full, smiles)
    rank_lig = assign_ranks(pred_lig, smiles)
    is_binder = smiles.map(lambda s: s in binder_can).to_numpy()

    out = lib.copy()
    out["pred_pkd"] = pred_full
    out["ligand_only_pred_pkd"] = pred_lig
    out["rank"] = rank_full
    out["ligand_only_rank"] = rank_lig
    out["valid_smiles"] = mask_full & mask_lig
    out["is_known_binder"] = is_binder
    out = out.sort_values("rank", ascending=True).reset_index(drop=True)
    cols = [
        "smiles",
        "pred_pkd",
        "rank",
        "ligand_only_pred_pkd",
        "ligand_only_rank",
        "mol_id",
        "source",
        "valid_smiles",
        "is_known_binder",
    ]
    ranked = out[cols]
    ranked_path = resolve(cfg, "ranked_csv")
    ranked_path.parent.mkdir(parents=True, exist_ok=True)
    ranked.to_csv(ranked_path, index=False)
    print(f"Wrote ranking n={len(ranked)} -> {ranked_path}")

    n = len(ranked)
    rec_p = recovery_block(ranked["rank"].to_numpy(), ranked["is_known_binder"].to_numpy(), n)
    rec_l = recovery_block(
        ranked["ligand_only_rank"].to_numpy(), ranked["is_known_binder"].to_numpy(), n
    )
    comparison = beats(rec_p, rec_l)

    detail = []
    for _, row in ranked[ranked["is_known_binder"]].iterrows():
        src = binder_can[row["smiles"]]
        detail.append(
            {
                "drug_id": str(src["drug_id"]),
                "smiles": row["smiles"],
                "true_pkd": float(src["pkd"]),
                "rank": int(row["rank"]),
                "ligand_only_rank": int(row["ligand_only_rank"]),
                "pred_pkd": float(row["pred_pkd"]),
                "ligand_only_pred_pkd": float(row["ligand_only_pred_pkd"]),
            }
        )
    detail = sorted(detail, key=lambda d: d["rank"])

    lib_sha = sha256_file(resolve(cfg, "library_csv"))
    top10 = ranked.head(10)[
        ["rank", "pred_pkd", "ligand_only_rank", "ligand_only_pred_pkd", "source", "is_known_binder", "smiles"]
    ].to_dict(orient="records")
    for row in top10:
        row["rank"] = int(row["rank"])
        row["ligand_only_rank"] = int(row["ligand_only_rank"])
        row["pred_pkd"] = float(row["pred_pkd"])
        row["ligand_only_pred_pkd"] = float(row["ligand_only_pred_pkd"])
        row["is_known_binder"] = bool(row["is_known_binder"])

    results = {
        "freeze_id": cfg["freeze_id"],
        "data_freeze_id": cfg.get("data_freeze_id"),
        "generated_at_america_toronto": datetime.now(TZ).isoformat(timespec="seconds"),
        "n_train": meta_full["n_train"],
        "n_train_ligand_only": meta_lig["n_train"],
        "n_held_out": int(len(held)),
        "models": {
            "full": "artifacts/models_v1.6/hgb_cold_screen_aac_dpc.joblib",
            "ligand_only": "artifacts/models_v1.6/hgb_cold_screen_ligand_only.joblib",
        },
        "models_dir": "artifacts/models_v1.6",
        "ranked_csv": str(ranked_path.relative_to(ROOT)),
        "screen": {
            "target_id": target_id,
            "seq_len": int(len(seq)),
            "n_held_out_pairs": int(len(held)),
            "n_davis_binders": int(len(binders)),
            "n_known_binders_in_library": int(is_binder.sum()),
            "n_library": n,
            "library_csv": str(resolve(cfg, "library_csv").relative_to(ROOT)),
            "library_sha256": lib_sha,
            "library_manifest_freeze_id": lib_manifest.get("freeze_id"),
            "selection": selection,
            "nearest_target_id": near["nearest_target_id"],
            "aac_dpc_cosine": near["aac_dpc_cosine"],
        },
        "recovery_protein": rec_p,
        "recovery_ligand_only": rec_l,
        "comparison": comparison,
        "known_binder_detail": detail,
        "top10": top10,
    }
    metrics_json = resolve(cfg, "metrics_json")
    metrics_json.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    write_metrics_md(results, cfg, resolve(cfg, "metrics_md"))
    print(f"Wrote {resolve(cfg, 'metrics_md')}")
    print(
        f"RECOVERY protein top1%={rec_p['fraction_top_1pct']:.4f} "
        f"top500={rec_p['fraction_top_500']:.4f} mean_rank={rec_p['mean_rank']:.2f} | "
        f"ligand top1%={rec_l['fraction_top_1pct']:.4f} "
        f"top500={rec_l['fraction_top_500']:.4f} mean_rank={rec_l['mean_rank']:.2f} | "
        f"beats={comparison['protein_beats_ligand_only']}"
    )


if __name__ == "__main__":
    main()
