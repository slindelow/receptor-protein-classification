#!/usr/bin/env python
"""v1.5 screening product: train HistGBM+ECFP+aac_dpc, screen a public library vs one DAVIS kinase."""
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
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from vs.config import load_config, resolve  # noqa: E402
from vs.features import build_pair_matrix, murcko_scaffold  # noqa: E402
from vs.model import save_artifact, train_model  # noqa: E402

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


def load_scaffold_train_val(pairs: pd.DataFrame, splits_dir: Path) -> pd.DataFrame:
    split_csv = splits_dir / "scaffold_split.csv"
    if not split_csv.exists():
        raise SystemExit(f"Missing split CSV (do not regenerate): {split_csv}")
    sid = pd.read_csv(split_csv)
    keep = [c for c in ["pair_id", "smiles", "sequence", "pkd", "drug_id", "target_id", "kd_nm"] if c in pairs.columns]
    merged = sid[["pair_id", "split"]].merge(pairs[keep], on="pair_id", how="inner")
    if len(merged) != len(sid):
        raise SystemExit(f"Split join size mismatch: sid={len(sid)} merged={len(merged)}")
    deploy = merged[merged["split"].isin(["train", "val"])].copy()
    return deploy


def _ok_smiles(s: str) -> bool:
    mol = Chem.MolFromSmiles(str(s))
    return mol is not None


def _canonical(s: str) -> str | None:
    mol = Chem.MolFromSmiles(str(s))
    if mol is None:
        return None
    try:
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        return None


def build_screen_library(cfg: dict, davis: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    screen = cfg["screen"]
    raw_dir = resolve(cfg, "raw_dir")
    url = screen["library_url"]
    hiv_path = raw_dir / "moleculenet_hiv.csv"
    if not hiv_path.exists():
        print(f"Downloading HIV library from {url} ...")
        urllib.request.urlretrieve(url, hiv_path)
    else:
        print(f"Using cached HIV file: {hiv_path}")

    hiv = pd.read_csv(hiv_path)
    if "smiles" not in hiv.columns:
        raise SystemExit(f"HIV file missing smiles column: {list(hiv.columns)}")

    # Canonicalize + validity filter
    rows = []
    for i, smi in enumerate(hiv["smiles"].astype(str)):
        can = _canonical(smi)
        if can is None:
            continue
        rows.append({"mol_id": f"hiv_{i}", "smiles": can, "source": "MoleculeNet_HIV"})
    hiv_clean = pd.DataFrame(rows).drop_duplicates("smiles")
    print(f"HIV valid unique canonical SMILES: {len(hiv_clean)}")

    davis_lig_set = set()
    davis_can_rows = []
    for _, r in davis.drop_duplicates("smiles").iterrows():
        can = _canonical(r["smiles"])
        if can is None:
            continue
        davis_lig_set.add(can)
        davis_can_rows.append(
            {
                "mol_id": f"davis_{r['drug_id']}",
                "smiles": can,
                "source": "DAVIS_spike_in",
            }
        )
    davis_lig = pd.DataFrame(davis_can_rows).drop_duplicates("smiles")

    # Prefer molecules NOT in DAVIS for the public pool; spike DAVIS back in for recovery.
    pool = hiv_clean[~hiv_clean["smiles"].isin(davis_lig_set)].copy()

    # Diverse subsample by Murcko scaffold (one mol per scaffold, then fill).
    n_target = int(screen["library_n"])
    rng = np.random.default_rng(cfg["seed"])
    scaffolds = []
    for smi in pool["smiles"]:
        scaffolds.append(murcko_scaffold(smi) or f"_noscaff_{hash(smi) & 0xFFFFFFFF:08x}")
    pool = pool.assign(scaffold=scaffolds)

    # Shuffle within each scaffold, take first of each scaffold, then pad.
    pool = pool.sample(frac=1.0, random_state=cfg["seed"]).reset_index(drop=True)
    first = pool.groupby("scaffold", sort=False).head(1)
    if len(first) >= n_target:
        chosen = first.iloc[:n_target].copy()
    else:
        rest = pool.drop(index=first.index)
        need = n_target - len(first)
        if need > 0 and len(rest) > 0:
            extra = rest.iloc[:need] if len(rest) <= need else rest.sample(n=need, random_state=cfg["seed"])
            chosen = pd.concat([first, extra], ignore_index=True)
        else:
            chosen = first.copy()
    chosen = chosen.drop(columns=["scaffold"], errors="ignore")

    if screen.get("spike_in_davis_ligands", True):
        lib = (
            pd.concat([chosen, davis_lig], ignore_index=True)
            .drop_duplicates("smiles")
            .sample(frac=1.0, random_state=cfg["seed"])
            .reset_index(drop=True)
        )
    else:
        lib = chosen.reset_index(drop=True)

    lib_path = resolve(cfg, "library_csv")
    lib_path.parent.mkdir(parents=True, exist_ok=True)
    lib[["mol_id", "smiles", "source"]].to_csv(lib_path, index=False)

    n_scaff = len({murcko_scaffold(s) or "" for s in lib["smiles"]})
    # empty scaffold string for failed: count unique non-empty + failed separately
    scaffs = [murcko_scaffold(s) for s in lib["smiles"]]
    n_unique_scaffolds = len(set(s for s in scaffs if s))

    manifest = {
        "freeze_id": cfg["freeze_id"],
        "generated_at_america_toronto": datetime.now(TZ).isoformat(timespec="seconds"),
        "library": {
            "n": int(len(lib)),
            "n_from_hiv": int((lib["source"] == "MoleculeNet_HIV").sum()),
            "n_davis_spike_in": int((lib["source"] == "DAVIS_spike_in").sum()),
            "n_unique_scaffolds": int(n_unique_scaffolds),
            "source_name": "MoleculeNet HIV (DeepChem S3 mirror)",
            "source_url": url,
            "source_note": "Public SMILES for screening only; HIV activity labels unused. DAVIS ligands spiked in for binder recovery analysis.",
            "hiv_raw_path": "data/raw/moleculenet_hiv.csv",
            "hiv_raw_sha256": sha256_file(hiv_path),
            "hiv_raw_bytes": hiv_path.stat().st_size,
            "library_path": str(lib_path.relative_to(ROOT)),
            "library_sha256": sha256_file(lib_path),
            "library_bytes": lib_path.stat().st_size,
            "seed": cfg["seed"],
            "requested_n": n_target,
        },
    }
    man_path = resolve(cfg, "library_manifest")
    man_path.parent.mkdir(parents=True, exist_ok=True)
    man_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"Screen library n={len(lib)} scaffolds={n_unique_scaffolds} -> {lib_path}")
    return lib, manifest


def train_screen_models(deploy: pd.DataFrame, cfg: dict) -> tuple[object, object, dict, dict]:
    radius = cfg["fingerprint"]["radius"]
    n_bits = cfg["fingerprint"]["n_bits"]
    ptype = cfg["protein_features"]["type"]

    print(f"Building features for n={len(deploy)} deploy pairs (ECFP + {ptype}) ...")
    X_full, m_full = build_pair_matrix(
        deploy["smiles"], deploy["sequence"], radius, n_bits, True, protein_type=ptype
    )
    y = deploy["pkd"].to_numpy()[m_full]
    X_full = X_full[m_full]
    print(f"Training full model on {len(y)} rows, dim={X_full.shape[1]} ...")
    model_full = train_model(X_full, y, cfg)

    X_lig, m_lig = build_pair_matrix(
        deploy["smiles"], None, radius, n_bits, False, protein_type=ptype
    )
    y_lig = deploy["pkd"].to_numpy()[m_lig]
    X_lig = X_lig[m_lig]
    print(f"Training ligand-only control on {len(y_lig)} rows ...")
    model_lig = train_model(X_lig, y_lig, cfg)

    meta_full = {
        "split": "scaffold_train_val",
        "mode": "full",
        "include_protein": True,
        "protein_type": ptype,
        "fingerprint": cfg["fingerprint"],
        "seed": cfg["seed"],
        "n_train": int(len(y)),
        "freeze_id": cfg["freeze_id"],
        "note": "Screening weights: HistGBM ECFP+aac_dpc fit on scaffold train+val (sofia-vs-v1.5-screen)",
    }
    meta_lig = {
        "split": "scaffold_train_val",
        "mode": "ligand_only",
        "include_protein": False,
        "protein_type": None,
        "fingerprint": cfg["fingerprint"],
        "seed": cfg["seed"],
        "n_train": int(len(y_lig)),
        "freeze_id": cfg["freeze_id"],
        "note": "Ligand-only control for screening (ECFP only)",
    }
    models_dir = resolve(cfg, "models_dir")
    save_artifact(models_dir / "hgb_screen_aac_dpc.joblib", model_full, meta_full)
    save_artifact(models_dir / "hgb_screen_ligand_only.joblib", model_lig, meta_lig)
    return model_full, model_lig, meta_full, meta_lig


def score_library(
    model,
    include_protein: bool,
    protein_type: str | None,
    lib: pd.DataFrame,
    target_seq: str,
    cfg: dict,
) -> np.ndarray:
    radius = cfg["fingerprint"]["radius"]
    n_bits = cfg["fingerprint"]["n_bits"]
    seqs = [target_seq] * len(lib) if include_protein else None
    ptype = protein_type or "aac"
    X, mask = build_pair_matrix(lib["smiles"], seqs, radius, n_bits, include_protein, protein_type=ptype)
    scores = np.full(len(lib), np.nan, dtype=float)
    if mask.any():
        scores[mask] = model.predict(X[mask])
    return scores, mask


def recovery_stats(
    ranked: pd.DataFrame,
    known_binders: set[str],
    known_all_for_target: pd.DataFrame,
) -> dict:
    """Rank of known binders (pKd>=7) among library; vs random baseline expectation."""
    n = len(ranked)
    binder_rows = ranked[ranked["smiles"].isin(known_binders)].copy()
    n_known_in_lib = int(len(binder_rows))
    ranks = binder_rows["rank"].astype(int).tolist() if n_known_in_lib else []
    # Mean reciprocal rank / enrichment in top-k
    out = {
        "n_library": n,
        "n_known_binders_total_for_target": int(len(known_binders)),
        "n_known_binders_in_library": n_known_in_lib,
        "known_binder_ranks": ranks,
        "known_binder_rank_mean": float(np.mean(ranks)) if ranks else None,
        "known_binder_rank_median": float(np.median(ranks)) if ranks else None,
        "n_in_top_50": int(sum(1 for r in ranks if r <= 50)),
        "n_in_top_100": int(sum(1 for r in ranks if r <= 100)),
        "n_in_top_500": int(sum(1 for r in ranks if r <= 500)),
    }
    # Random baseline: expected rank of a random item is (n+1)/2; expected count in top-k is k * n_known / n
    if n_known_in_lib and n:
        out["random_expected_rank_mean"] = float((n + 1) / 2.0)
        out["random_expected_in_top_50"] = float(50.0 * n_known_in_lib / n)
        out["random_expected_in_top_100"] = float(100.0 * n_known_in_lib / n)
        out["random_expected_in_top_500"] = float(500.0 * n_known_in_lib / n)
        # Enrichment vs random for top-500
        exp500 = out["random_expected_in_top_500"]
        out["enrichment_top_500_vs_random"] = (
            float(out["n_in_top_500"] / exp500) if exp500 > 0 else None
        )
    # Also list each known binder with true pKd and rank
    detail = []
    for _, row in known_all_for_target.iterrows():
        can = _canonical(row["smiles"])
        if can is None or can not in known_binders:
            continue
        hit = ranked[ranked["smiles"] == can]
        if len(hit) == 0:
            continue
        detail.append(
            {
                "drug_id": str(row["drug_id"]),
                "smiles": can,
                "true_pkd": float(row["pkd"]),
                "rank": int(hit.iloc[0]["rank"]),
                "pred_pkd": float(hit.iloc[0]["pred_pkd"]),
            }
        )
    detail = sorted(detail, key=lambda d: d["rank"])
    out["known_binder_detail"] = detail
    return out


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
    rec = results["recovery"]
    top20 = results["top20"]
    lines = [
        "# Virtual Screening v1.5 - Screening product metrics",
        "",
        f"**Experiment freeze:** `{results['freeze_id']}`  ",
        f"**Data freeze:** `{results['data_freeze_id']}`  ",
        f"**Generated (America/Toronto):** {results['generated_at_america_toronto']}  ",
        f"**Model:** HistGBM + ECFP(r={cfg['fingerprint']['radius']},{cfg['fingerprint']['n_bits']}) + aac_dpc  ",
        f"**Weights:** `{results['models']['full']}` (scaffold train+val, n_train={results['n_train']})  ",
        f"**Ligand-only control:** `{results['models']['ligand_only']}`  ",
        "",
        "Honest scope: these are **model scores** (predicted pKd), not wet-lab hits. No prospective claim.",
        "",
        "## Demo target",
        "",
        f"- **target_id:** `{sc['target_id']}` (DAVIS kinase; present in scaffold train, not a cold protein)",
        f"- **sequence length:** {sc['seq_len']}",
        f"- **DAVIS labeled pairs for target:** {sc['n_davis_pairs']}",
        f"- **DAVIS binders (pKd >= {cfg['binder']['pkd_threshold']}):** {sc['n_davis_binders']}",
        "",
        "## Screening library",
        "",
        f"- **n molecules:** {sc['n_library']}",
        f"- **n unique Murcko scaffolds:** {sc['n_unique_scaffolds']}",
        f"- **source:** {sc['library_source']} ({sc['library_url']})",
        f"- **raw SHA256:** `{sc['library_raw_sha256']}`",
        f"- **library CSV:** `{sc['library_csv']}`",
        f"- **HIV pool (non-DAVIS) + DAVIS ligand spike-in** for recovery analysis",
        "",
        "## Top-20 predicted pKd",
        "",
        "| rank | pred_pkd | ligand_only_pred_pkd | source | smiles |",
        "|---:|---:|---:|---|---|",
    ]
    for row in top20:
        lines.append(
            f"| {row['rank']} | {fmt(row['pred_pkd'])} | {fmt(row['ligand_only_pred_pkd'])} | {row['source']} | `{row['smiles']}` |"
        )
    lines += [
        "",
        "## Known-binder recovery (DAVIS labels for this target overlapping the library)",
        "",
        f"- Known binders in library: **{rec['n_known_binders_in_library']}** / {rec['n_known_binders_total_for_target']}",
        f"- Mean rank of known binders: **{fmt(rec.get('known_binder_rank_mean'))}** (random expected mean rank ≈ {fmt(rec.get('random_expected_rank_mean'))})",
        f"- Median rank: **{fmt(rec.get('known_binder_rank_median'))}**",
        f"- In top-50: {rec.get('n_in_top_50')} (random expected ≈ {fmt(rec.get('random_expected_in_top_50'))})",
        f"- In top-100: {rec.get('n_in_top_100')} (random expected ≈ {fmt(rec.get('random_expected_in_top_100'))})",
        f"- In top-500: {rec.get('n_in_top_500')} (random expected ≈ {fmt(rec.get('random_expected_in_top_500'))}); enrichment ≈ {fmt(rec.get('enrichment_top_500_vs_random'))}x",
        "",
        "Note: screening weights were fit on scaffold train+val, which includes some LCK-labeled pairs whose scaffolds were not held out. Recovery is a sanity check (known actives should rank high), not a cold-ligand prospective test.",
        "",
        "### Known binder detail (by rank)",
        "",
        "| rank | true_pkd | pred_pkd | drug_id | smiles |",
        "|---:|---:|---:|---|---|",
    ]
    for d in rec.get("known_binder_detail", []):
        lines.append(
            f"| {d['rank']} | {fmt(d['true_pkd'])} | {fmt(d['pred_pkd'])} | {d['drug_id']} | `{d['smiles']}` |"
        )
    lines += [
        "",
        "## Artifacts",
        "",
        f"- Ranked CSV: `{results['ranked_csv']}` (columns: smiles, pred_pkd, ligand_only_pred_pkd, rank, ...)",
        f"- Models dir: `{results['models_dir']}`",
        f"- Library manifest: `{sc['library_manifest']}`",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines))


def main(argv: list[str] | None = None) -> None:
    cfg_path = Path(argv[0]) if argv else ROOT / "configs" / "v1.5_screen.yaml"
    if len(sys.argv) > 1:
        cfg_path = Path(sys.argv[1])
    cfg = load_config(cfg_path)
    print(f"Config: {cfg_path}")
    verify_data_freeze(cfg)

    raw_dir = resolve(cfg, "raw_dir")
    davis = pd.read_csv(raw_dir / "davis_pairs.csv")
    deploy = load_scaffold_train_val(davis, resolve(cfg, "splits_dir"))
    print(f"Deploy set (scaffold train+val): n={len(deploy)}")

    lib, lib_manifest = build_screen_library(cfg, davis)

    model_full, model_lig, meta_full, meta_lig = train_screen_models(deploy, cfg)

    tid = cfg["screen"]["target_id"]
    tsub = davis[davis["target_id"] == tid]
    if tsub.empty:
        raise SystemExit(f"target_id {tid!r} not in DAVIS")
    seq = tsub["sequence"].iloc[0]
    binders = tsub[tsub["pkd"] >= cfg["binder"]["pkd_threshold"]].copy()
    binder_cans = set()
    for smi in binders["smiles"]:
        can = _canonical(smi)
        if can:
            binder_cans.add(can)

    print(f"Scoring library vs {tid} (seq_len={len(seq)}) ...")
    pred_full, mask_full = score_library(model_full, True, "aac_dpc", lib, seq, cfg)
    pred_lig, mask_lig = score_library(model_lig, False, None, lib, seq, cfg)

    out = lib.copy()
    out["pred_pkd"] = pred_full
    out["ligand_only_pred_pkd"] = pred_lig
    out["valid_smiles"] = mask_full & mask_lig
    out = out.sort_values("pred_pkd", ascending=False, na_position="last").reset_index(drop=True)
    out.insert(0, "rank", np.arange(1, len(out) + 1))
    # column order for export
    cols = ["smiles", "pred_pkd", "rank", "ligand_only_pred_pkd", "mol_id", "source", "valid_smiles"]
    ranked = out[cols]

    ranked_path = resolve(cfg, "ranked_csv")
    ranked_path.parent.mkdir(parents=True, exist_ok=True)
    ranked.to_csv(ranked_path, index=False)
    print(f"Wrote full ranking n={len(ranked)} -> {ranked_path}")

    # scaffolds on final library
    n_scaff = len({s for s in (murcko_scaffold(x) for x in ranked["smiles"]) if s})

    rec = recovery_stats(ranked, binder_cans, binders)
    top20 = ranked.head(20)[
        ["rank", "pred_pkd", "ligand_only_pred_pkd", "source", "smiles"]
    ].to_dict(orient="records")

    results = {
        "freeze_id": cfg["freeze_id"],
        "data_freeze_id": cfg.get("data_freeze_id"),
        "generated_at_america_toronto": datetime.now(TZ).isoformat(timespec="seconds"),
        "n_train": meta_full["n_train"],
        "models": {
            "full": "artifacts/models_v1.5/hgb_screen_aac_dpc.joblib",
            "ligand_only": "artifacts/models_v1.5/hgb_screen_ligand_only.joblib",
        },
        "models_dir": "artifacts/models_v1.5",
        "ranked_csv": str(ranked_path.relative_to(ROOT)),
        "screen": {
            "target_id": tid,
            "seq_len": len(seq),
            "n_davis_pairs": int(len(tsub)),
            "n_davis_binders": int(len(binders)),
            "n_library": int(len(ranked)),
            "n_unique_scaffolds": int(n_scaff),
            "library_source": "MoleculeNet HIV (DeepChem S3)",
            "library_url": cfg["screen"]["library_url"],
            "library_raw_sha256": lib_manifest["library"]["hiv_raw_sha256"],
            "library_csv": str(resolve(cfg, "library_csv").relative_to(ROOT)),
            "library_manifest": str(resolve(cfg, "library_manifest").relative_to(ROOT)),
        },
        "recovery": rec,
        "top20": top20,
        "top5": ranked.head(5)[["smiles", "pred_pkd", "rank"]].to_dict(orient="records"),
    }

    metrics_json = resolve(cfg, "metrics_json")
    metrics_json.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    write_metrics_md(results, cfg, resolve(cfg, "metrics_md"))
    print(f"Wrote {resolve(cfg, 'metrics_md')}")
    print("Top-5:")
    for r in results["top5"]:
        print(f"  rank={r['rank']} pred_pkd={r['pred_pkd']:.4f} {r['smiles'][:80]}")
    print(
        f"Recovery: {rec['n_known_binders_in_library']} binders; "
        f"mean_rank={rec.get('known_binder_rank_mean')}; "
        f"top500={rec.get('n_in_top_500')} (random~{rec.get('random_expected_in_top_500')})"
    )


if __name__ == "__main__":
    main()
