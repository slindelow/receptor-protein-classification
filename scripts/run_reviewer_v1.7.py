#!/usr/bin/env python
"""Sequence-identity holdouts, censored DAVIS labels, and an additive baseline.

Thresholds, seeds, and the censoring rule are read from configs/v1.7_reviewer.yaml.
Metrics are written to artifacts/metrics_v1.7.json. The prose in the README and
RESEARCH_REPORT.md is rendered from that file.

This script does not rewrite artifacts/metrics_v1.1.json, artifacts/metrics_uncertainty.json,
or artifacts/scores_histgbm_cold_protein.csv.
"""
from __future__ import annotations

import os

os.environ["OMP_NUM_THREADS"] = "4"
os.environ["MKL_NUM_THREADS"] = "4"
os.environ["OPENBLAS_NUM_THREADS"] = "4"

import argparse
import copy
import hashlib
import importlib.util
import json
import shutil
import subprocess
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from vs.data import load_davis_tab  # noqa: E402
from vs.features import (  # noqa: E402
    build_pair_matrix,
    sequence_to_aac_dpc,
    smiles_to_ecfp,
)
from vs.model import save_artifact, train_model  # noqa: E402
from vs.reviewer_stats import (  # noqa: E402
    balanced_two_way_shares,
    bin_index,
    bin_label,
    bootstrap_mean,
    dummy_r2,
    left_censored_latent_mean,
)
from vs.splits import _assign_by_group  # noqa: E402

CFG_PATH = ROOT / "configs" / "v1.7_reviewer.yaml"
PAIRS_CSV = ROOT / "data" / "raw" / "davis_pairs.csv"
TAB_PATH = ROOT / "data" / "raw" / "davis_tdc.tab"
MANIFEST = ROOT / "data" / "MANIFEST.json"
COLD_SPLIT = ROOT / "data" / "splits" / "cold_protein_split.csv"
SAVED_SCORES = ROOT / "artifacts" / "scores_histgbm_cold_protein.csv"
UNCERTAINTY = ROOT / "artifacts" / "metrics_uncertainty.json"
ART = ROOT / "artifacts" / "v1.7"
MODELS = ROOT / "artifacts" / "models_v1.7"
METRICS_JSON = ROOT / "artifacts" / "metrics_v1.7.json"
METRICS_MD = ROOT / "artifacts" / "metrics_v1.7.md"
DAVIS_URL = "https://dataverse.harvard.edu/api/access/datafile/5219748"
COVERAGE = 0.8


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rel(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT.resolve()))
    except ValueError:
        return path.name


def load_cfg() -> dict:
    with CFG_PATH.open() as handle:
        cfg = yaml.safe_load(handle)
    if int(cfg["omp_num_threads"]) != 4:
        raise SystemExit("configs/v1.7_reviewer.yaml omp_num_threads must stay 4")
    return cfg


def load_uncertainty_module():
    path = ROOT / "scripts" / "uncertainty_bootstrap.py"
    spec = importlib.util.spec_from_file_location("uncertainty_bootstrap_v17", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "kinase-binding-study"})
    with urllib.request.urlopen(request, timeout=120) as response, dest.open("wb") as handle:
        shutil.copyfileobj(response, handle)


def ensure_davis(manifest: dict) -> pd.DataFrame:
    expected_tab = manifest["files"]["davis_tdc.tab"]["sha256"]
    expected_csv = manifest["files"]["davis_pairs.csv"]["sha256"]
    if not TAB_PATH.exists() or sha256_file(TAB_PATH) != expected_tab:
        download(DAVIS_URL, TAB_PATH)
    got_tab = sha256_file(TAB_PATH)
    if got_tab != expected_tab:
        raise SystemExit(f"DAVIS tab sha256 {got_tab} does not match data/MANIFEST.json")
    if not (PAIRS_CSV.exists() and sha256_file(PAIRS_CSV) == expected_csv):
        frame = load_davis_tab(TAB_PATH)
        PAIRS_CSV.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(PAIRS_CSV, index=False)
    got_csv = sha256_file(PAIRS_CSV)
    if got_csv != expected_csv:
        raise SystemExit(f"DAVIS pairs sha256 {got_csv} does not match data/MANIFEST.json")
    frame = pd.read_csv(PAIRS_CSV)
    for column in ("pair_id", "drug_id", "target_id", "smiles", "sequence"):
        frame[column] = frame[column].astype(str)
    frame["kd_nm"] = frame["kd_nm"].astype(float)
    frame["pkd"] = frame["pkd"].astype(float)
    if frame["pair_id"].duplicated().any():
        raise SystemExit("duplicate pair_id in DAVIS pairs")
    nseq = frame.groupby("target_id")["sequence"].nunique()
    if int(nseq.min()) != 1 or int(nseq.max()) != 1:
        raise SystemExit("a target_id maps to more than one sequence")
    return frame


def ensure_mmseqs(cfg: dict) -> Path:
    binary = ROOT / "tools" / "mmseqs" / "bin" / "mmseqs"
    if binary.exists():
        return binary
    tarball = ROOT / "tools" / "mmseqs-linux-avx2.tar.gz"
    download(cfg["mmseqs"]["url"], tarball)
    digest = sha256_file(tarball)
    if digest != cfg["mmseqs"]["sha256"]:
        raise SystemExit(f"MMseqs2 tarball sha256 {digest} does not match the config pin")
    shutil.unpack_archive(tarball, ROOT / "tools")
    if not binary.exists():
        raise SystemExit("MMseqs2 binary was not found after unpacking the tarball")
    return binary


def mmseqs_version(binary: Path) -> str:
    out = subprocess.check_output([str(binary), "version"], text=True).strip()
    return out.splitlines()[-1].strip()


def write_fasta(pairs: pd.DataFrame, path: Path) -> None:
    seqs = pairs.drop_duplicates("target_id")[["target_id", "sequence"]]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for target_id, sequence in seqs.itertuples(index=False):
            if any(ch.isspace() for ch in target_id) or ">" in target_id:
                raise SystemExit(f"target id is not a safe FASTA header: {target_id}")
            handle.write(f">{target_id}\n{sequence}\n")


def run_cmd(argv: list[str]) -> None:
    printable = " ".join(rel(Path(part)) if str(part).startswith(str(ROOT)) else str(part) for part in argv)
    print(printable, flush=True)
    subprocess.check_call(argv)


def search_alignments(binary: Path, fasta: Path, cfg: dict, dest: Path) -> None:
    tmp = ROOT / "tools" / "tmp" / "search"
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True, exist_ok=True)
    mm = cfg["mmseqs"]
    argv = [
        str(binary),
        "easy-search",
        str(fasta),
        str(fasta),
        str(dest),
        str(tmp),
        "--min-seq-id",
        "0.0",
        "-c",
        "0.0",
        "--cov-mode",
        "0",
        "-s",
        str(mm["sensitivity"]),
        "--max-seqs",
        str(mm["max_seqs_search"]),
        "-e",
        str(mm["e_value"]),
        "--alignment-mode",
        str(mm["alignment_mode"]),
        "-a",
        str(mm["add_backtrace"]),
        "--seq-id-mode",
        "0",
        "--mask",
        str(mm["mask"]),
        "--threads",
        str(mm["threads"]),
        "--format-output",
        "query,target,fident,alnlen,qlen,tlen,qcov,tcov,nident",
    ]
    run_cmd(argv)
    shutil.rmtree(tmp, ignore_errors=True)


def cluster_kinases(binary: Path, fasta: Path, cfg: dict, threshold: float, dest_tsv: Path) -> None:
    tag = threshold_tag(threshold)
    prefix_dir = ART / f"clu_{tag}"
    if prefix_dir.exists():
        shutil.rmtree(prefix_dir)
    prefix_dir.mkdir(parents=True, exist_ok=True)
    tmp = ROOT / "tools" / "tmp" / f"clu_{tag}"
    if tmp.exists():
        shutil.rmtree(tmp)
    mm = cfg["mmseqs"]
    prefix = prefix_dir / "clu"
    argv = [
        str(binary),
        "easy-cluster",
        str(fasta),
        str(prefix),
        str(tmp),
        "--min-seq-id",
        f"{threshold:.2f}",
        "-c",
        str(mm["coverage"]),
        "--cov-mode",
        str(mm["cov_mode"]),
        "--cluster-mode",
        str(mm["cluster_mode"]),
        "--seq-id-mode",
        str(mm["seq_id_mode"]),
        "-s",
        str(mm["sensitivity"]),
        "--max-seqs",
        str(mm["max_seqs_cluster"]),
        "--single-step-clustering",
        str(mm["single_step_clustering"]),
        "-e",
        str(mm["e_value"]),
        "--mask",
        str(mm["mask"]),
        "--threads",
        str(mm["threads"]),
    ]
    run_cmd(argv)
    produced = Path(str(prefix) + "_cluster.tsv")
    if not produced.exists():
        raise SystemExit(f"MMseqs2 did not write {produced.name}")
    dest_tsv.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(produced, dest_tsv)
    shutil.rmtree(tmp, ignore_errors=True)


def threshold_tag(threshold: float) -> str:
    return f"c{int(round(threshold * 100)):02d}"


def split_csv_for(stem: str) -> Path:
    if stem == "cold_protein":
        return COLD_SPLIT
    if stem.startswith("seqid_c"):
        return ROOT / "data" / "splits" / f"seqid_cold_{stem.split('_', 1)[1]}.csv"
    raise SystemExit(f"no split csv for {stem}")


def as_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.strip().str.lower().isin(["true", "1", "1.0"])


def load_alignments(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["query", "target", "fident", "alnlen", "qlen", "tlen", "qcov", "tcov", "nident"],
        dtype={"query": str, "target": str},
    )
    return frame


def build_neighbors(frame: pd.DataFrame) -> dict[str, list[tuple]]:
    neighbors: dict[str, list[tuple]] = defaultdict(list)
    for row in frame.itertuples(index=False):
        if row.query == row.target:
            continue
        shorter = min(float(row.qlen), float(row.tlen))
        alnlen = float(row.alnlen)
        if shorter <= 0 or alnlen <= 0:
            continue
        nident = float(row.nident)
        id_short = nident / shorter
        fident = nident / alnlen
        both80 = bool(float(row.qcov) >= COVERAGE and float(row.tcov) >= COVERAGE)
        record = (id_short, fident, both80)
        neighbors[str(row.query)].append((str(row.target),) + record)
        neighbors[str(row.target)].append((str(row.query),) + record)
    return neighbors


def max_identity(kinase: str, pool: set[str], neighbors: dict) -> tuple[float, float]:
    best_short = 0.0
    best_fident = 0.0
    for other, id_short, fident, both80 in neighbors.get(kinase, ()):
        if other == kinase or other not in pool:
            continue
        if fident > best_fident:
            best_fident = fident
        if both80 and id_short > best_short:
            best_short = id_short
    return best_short, best_fident


def identity_map(kinases: list[str], pool: set[str], neighbors: dict) -> pd.DataFrame:
    rows = []
    for kinase in kinases:
        short, fident = max_identity(kinase, pool, neighbors)
        rows.append(
            {
                "target_id": kinase,
                "max_id_short_both80": short,
                "max_fident_any": fident,
            }
        )
    return pd.DataFrame(rows)


def write_sequence_split(pairs: pd.DataFrame, cluster_tsv: Path, threshold: float, cfg: dict) -> Path:
    tag = threshold_tag(threshold)
    clusters = pd.read_csv(cluster_tsv, sep="\t", header=None, names=["rep", "member"], dtype=str)
    mapping = dict(zip(clusters["member"], clusters["rep"]))
    missing = sorted(set(pairs["target_id"]) - set(mapping))
    if missing:
        raise SystemExit(f"{tag} cluster file is missing {len(missing)} kinases, including {missing[:3]}")
    work = pairs.copy().reset_index(drop=True)
    work["cluster_id"] = work["target_id"].map(mapping)
    assign = _assign_by_group(
        work["cluster_id"].astype(str).tolist(),
        float(cfg["splits"]["test_frac"]),
        float(cfg["splits"]["val_frac"]),
        int(cfg["seed"]),
    )
    work["split"] = assign
    work["split_protocol"] = f"seqid_cold_{tag}"
    dest = ROOT / "data" / "splits" / f"seqid_cold_{tag}.csv"
    columns = ["split_protocol", "pair_id", "drug_id", "target_id", "smiles", "split", "cluster_id"]
    dest.parent.mkdir(parents=True, exist_ok=True)
    work[columns].to_csv(dest, index=False)
    return dest


def censored_mask(kd: pd.Series, cap: float) -> np.ndarray:
    return np.isclose(kd.to_numpy(dtype=float), cap, rtol=0.0, atol=1e-6)


def feature_caches(pairs: pd.DataFrame, cfg: dict) -> tuple[dict, dict]:
    radius = int(cfg["fingerprint"]["radius"])
    n_bits = int(cfg["fingerprint"]["n_bits"])
    fingerprints = {}
    for smiles in pairs["smiles"].unique():
        fingerprints[smiles] = smiles_to_ecfp(smiles, radius=radius, n_bits=n_bits)
    proteins = {}
    for sequence in pairs["sequence"].unique():
        proteins[sequence] = sequence_to_aac_dpc(sequence)
    check_smiles = next(iter(fingerprints))
    check_seq = next(iter(proteins))
    built, mask = build_pair_matrix(
        [check_smiles],
        [check_seq],
        radius,
        n_bits,
        True,
        protein_type="aac_dpc",
    )
    manual = np.concatenate([fingerprints[check_smiles], proteins[check_seq]])
    if not mask[0] or not np.allclose(built[0], manual):
        raise SystemExit("feature construction does not match build_pair_matrix")
    return fingerprints, proteins


def design_matrix(frame: pd.DataFrame, fingerprints: dict, proteins: dict, include_protein: bool) -> np.ndarray:
    n_bits = len(next(iter(fingerprints.values())))
    width = n_bits + (420 if include_protein else 0)
    matrix = np.zeros((len(frame), width), dtype=np.float32)
    smiles = frame["smiles"].to_numpy()
    sequences = frame["sequence"].to_numpy()
    for i, (one_smiles, sequence) in enumerate(zip(smiles, sequences)):
        fingerprint = fingerprints[one_smiles]
        if fingerprint is None:
            raise SystemExit("invalid SMILES in a fit row")
        if include_protein:
            matrix[i] = np.concatenate([fingerprint, proteins[sequence]])
        else:
            matrix[i] = fingerprint
    return matrix


def protein_by_target(pairs: pd.DataFrame) -> dict[str, np.ndarray]:
    seqs = pairs.drop_duplicates("target_id")
    return {row.target_id: sequence_to_aac_dpc(row.sequence) for row in seqs.itertuples(index=False)}


def fit_hgb(X: np.ndarray, y: np.ndarray, cfg: dict, path: Path, meta: dict):
    if path.exists():
        from vs.model import load_artifact

        loaded = load_artifact(path)
        model = loaded["model"] if isinstance(loaded, dict) else loaded
        return model, int(getattr(model, "n_iter_", -1)), True
    model = train_model(X, y, cfg)
    meta = dict(meta)
    meta["n_iter"] = int(model.n_iter_)
    save_artifact(path, model, meta)
    return model, int(model.n_iter_), False


def ridge_alphas(cfg: dict) -> np.ndarray:
    add = cfg["additive"]
    return np.logspace(float(add["alpha_min_log10"]), float(add["alpha_max_log10"]), int(add["alpha_count"]))


def fit_additive(train: pd.DataFrame, proteins: dict[str, np.ndarray], alphas: np.ndarray) -> dict:
    mu = float(train["pkd"].mean())
    ligand_mean = train.groupby("drug_id")["pkd"].mean()
    alpha = ligand_mean - mu
    residual = train["pkd"].to_numpy(dtype=float) - mu - train["drug_id"].map(alpha).to_numpy(dtype=float)
    beta = pd.Series(residual, index=train.index).groupby(train["target_id"]).mean()
    ids = [str(i) for i in beta.index]
    matrix = np.vstack([proteins[i] for i in ids])
    target = beta.to_numpy(dtype=float)
    scaler = StandardScaler()
    scaled = scaler.fit_transform(matrix)
    ridge = RidgeCV(alphas=alphas)
    ridge.fit(scaled, target)
    chosen = float(np.ravel(ridge.alpha_)[0])
    train_r2 = float(ridge.score(scaled, target))
    return {
        "mu": mu,
        "alpha": alpha,
        "scaler": scaler,
        "ridge": ridge,
        "alpha_chosen": chosen,
        "train_kinase_r2": train_r2,
    }


def predict_additive(model: dict, frame: pd.DataFrame, proteins: dict[str, np.ndarray]) -> tuple[np.ndarray, int]:
    ids = list(dict.fromkeys(frame["target_id"].astype(str)))
    matrix = np.vstack([proteins[i] for i in ids])
    scaled = model["scaler"].transform(matrix)
    beta_hat = model["ridge"].predict(scaled)
    beta_map = dict(zip(ids, beta_hat))
    alpha = frame["drug_id"].map(model["alpha"])
    missing = int(alpha.isna().sum())
    prediction = (
        model["mu"]
        + alpha.fillna(0.0).to_numpy(dtype=float)
        + frame["target_id"].astype(str).map(beta_map).to_numpy(dtype=float)
    )
    return prediction, missing


def oracle_additive(train: pd.DataFrame, test: pd.DataFrame) -> np.ndarray:
    """Additive ceiling that uses each test kinase's own labels. Not a held-out score."""
    mu = float(train["pkd"].mean())
    alpha = train.groupby("drug_id")["pkd"].mean() - mu
    residual = test["pkd"].to_numpy(dtype=float) - mu - test["drug_id"].map(alpha).to_numpy(dtype=float)
    beta = pd.Series(residual, index=test.index).groupby(test["target_id"]).mean()
    return mu + test["drug_id"].map(alpha).to_numpy(dtype=float) + test["target_id"].map(beta).to_numpy(dtype=float)


def tobit_fit(X: np.ndarray, y: np.ndarray, censored: np.ndarray, cfg: dict, rounds: int, cap: float, path: Path, meta: dict):
    if path.exists():
        from vs.model import load_artifact

        loaded = load_artifact(path)
        model = loaded["model"] if isinstance(loaded, dict) else loaded
        history = loaded.get("meta", {}).get("tobit_history", []) if isinstance(loaded, dict) else []
        return model, history, True
    y_work = y.copy()
    history = []
    model = None
    for iteration in range(rounds):
        model = train_model(X, y_work, cfg)
        mu = np.asarray(model.predict(X), dtype=float)
        observed = ~censored
        if int(observed.sum()) < 2:
            raise SystemExit("Tobit fit has fewer than 2 uncensored training rows")
        sigma = float(np.std(y[observed] - mu[observed], ddof=1))
        sigma = max(sigma, 1e-3)
        y_work = y.copy()
        y_work[censored] = left_censored_latent_mean(mu[censored], sigma, cap)
        history.append(
            {
                "fit": iteration + 1,
                "sigma": sigma,
                "mean_imputed_censored": float(y_work[censored].mean()) if censored.any() else None,
                "n_iter": int(model.n_iter_),
            }
        )
    meta = dict(meta)
    meta["n_iter"] = int(model.n_iter_)
    meta["tobit_history"] = history
    meta["tobit_em_fits"] = rounds
    save_artifact(path, model, meta)
    return model, history, False


def joined_split(pairs: pd.DataFrame, split_csv: Path) -> pd.DataFrame:
    sid = pd.read_csv(split_csv)
    sid["pair_id"] = sid["pair_id"].astype(str)
    sid["split"] = sid["split"].astype(str)
    keep = ["pair_id", "split"]
    if "cluster_id" in sid.columns:
        sid["cluster_id"] = sid["cluster_id"].astype(str)
        keep.append("cluster_id")
    merged = sid[keep].merge(pairs, on="pair_id", how="inner")
    if len(merged) != len(sid):
        raise SystemExit(f"split join mismatch for {split_csv.name}: {len(merged)} vs {len(sid)}")
    return merged


def attach_identity(frame: pd.DataFrame, train_ids: set[str], neighbors: dict) -> pd.DataFrame:
    kinases = frame["target_id"].astype(str).unique().tolist()
    mapped = identity_map(kinases, train_ids, neighbors)
    out = frame.merge(mapped, on="target_id", how="left")
    if out["max_id_short_both80"].isna().any():
        raise SystemExit("identity merge failed")
    return out


def within_drug_range(frame: pd.DataFrame, column: str) -> float:
    spans = frame.groupby("drug_id")[column].agg(lambda s: float(s.max() - s.min()))
    return float(spans.max()) if len(spans) else float("nan")


def flatten_delta(block: dict) -> dict:
    return {
        "point": float(block["point"]),
        "low": float(block["ci95"]["low"]) if block["ci95"]["low"] is not None else None,
        "high": float(block["ci95"]["high"]) if block["ci95"]["high"] is not None else None,
        "includes_0": block["ci95"]["includes_0"],
        "p_value": float(block["test"]["p_value"]) if block.get("test") and block["test"].get("p_value") is not None else None,
    }


def flatten_contrast(raw: dict) -> dict:
    return {
        "n_rows": int(raw["n_rows"]),
        "n_kinases": int(raw["n_kinases"]),
        "spearman_a": float(raw["spearman_protein"]["point"]),
        "spearman_b": float(raw["spearman_ligand_only"]["point"]),
        "delta_spearman": flatten_delta(raw["delta_spearman_protein_minus_ligand"]),
        "ef_a": float(raw["ef_at_1pct_protein"]["point"]),
        "ef_b": float(raw["ef_at_1pct_ligand_only"]["point"]),
        "delta_ef": flatten_delta(raw["delta_ef_at_1pct_protein_minus_ligand"]),
    }


def contrast_or_point(ub, y, a, b, kinases, n_boot: int, seed: int, min_kinases: int) -> dict:
    y = np.asarray(y, dtype=float)
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    kinases = np.asarray(kinases)
    if len(y) < 2:
        return {
            "n_rows": int(len(y)),
            "n_kinases": int(pd.unique(kinases).size) if len(kinases) else 0,
            "spearman_a": None,
            "spearman_b": None,
            "delta_spearman": {"point": None, "low": None, "high": None, "includes_0": None, "p_value": None},
            "ef_a": None,
            "ef_b": None,
            "delta_ef": {"point": None, "low": None, "high": None, "includes_0": None, "p_value": None},
            "interval_computed": False,
        }
    n_kinases = int(pd.unique(kinases).size)
    if n_kinases >= min_kinases:
        raw = ub.bootstrap_kinase_cluster(y, a, b, kinases, n_boot, seed)
        out = flatten_contrast(raw)
        out["interval_computed"] = True
        return out
    rho_a = ub.spearman_rho(y, a)
    rho_b = ub.spearman_rho(y, b)
    from vs.metrics_lib import enrichment_factor

    y_bin = (y >= 7.0).astype(int)
    ef_a = float(enrichment_factor(y_bin, a, 0.01))
    ef_b = float(enrichment_factor(y_bin, b, 0.01))
    delta = None if not (np.isfinite(rho_a) and np.isfinite(rho_b)) else float(rho_a - rho_b)
    delta_ef = None if not (np.isfinite(ef_a) and np.isfinite(ef_b)) else float(ef_a - ef_b)
    return {
        "n_rows": int(len(y)),
        "n_kinases": n_kinases,
        "spearman_a": None if not np.isfinite(rho_a) else float(rho_a),
        "spearman_b": None if not np.isfinite(rho_b) else float(rho_b),
        "delta_spearman": {"point": delta, "low": None, "high": None, "includes_0": None, "p_value": None},
        "ef_a": None if not np.isfinite(ef_a) else ef_a,
        "ef_b": None if not np.isfinite(ef_b) else ef_b,
        "delta_ef": {"point": delta_ef, "low": None, "high": None, "includes_0": None, "p_value": None},
        "interval_computed": False,
    }


def within_delta(ub, frame: pd.DataFrame, col_a: str, col_b: str, n_boot: int, seed: int, min_kinases: int) -> dict:
    deltas = []
    for _, sub in frame.groupby("target_id", sort=False):
        rho_a = ub.spearman_rho(sub["pkd"].to_numpy(dtype=float), sub[col_a].to_numpy(dtype=float))
        rho_b = ub.spearman_rho(sub["pkd"].to_numpy(dtype=float), sub[col_b].to_numpy(dtype=float))
        if np.isfinite(rho_a) and np.isfinite(rho_b):
            deltas.append(float(rho_a - rho_b))
        else:
            deltas.append(np.nan)
    values = np.asarray(deltas, dtype=float)
    n_finite = int(np.isfinite(values).sum())
    if n_finite < min_kinases:
        finite = values[np.isfinite(values)]
        point = float(np.mean(finite)) if finite.size else None
        return {
            "point": point,
            "n_finite": n_finite,
            "n_kinases": int(frame["target_id"].nunique()),
            "ci95": {"low": None, "high": None, "includes_0": None},
            "p_value": None,
            "interval_computed": False,
        }
    out = bootstrap_mean(values, n_boot, seed)
    out["n_kinases"] = int(frame["target_id"].nunique())
    out["interval_computed"] = True
    return out


def bin_table(ub, frame: pd.DataFrame, column: str, edges: list[float], col_a: str, col_b: str, n_boot: int, seed: int, min_kinases: int) -> list[dict]:
    labels = frame.groupby("target_id", sort=False)[column].first()
    rows = []
    for index in range(len(edges) - 1):
        label = bin_label(edges, index)
        members = [kinase for kinase, value in labels.items() if bin_index(float(value), edges) == index]
        sub = frame[frame["target_id"].isin(members)]
        record = {
            "label": label,
            "n_kinases": int(len(members)),
            "n_pairs": int(len(sub)),
            "contrast": None,
        }
        if len(sub) >= 2 and len(members) >= 1:
            record["contrast"] = contrast_or_point(
                ub,
                sub["pkd"],
                sub[col_a],
                sub[col_b],
                sub["target_id"].astype(str),
                n_boot,
                seed,
                min_kinases,
            )
        rows.append(record)
    return rows


def identity_association(ub, frame: pd.DataFrame, id_col: str, col_a: str, col_b: str, n_boot: int, seed: int, min_kinases: int) -> dict:
    rows = []
    for kinase, sub in frame.groupby("target_id", sort=False):
        rho_a = ub.spearman_rho(sub["pkd"].to_numpy(dtype=float), sub[col_a].to_numpy(dtype=float))
        rho_b = ub.spearman_rho(sub["pkd"].to_numpy(dtype=float), sub[col_b].to_numpy(dtype=float))
        if not (np.isfinite(rho_a) and np.isfinite(rho_b)):
            continue
        rows.append((float(sub[id_col].iloc[0]), float(rho_a - rho_b)))
    if len(rows) < 3:
        return {"n_kinases": len(rows), "spearman": None, "low": None, "high": None, "includes_0": None, "p_value": None}
    ident = np.array([r[0] for r in rows], dtype=float)
    delta = np.array([r[1] for r in rows], dtype=float)
    point = ub.spearman_rho(ident, delta)
    if len(rows) < min_kinases or not np.isfinite(point):
        return {
            "n_kinases": len(rows),
            "spearman": None if not np.isfinite(point) else float(point),
            "low": None,
            "high": None,
            "includes_0": None,
            "p_value": None,
        }
    rng = np.random.default_rng(seed)
    samples = np.empty(n_boot, dtype=float)
    n = len(rows)
    for draw in range(n_boot):
        take = rng.integers(0, n, size=n)
        samples[draw] = ub.spearman_rho(ident[take], delta[take])
    finite = samples[np.isfinite(samples)]
    low, high = np.percentile(finite, [2.5, 97.5])
    recentered = finite - point
    count = int(np.sum(np.abs(recentered) >= abs(point)))
    return {
        "n_kinases": len(rows),
        "spearman": float(point),
        "low": float(low),
        "high": float(high),
        "includes_0": bool(low <= 0.0 <= high),
        "p_value": float((1 + count) / (len(finite) + 1)),
    }


def prepare_frame(
    pairs: pd.DataFrame,
    split_csv: Path,
    neighbors: dict,
    cap: float,
) -> tuple[pd.DataFrame, pd.DataFrame, set[str]]:
    merged = joined_split(pairs, split_csv)
    train_ids = set(merged.loc[merged["split"] == "train", "target_id"].astype(str))
    test = merged[merged["split"] == "test"].sort_values("pair_id").reset_index(drop=True)
    train = merged[merged["split"] == "train"].sort_values("pair_id").reset_index(drop=True)
    test = attach_identity(test, train_ids, neighbors)
    test["censored"] = censored_mask(test["kd_nm"], cap)
    train["censored"] = censored_mask(train["kd_nm"], cap)
    return train, test, train_ids


def fit_standard_pair(
    train: pd.DataFrame,
    test: pd.DataFrame,
    fingerprints: dict,
    proteins: dict,
    cfg: dict,
    stem: str,
    early_stopping: bool,
    fit_log: list,
) -> tuple[np.ndarray, np.ndarray, dict]:
    local = copy.deepcopy(cfg)
    local["model"]["early_stopping"] = bool(early_stopping)
    Xtr_l = design_matrix(train, fingerprints, proteins, False)
    Xtr_p = design_matrix(train, fingerprints, proteins, True)
    Xte_l = design_matrix(test, fingerprints, proteins, False)
    Xte_p = design_matrix(test, fingerprints, proteins, True)
    y = train["pkd"].to_numpy(dtype=float)
    info = {}
    preds = {}
    for name, Xtr, Xte, protein in (
        ("ligand_only", Xtr_l, Xte_l, "none"),
        ("aac_dpc", Xtr_p, Xte_p, "aac_dpc"),
    ):
        path = MODELS / f"hgb_{stem}_{name}.joblib"
        model, n_iter, loaded = fit_hgb(
            Xtr,
            y,
            local,
            path,
            {
                "split": stem,
                "protein": protein,
                "early_stopping": bool(early_stopping),
                "sorted_by": "pair_id",
                "freeze_id": cfg["freeze_id"],
            },
        )
        pred = np.asarray(model.predict(Xte), dtype=float)
        preds[name] = pred
        record = {
            "name": f"{stem}_{name}",
            "n_train": int(len(train)),
            "n_test": int(len(test)),
            "n_iter": n_iter,
            "early_stopping": bool(early_stopping),
            "loaded_existing_model": loaded,
            "model_file": rel(path),
        }
        fit_log.append(record)
        info[name] = record
        print(f"fit {stem} {name} trees {n_iter} loaded {loaded}", flush=True)
    return preds["ligand_only"], preds["aac_dpc"], info


def evaluate_split(ub, frame: pd.DataFrame, cfg: dict, col_a: str = "pred_aac_dpc", col_b: str = "pred_ligand_only") -> dict:
    n_boot = int(cfg["n_bootstrap"])
    seed = int(cfg["seed"])
    min_kinases = int(cfg["min_kinases_for_interval"])
    edges = [float(v) for v in cfg["identity_bins"]["edges"]]
    full = contrast_or_point(
        ub, frame["pkd"], frame[col_a], frame[col_b], frame["target_id"].astype(str), n_boot, seed, min_kinases
    )
    unc = frame.loc[~frame["censored"].to_numpy()]
    uncensored = contrast_or_point(
        ub, unc["pkd"], unc[col_a], unc[col_b], unc["target_id"].astype(str), n_boot, seed, min_kinases
    )
    return {
        "full_vs_ligand": full,
        "uncensored_vs_ligand": uncensored,
        "within_full_vs_ligand": within_delta(ub, frame, col_a, col_b, n_boot, seed, min_kinases),
        "within_uncensored_vs_ligand": within_delta(ub, unc, col_a, col_b, n_boot, seed, min_kinases),
        "bins_id_short": bin_table(ub, frame, "max_id_short_both80", edges, col_a, col_b, n_boot, seed, min_kinases),
        "bins_fident": bin_table(ub, frame, "max_fident_any", edges, col_a, col_b, n_boot, seed, min_kinases),
        "identity_association_id_short": identity_association(
            ub, frame, "max_id_short_both80", col_a, col_b, n_boot, seed, min_kinases
        ),
        "identity_association_fident": identity_association(
            ub, frame, "max_fident_any", col_a, col_b, n_boot, seed, min_kinases
        ),
    }


def split_census(train: pd.DataFrame, test: pd.DataFrame, threshold: float | None, neighbors: dict) -> dict:
    def cluster_count(frame: pd.DataFrame) -> int | None:
        if "cluster_id" not in frame.columns:
            return None
        return int(frame["cluster_id"].nunique())

    test_ids = test["target_id"].astype(str).unique().tolist()
    train_ids = set(train["target_id"].astype(str))
    ident = test.groupby("target_id")["max_id_short_both80"].first()
    ident_f = test.groupby("target_id")["max_fident_any"].first()
    n_leak = None
    n_leak_f = None
    if threshold is not None:
        n_leak = int((ident.to_numpy(dtype=float) >= float(threshold) - 1e-12).sum())
        n_leak_f = int((ident_f.to_numpy(dtype=float) >= float(threshold) - 1e-12).sum())
    sizes = None
    if "cluster_id" in pd.concat([train, test], ignore_index=True).columns:
        both = pd.concat([train, test], ignore_index=True)
        if "cluster_id" in both.columns:
            sizes = both.drop_duplicates("target_id").groupby("cluster_id").size()
    return {
        "threshold": threshold,
        "n_train_pairs": int(len(train)),
        "n_test_pairs": int(len(test)),
        "n_train_kinases": int(train["target_id"].nunique()),
        "n_test_kinases": int(test["target_id"].nunique()),
        "n_train_ligands": int(train["drug_id"].nunique()),
        "n_test_ligands": int(test["drug_id"].nunique()),
        "n_train_clusters": cluster_count(train),
        "n_test_clusters": cluster_count(test),
        "n_clusters": None if sizes is None else int(sizes.size),
        "largest_cluster": None if sizes is None else int(sizes.max()),
        "n_test_kinases_with_train_identity_at_least_threshold": n_leak,
        "n_test_kinases_with_train_fident_at_least_threshold": n_leak_f,
        "median_max_id_short_both80": float(ident.median()),
        "max_of_max_id_short_both80": float(ident.max()),
        "median_max_fident_any": float(ident_f.median()),
        "max_of_max_fident_any": float(ident_f.max()),
        "n_test_kinases_scored_for_identity": int(len(test_ids)),
        "n_train_kinases_in_identity_pool": int(len(train_ids)),
        "neighbor_index_kinases": int(len(neighbors)),
    }


def py(obj):
    if isinstance(obj, dict):
        return {str(k): py(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [py(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return py(obj.tolist())
    if isinstance(obj, (np.floating, float)):
        value = float(obj)
        if not np.isfinite(value):
            return None
        return value
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if obj is None or isinstance(obj, (str, int)):
        return obj
    return obj


def published_cold(ub) -> dict:
    payload = json.loads(UNCERTAINTY.read_text())
    match = None
    for comp in payload["comparisons"]:
        if comp.get("name", "").startswith("DAVIS cold-protein, HistGBM"):
            match = comp
            break
    if match is None or match.get("status") != "computed":
        raise SystemExit("published HistGBM cold-protein interval is not in artifacts/metrics_uncertainty.json")
    result = match["result"]
    scores = pd.read_csv(SAVED_SCORES)
    recomputed = ub.bootstrap_kinase_cluster(
        scores["pkd"].to_numpy(dtype=float),
        scores["pred_aac_dpc"].to_numpy(dtype=float),
        scores["pred_ligand_only"].to_numpy(dtype=float),
        scores["target_id"].astype(str).to_numpy(),
        int(payload["n_bootstrap"]),
        int(payload["seed"]),
    )
    published_point = float(result["delta_spearman_protein_minus_ligand"]["point"])
    recomputed_point = float(recomputed["delta_spearman_protein_minus_ligand"]["point"])
    published_low = float(result["delta_spearman_protein_minus_ligand"]["ci95"]["low"])
    recomputed_low = float(recomputed["delta_spearman_protein_minus_ligand"]["ci95"]["low"])
    return {
        "source": "artifacts/metrics_uncertainty.json",
        "recompute_abs_diff_point": abs(published_point - recomputed_point),
        "recompute_abs_diff_low": abs(published_low - recomputed_low),
        "recompute_matches": bool(
            abs(published_point - recomputed_point) <= 1e-9 and abs(published_low - recomputed_low) <= 1e-9
        ),
        "spearman_protein": float(result["spearman_protein"]["point"]),
        "spearman_ligand_only": float(result["spearman_ligand_only"]["point"]),
        "delta_spearman": flatten_delta(result["delta_spearman_protein_minus_ligand"]),
        "ef_protein": float(result["ef_at_1pct_protein"]["point"]),
        "ef_ligand_only": float(result["ef_at_1pct_ligand_only"]["point"]),
        "delta_ef": flatten_delta(result["delta_ef_at_1pct_protein_minus_ligand"]),
        "n_rows": int(result["n_rows"]),
        "n_kinases": int(result["n_kinases"]),
    }


def save_pred(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)


def build_predictions(cfg: dict, pairs: pd.DataFrame, neighbors: dict, force: bool) -> None:
    ART.mkdir(parents=True, exist_ok=True)
    MODELS.mkdir(parents=True, exist_ok=True)
    fingerprints, prot_by_seq = feature_caches(pairs, cfg)
    prot_by_id = protein_by_target(pairs)
    alphas = ridge_alphas(cfg)
    cap = float(cfg["censoring"]["pkd_cap"])
    kd_cap = float(cfg["censoring"]["kd_cap_nm"])
    fit_log = []
    if (ART / "fit_log.json").exists() and not force:
        fit_log = json.loads((ART / "fit_log.json").read_text())

    # Original cold-protein test uses the saved scores. New fits are only the secondary models.
    cold_train, cold_test, _ = prepare_frame(pairs, COLD_SPLIT, neighbors, kd_cap)
    cold_path = ART / "preds_cold_protein.csv"
    if force or not cold_path.exists():
        saved = pd.read_csv(SAVED_SCORES)
        saved["pair_id"] = saved["pair_id"].astype(str)
        saved["target_id"] = saved["target_id"].astype(str)
        base = cold_test.merge(
            saved[["pair_id", "pred_ligand_only", "pred_aac_dpc", "pkd"]].rename(columns={"pkd": "pkd_saved"}),
            on="pair_id",
            how="inner",
        )
        if len(base) != len(cold_test):
            raise SystemExit("saved cold-protein scores do not cover the frozen test split")
        if float(np.max(np.abs(base["pkd"] - base["pkd_saved"]))) > 1e-6:
            raise SystemExit("saved scores pkd does not match DAVIS pairs")
        additive_model = fit_additive(cold_train, prot_by_id, alphas)
        pred_add, missing = predict_additive(additive_model, base, prot_by_id)
        base["pred_additive"] = pred_add
        base["additive_alpha"] = additive_model["alpha_chosen"]
        base["additive_train_kinase_r2"] = additive_model["train_kinase_r2"]
        base["n_test_ligand_rows_missing_from_train"] = missing
        oracle = oracle_additive(cold_train, base)
        base["pred_oracle_additive"] = oracle
        save_pred(base, cold_path)
        (ART / "additive_cold_protein.json").write_text(
            json.dumps(
                {
                    "alpha_chosen": additive_model["alpha_chosen"],
                    "train_kinase_r2": additive_model["train_kinase_r2"],
                    "n_test_ligand_rows_missing_from_train": missing,
                },
                indent=2,
            )
            + "\n"
        )
        print("wrote saved-score cold-protein frame", flush=True)

    done_standard = set()
    for threshold in cfg["splits"]["thresholds"]:
        tag = threshold_tag(float(threshold))
        stem = f"seqid_{tag}"
        pred_path = ART / f"preds_{stem}.csv"
        if pred_path.exists() and not force:
            done_standard.add(stem)
            continue
        split_csv = ROOT / "data" / "splits" / f"seqid_cold_{tag}.csv"
        train, test, _ = prepare_frame(pairs, split_csv, neighbors, kd_cap)
        pred_l, pred_p, _ = fit_standard_pair(train, test, fingerprints, prot_by_seq, cfg, stem, True, fit_log)
        if within_drug_range(pd.DataFrame({"drug_id": test["drug_id"], "pred": pred_l}), "pred") > 1e-4:
            raise SystemExit(f"{stem} ligand-only predictions vary within a drug")
        additive_model = fit_additive(train, prot_by_id, alphas)
        pred_add, missing = predict_additive(additive_model, test, prot_by_id)
        test = test.copy()
        test["pred_ligand_only"] = pred_l
        test["pred_aac_dpc"] = pred_p
        test["pred_additive"] = pred_add
        test["pred_oracle_additive"] = oracle_additive(train, test)
        test["additive_alpha"] = additive_model["alpha_chosen"]
        test["additive_train_kinase_r2"] = additive_model["train_kinase_r2"]
        test["n_test_ligand_rows_missing_from_train"] = missing
        save_pred(test, pred_path)
        done_standard.add(stem)

    rounds = int(cfg["censoring"]["tobit_em_fits"])
    for stem in cfg["censoring"]["tobit_splits"]:
        pred_path = ART / f"preds_{stem}_tobit.csv"
        if pred_path.exists() and not force:
            continue
        train, test, _ = prepare_frame(pairs, split_csv_for(stem), neighbors, kd_cap)
        Xtr_l = design_matrix(train, fingerprints, prot_by_seq, False)
        Xtr_p = design_matrix(train, fingerprints, prot_by_seq, True)
        Xte_l = design_matrix(test, fingerprints, prot_by_seq, False)
        Xte_p = design_matrix(test, fingerprints, prot_by_seq, True)
        y = train["pkd"].to_numpy(dtype=float)
        cens = train["censored"].to_numpy()
        out = test[["pair_id", "target_id", "pkd", "censored"]].copy()
        for name, Xtr, Xte in (("ligand_only", Xtr_l, Xte_l), ("aac_dpc", Xtr_p, Xte_p)):
            path = MODELS / f"hgb_{stem}_tobit_{name}.joblib"
            model, history, loaded = tobit_fit(
                Xtr,
                y,
                cens,
                cfg,
                rounds,
                cap,
                path,
                {"split": stem, "protein": name, "task": "tobit", "freeze_id": cfg["freeze_id"]},
            )
            out[f"pred_tobit_{name}"] = np.asarray(model.predict(Xte), dtype=float)
            fit_log.append(
                {
                    "name": f"{stem}_tobit_{name}",
                    "n_train": int(len(train)),
                    "n_censored_train": int(cens.sum()),
                    "n_iter": int(model.n_iter_),
                    "loaded_existing_model": loaded,
                    "tobit_history": history,
                    "model_file": rel(path),
                }
            )
            print(f"tobit {stem} {name} loaded {loaded}", flush=True)
        save_pred(out, pred_path)

    for stem in cfg["censoring"]["drop_censored_splits"]:
        pred_path = ART / f"preds_{stem}_drop_censored.csv"
        if pred_path.exists() and not force:
            continue
        train, test, _ = prepare_frame(pairs, split_csv_for(stem), neighbors, kd_cap)
        kept = train.loc[~train["censored"]].copy()
        pred_l, pred_p, _ = fit_standard_pair(
            kept, test, fingerprints, prot_by_seq, cfg, f"{stem}_drop_censored", True, fit_log
        )
        out = test[["pair_id", "target_id", "pkd", "censored"]].copy()
        out["pred_drop_ligand_only"] = pred_l
        out["pred_drop_aac_dpc"] = pred_p
        save_pred(out, pred_path)

    for stem in cfg["censoring"]["early_stopping_off_splits"]:
        pred_path = ART / f"preds_{stem}_no_early_stop.csv"
        if pred_path.exists() and not force:
            continue
        train, test, _ = prepare_frame(pairs, split_csv_for(stem), neighbors, kd_cap)
        pred_l, pred_p, _ = fit_standard_pair(
            train, test, fingerprints, prot_by_seq, cfg, f"{stem}_no_early_stop", False, fit_log
        )
        out = test[["pair_id", "target_id", "pkd"]].copy()
        out["pred_noes_ligand_only"] = pred_l
        out["pred_noes_aac_dpc"] = pred_p
        save_pred(out, pred_path)

    (ART / "fit_log.json").write_text(json.dumps(py(fit_log), indent=2) + "\n")


def assemble(cfg: dict, pairs: pd.DataFrame, neighbors: dict, mmseqs_meta: dict) -> dict:
    ub = load_uncertainty_module()
    cap = float(cfg["censoring"]["kd_cap_nm"])
    cens = censored_mask(pairs["kd_nm"], cap)
    edges = [float(v) for v in cfg["identity_bins"]["edges"]]
    n_boot = int(cfg["n_bootstrap"])
    seed = int(cfg["seed"])
    min_kinases = int(cfg["min_kinases_for_interval"])
    shares = balanced_two_way_shares(pairs["pkd"], pairs["drug_id"], pairs["target_id"])
    unc_pairs = pairs.loc[~cens]
    r2_ligand = dummy_r2(unc_pairs["pkd"].to_numpy(), [unc_pairs["drug_id"].to_numpy()])
    r2_kinase = dummy_r2(unc_pairs["pkd"].to_numpy(), [unc_pairs["target_id"].to_numpy()])
    r2_add = dummy_r2(
        unc_pairs["pkd"].to_numpy(),
        [unc_pairs["drug_id"].to_numpy(), unc_pairs["target_id"].to_numpy()],
    )
    vectors = np.vstack([sequence_to_aac_dpc(seq) for seq in pairs.drop_duplicates("target_id")["sequence"]])
    # Round to avoid float noise creating fake uniqueness.
    rounded = np.round(vectors.astype(float), 8)
    unique_vectors = len({row.tobytes() for row in rounded})

    published = published_cold(ub)
    splits = {}

    def one_split(name: str, pred_path: Path, threshold: float | None, split_csv: Path) -> dict:
        frame = pd.read_csv(pred_path)
        frame["target_id"] = frame["target_id"].astype(str)
        frame["pair_id"] = frame["pair_id"].astype(str)
        frame["censored"] = as_bool(frame["censored"])
        train, test, train_ids = prepare_frame(pairs, split_csv, neighbors, cap)
        # Identity and census come from the joined split, predictions from the saved frame.
        # Reattach identity from this run so a stale column cannot drift.
        if set(frame["pair_id"].astype(str)) != set(test["pair_id"].astype(str)) or len(frame) != len(test):
            raise SystemExit(f"{name} prediction pair_ids do not match the test split")
        frame = frame.drop(columns=[c for c in ("max_id_short_both80", "max_fident_any") if c in frame.columns])
        ident = identity_map(frame["target_id"].unique().tolist(), train_ids, neighbors)
        frame = frame.merge(ident, on="target_id", how="left")
        block = split_census(train, test, threshold, neighbors)
        block["prediction_file"] = rel(pred_path)
        block["prediction_sha256"] = sha256_file(pred_path)
        block["source"] = (
            "artifacts/scores_histgbm_cold_protein.csv"
            if name == "cold_protein"
            else "refit on this split, rows sorted by pair_id"
        )
        evaluated = evaluate_split(ub, frame, cfg)
        block.update(evaluated)
        if "pred_additive" in frame.columns:
            block["full_vs_additive"] = contrast_or_point(
                ub,
                frame["pkd"],
                frame["pred_aac_dpc"],
                frame["pred_additive"],
                frame["target_id"].astype(str),
                n_boot,
                seed,
                min_kinases,
            )
            block["additive_vs_ligand"] = contrast_or_point(
                ub,
                frame["pkd"],
                frame["pred_additive"],
                frame["pred_ligand_only"],
                frame["target_id"].astype(str),
                n_boot,
                seed,
                min_kinases,
            )
            block["within_full_vs_additive"] = within_delta(
                ub, frame, "pred_aac_dpc", "pred_additive", n_boot, seed, min_kinases
            )
            block["additive_alpha"] = float(frame["additive_alpha"].iloc[0])
            block["additive_train_kinase_r2"] = float(frame["additive_train_kinase_r2"].iloc[0])
            block["n_test_ligand_rows_missing_from_train"] = int(frame["n_test_ligand_rows_missing_from_train"].iloc[0])
            block["oracle_spearman_uses_test_labels"] = float(
                ub.spearman_rho(frame["pkd"].to_numpy(dtype=float), frame["pred_oracle_additive"].to_numpy(dtype=float))
            )
        block["ligand_only_max_within_drug_range"] = within_drug_range(frame, "pred_ligand_only")
        tobit_path = ART / f"preds_{name}_tobit.csv"
        if tobit_path.exists():
            tobit = pd.read_csv(tobit_path)
            tobit["target_id"] = tobit["target_id"].astype(str)
            tobit["censored"] = as_bool(tobit["censored"])
            block["tobit_vs_tobit_ligand"] = contrast_or_point(
                ub,
                tobit["pkd"],
                tobit["pred_tobit_aac_dpc"],
                tobit["pred_tobit_ligand_only"],
                tobit["target_id"].astype(str),
                n_boot,
                seed,
                min_kinases,
            )
            unc = tobit.loc[~tobit["censored"].to_numpy()]
            block["tobit_uncensored_vs_tobit_ligand"] = contrast_or_point(
                ub,
                unc["pkd"],
                unc["pred_tobit_aac_dpc"],
                unc["pred_tobit_ligand_only"],
                unc["target_id"].astype(str),
                n_boot,
                seed,
                min_kinases,
            )
        drop_path = ART / f"preds_{name}_drop_censored.csv"
        if drop_path.exists():
            dropped = pd.read_csv(drop_path)
            dropped["target_id"] = dropped["target_id"].astype(str)
            dropped["censored"] = as_bool(dropped["censored"])
            unc = dropped.loc[~dropped["censored"].to_numpy()]
            block["drop_censored_uncensored"] = contrast_or_point(
                ub,
                unc["pkd"],
                unc["pred_drop_aac_dpc"],
                unc["pred_drop_ligand_only"],
                unc["target_id"].astype(str),
                n_boot,
                seed,
                min_kinases,
            )
            block["drop_censored_full_labels"] = contrast_or_point(
                ub,
                dropped["pkd"],
                dropped["pred_drop_aac_dpc"],
                dropped["pred_drop_ligand_only"],
                dropped["target_id"].astype(str),
                n_boot,
                seed,
                min_kinases,
            )
        noes_path = ART / f"preds_{name}_no_early_stop.csv"
        if noes_path.exists():
            noes = pd.read_csv(noes_path)
            noes["target_id"] = noes["target_id"].astype(str)
            # pair with labels from the main frame via pair_id if pkd present
            block["no_early_stop"] = contrast_or_point(
                ub,
                noes["pkd"],
                noes["pred_noes_aac_dpc"],
                noes["pred_noes_ligand_only"],
                noes["target_id"].astype(str),
                n_boot,
                seed,
                min_kinases,
            )
        return block

    splits["cold_protein"] = one_split("cold_protein", ART / "preds_cold_protein.csv", None, COLD_SPLIT)
    for threshold in cfg["splits"]["thresholds"]:
        tag = threshold_tag(float(threshold))
        name = f"seqid_{tag}"
        splits[name] = one_split(
            name,
            ART / f"preds_{name}.csv",
            float(threshold),
            ROOT / "data" / "splits" / f"seqid_cold_{tag}.csv",
        )

    fit_log = []
    if (ART / "fit_log.json").exists():
        fit_log = json.loads((ART / "fit_log.json").read_text())
    payload = {
        "freeze_id": cfg["freeze_id"],
        "data_freeze_id": cfg["data_freeze_id"],
        "seed": int(cfg["seed"]),
        "n_bootstrap": int(cfg["n_bootstrap"]),
        "omp_num_threads": 4,
        "config": "configs/v1.7_reviewer.yaml",
        "config_sha256": sha256_file(CFG_PATH),
        "identity_bin_edges": edges,
        "min_kinases_for_interval": int(cfg["min_kinases_for_interval"]),
        "censoring": {
            "kd_cap_nm": cap,
            "pkd_cap": float(cfg["censoring"]["pkd_cap"]),
            "n_pairs": int(len(pairs)),
            "n_censored": int(cens.sum()),
            "fraction_censored": float(cens.mean()),
            "n_kd_above_cap": int((pairs["kd_nm"] > cap + 1e-6).sum()),
            "n_uncensored": int((~cens).sum()),
            "definition": (
                "A pair is censored when Kd equals 10000 nM within an absolute tolerance of 1e-6. "
                "Kd is right-censored at that cap. pKd is left-censored at 5."
            ),
        },
        "variance_all_cells": shares,
        "variance_uncensored_r2": {
            "n_pairs": int(len(unc_pairs)),
            "r2_ligand": r2_ligand,
            "r2_kinase": r2_kinase,
            "r2_ligand_plus_kinase": r2_add,
            "share_outside_additive": float(1.0 - r2_add),
            "note": (
                "Uncensored rows are not a complete matrix. R^2 values are from least squares "
                "with an intercept and drop-first dummies. The share outside the additive fit "
                "is interaction plus residual error."
            ),
        },
        "composition_uniqueness": {
            "n_kinases": int(pairs["target_id"].nunique()),
            "n_unique_sequences": int(pairs.drop_duplicates("target_id")["sequence"].nunique()),
            "n_unique_aac_dpc_rounded_1e8": int(unique_vectors),
        },
        "published_cold_protein": published,
        "mmseqs": mmseqs_meta,
        "splits": splits,
        "fits": fit_log,
        "models_refit_for_original_cold_protein_headline": False,
    }
    return py(payload)


def prepare_inputs(cfg: dict, force: bool) -> tuple[pd.DataFrame, dict, dict]:
    manifest = json.loads(MANIFEST.read_text())
    if manifest.get("freeze_id") != cfg["data_freeze_id"]:
        raise SystemExit("MANIFEST freeze_id does not match the v1.7 config")
    pairs = ensure_davis(manifest)
    binary = ensure_mmseqs(cfg)
    version = mmseqs_version(binary)
    if cfg["mmseqs"]["version_commit"] not in version:
        raise SystemExit(f"MMseqs2 version {version} does not contain the pinned commit")
    ART.mkdir(parents=True, exist_ok=True)
    fasta = ART / "kinases.fasta"
    write_fasta(pairs, fasta)
    alignment = ART / "mmseqs_pairs.tsv"
    if force or not alignment.exists():
        search_alignments(binary, fasta, cfg, alignment)
    for threshold in cfg["splits"]["thresholds"]:
        dest = ART / f"clusters_{threshold_tag(float(threshold))}.tsv"
        split_csv = ROOT / "data" / "splits" / f"seqid_cold_{threshold_tag(float(threshold))}.csv"
        if force or not dest.exists():
            cluster_kinases(binary, fasta, cfg, float(threshold), dest)
        if force or not split_csv.exists():
            write_sequence_split(pairs, dest, float(threshold), cfg)
    align = load_alignments(alignment)
    self_hits = align.loc[align["query"].astype(str) == align["target"].astype(str)].copy()
    kinase_ids = set(pairs["target_id"].astype(str))
    self_ids = set(self_hits["query"].astype(str))
    missing_self = sorted(kinase_ids - self_ids)
    if missing_self:
        raise SystemExit(f"alignment is missing self-hits for {len(missing_self)} kinases, including {missing_self[:3]}")
    if self_hits.empty or float(self_hits["nident"].min()) <= 0:
        raise SystemExit("alignment self-hits have nident 0; backtrace was not applied")
    bad_self = self_hits.loc[
        (self_hits["nident"].to_numpy(dtype=float) != self_hits["qlen"].to_numpy(dtype=float))
        | (self_hits["nident"].to_numpy(dtype=float) != self_hits["tlen"].to_numpy(dtype=float))
    ]
    if len(bad_self):
        raise SystemExit(f"{len(bad_self)} self-hits have nident different from the sequence length")
    neighbors = build_neighbors(align)
    meta = {
        "release": cfg["mmseqs"]["release"],
        "version": version,
        "cluster_mode": int(cfg["mmseqs"]["cluster_mode"]),
        "seq_id_mode": int(cfg["mmseqs"]["seq_id_mode"]),
        "coverage": float(cfg["mmseqs"]["coverage"]),
        "cov_mode": int(cfg["mmseqs"]["cov_mode"]),
        "sensitivity": float(cfg["mmseqs"]["sensitivity"]),
        "e_value": str(cfg["mmseqs"]["e_value"]),
        "mask": int(cfg["mmseqs"]["mask"]),
        "single_step_clustering": int(cfg["mmseqs"]["single_step_clustering"]),
        "alignment_file": rel(alignment),
        "alignment_sha256": sha256_file(alignment),
        "n_self_hits": int(len(self_hits)),
        "min_self_nident": float(self_hits["nident"].min()),
        "self_hit_nident_equals_length": True,
        "identity_for_bins": (
            "Shorter-sequence identity is nident / min(query length, target length). "
            "A hit counts toward max identity to the training set only when query coverage "
            "and target coverage are both at least 0.8. The unrestricted maximum uses "
            "nident / alignment length with no coverage filter. Both directions of the search "
            "are considered, and the maximum is kept. The cluster threshold inside MMseqs2 "
            "uses --seq-id-mode 1, which is identity over the shorter sequence, and "
            "--cov-mode 0 with -c 0.8."
        ),
    }
    return pairs, neighbors, meta


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=["all", "metrics", "report"], default="all")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    cfg = load_cfg()
    if args.stage == "report":
        payload = json.loads(METRICS_JSON.read_text())
    else:
        pairs, neighbors, mmseqs_meta = prepare_inputs(cfg, args.force)
        if args.stage == "all":
            build_predictions(cfg, pairs, neighbors, args.force)
        payload = assemble(cfg, pairs, neighbors, mmseqs_meta)
        METRICS_JSON.write_text(json.dumps(payload, indent=2) + "\n")
        print(f"wrote {rel(METRICS_JSON)}", flush=True)
    from vs.v17_render import write_documents

    write_documents(payload, ROOT)
    cold = payload["splits"]["seqid_c40"]["full_vs_ligand"]["delta_spearman"]
    print(
        "seqid_c40 delta",
        cold["point"],
        "interval",
        cold["low"],
        cold["high"],
        "includes0",
        cold["includes_0"],
        flush=True,
    )


if __name__ == "__main__":
    main()
