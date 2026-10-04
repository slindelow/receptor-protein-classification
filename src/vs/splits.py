from __future__ import annotations

from collections import defaultdict

import numpy as np
import pandas as pd

from .features import murcko_scaffold


def _assign_by_group(
    groups: list[str],
    test_frac: float,
    val_frac: float,
    seed: int,
) -> np.ndarray:
    """Assign each unique group to train/val/test by group size (greedy fill)."""
    rng = np.random.default_rng(seed)
    unique = sorted(set(groups))
    rng.shuffle(unique)
    # group -> indices
    g2idx = defaultdict(list)
    for i, g in enumerate(groups):
        g2idx[g].append(i)
    sizes = {g: len(g2idx[g]) for g in unique}
    n = len(groups)
    n_test = int(round(test_frac * n))
    n_val = int(round(val_frac * n))
    assign = np.full(n, "train", dtype=object)
    # fill test then val by walking shuffled groups
    for split_name, budget in (("test", n_test), ("val", n_val)):
        filled = 0
        remaining = [g for g in unique if all(assign[i] == "train" for i in g2idx[g])]
        # sort remaining large-first after shuffle order for stability
        for g in remaining:
            if filled >= budget:
                break
            for i in g2idx[g]:
                assign[i] = split_name
            filled += sizes[g]
    return assign


def scaffold_split(df: pd.DataFrame, test_frac: float, val_frac: float, seed: int) -> pd.DataFrame:
    scaffolds = [murcko_scaffold(s) or f"INVALID::{i}" for i, s in enumerate(df["smiles"])]
    assign = _assign_by_group(scaffolds, test_frac, val_frac, seed)
    out = df.copy()
    out["scaffold"] = scaffolds
    out["split"] = assign
    return out


def cold_protein_split(
    df: pd.DataFrame, test_frac: float, val_frac: float, seed: int
) -> pd.DataFrame:
    groups = df["target_id"].astype(str).tolist()
    assign = _assign_by_group(groups, test_frac, val_frac, seed)
    out = df.copy()
    out["split"] = assign
    return out


def write_split_ids(df: pd.DataFrame, path, split_name: str) -> None:
    cols = ["pair_id", "drug_id", "target_id", "smiles", "split"]
    if "scaffold" in df.columns:
        cols.append("scaffold")
    sub = df[cols].copy()
    sub.insert(0, "split_protocol", split_name)
    path.parent.mkdir(parents=True, exist_ok=True)
    sub.to_csv(path, index=False)
