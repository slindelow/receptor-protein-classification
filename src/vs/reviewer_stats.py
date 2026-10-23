"""Pure helpers for the v1.7 reviewer checks.

No file paths and no model fits. The training script imports these so the
ANOVA, censoring imputation, bins, and within-kinase interval can be tested
without DAVIS or MMseqs2.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import norm


def balanced_two_way_shares(y: np.ndarray, ligand: np.ndarray, kinase: np.ndarray) -> dict:
    """Variance shares for a complete ligand-by-kinase matrix with one observation per cell.

    The residual share is interaction plus measurement error. With one observation
    per cell those two pieces are not separable.
    """
    frame = pd.DataFrame(
        {
            "y": np.asarray(y, dtype=float),
            "ligand": np.asarray(ligand).astype(str),
            "kinase": np.asarray(kinase).astype(str),
        }
    )
    counts = frame.groupby(["ligand", "kinase"], sort=False).size()
    if int(counts.min()) != 1 or int(counts.max()) != 1:
        raise ValueError("balanced ANOVA requires one observation per ligand-kinase cell")
    pivot = frame.pivot(index="ligand", columns="kinase", values="y")
    if pivot.isna().any().any():
        raise ValueError("balanced ANOVA requires a complete matrix")
    matrix = pivot.to_numpy(dtype=float)
    n_ligand, n_kinase = matrix.shape
    grand = float(matrix.mean())
    ligand_effect = matrix.mean(axis=1) - grand
    kinase_effect = matrix.mean(axis=0) - grand
    residual = matrix - grand - ligand_effect[:, None] - kinase_effect[None, :]
    ss_total = float(np.sum((matrix - grand) ** 2))
    ss_ligand = float(n_kinase * np.sum(ligand_effect ** 2))
    ss_kinase = float(n_ligand * np.sum(kinase_effect ** 2))
    ss_residual = float(np.sum(residual ** 2))
    if ss_total <= 0:
        raise ValueError("total sum of squares is zero")
    return {
        "n_ligands": int(n_ligand),
        "n_kinases": int(n_kinase),
        "n_cells": int(n_ligand * n_kinase),
        "ss_total": ss_total,
        "ss_ligand": ss_ligand,
        "ss_kinase": ss_kinase,
        "ss_residual": ss_residual,
        "share_ligand": ss_ligand / ss_total,
        "share_kinase": ss_kinase / ss_total,
        "share_residual": ss_residual / ss_total,
        "share_sum": (ss_ligand + ss_kinase + ss_residual) / ss_total,
    }


def dummy_r2(y: np.ndarray, groups: list[np.ndarray]) -> float:
    """R^2 of an intercept plus drop-first dummies for each grouping factor."""
    y = np.asarray(y, dtype=float)
    pieces = [np.ones((len(y), 1), dtype=float)]
    frame = pd.DataFrame({f"g{i}": np.asarray(g).astype(str) for i, g in enumerate(groups)})
    for column in frame.columns:
        dummies = pd.get_dummies(frame[column], drop_first=True)
        if dummies.shape[1] == 0:
            continue
        pieces.append(dummies.to_numpy(dtype=float))
    design = np.column_stack(pieces)
    coef, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    fitted = design @ coef
    ss_total = float(np.sum((y - y.mean()) ** 2))
    if ss_total <= 0:
        return float("nan")
    ss_residual = float(np.sum((y - fitted) ** 2))
    return 1.0 - ss_residual / ss_total


def left_censored_latent_mean(mu: np.ndarray, sigma: float, censor_at: float) -> np.ndarray:
    """E[Y | Y <= censor_at] for Y ~ Normal(mu, sigma^2), elementwise in mu."""
    mu = np.asarray(mu, dtype=float)
    sigma = float(sigma)
    if not np.isfinite(sigma) or sigma <= 0:
        raise ValueError("sigma must be positive")
    z = (float(censor_at) - mu) / sigma
    mills = norm.pdf(z) / np.clip(norm.cdf(z), 1e-12, None)
    return mu - sigma * mills


def bin_index(value: float, edges: list[float]) -> int:
    """Bin index for edges. The last bin is closed on the right. Earlier bins are right-open."""
    edges = [float(e) for e in edges]
    value = float(value)
    if len(edges) < 2:
        raise ValueError("edges need at least two cut points")
    if value < edges[0] or value > edges[-1]:
        raise ValueError(f"value {value} is outside [{edges[0]}, {edges[-1]}]")
    last = len(edges) - 2
    for i in range(last):
        if edges[i] <= value < edges[i + 1]:
            return i
    return last


def bin_label(edges: list[float], index: int) -> str:
    left = edges[index]
    right = edges[index + 1]
    if index == len(edges) - 2:
        return f"[{left:.2f}, {right:.2f}]"
    return f"[{left:.2f}, {right:.2f})"


def bootstrap_mean(values: np.ndarray, n_boot: int, seed: int) -> dict:
    """Percentile interval for the mean of a kinase-level sample. NaN values are dropped."""
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    n = int(finite.size)
    if n == 0:
        return {
            "point": None,
            "n_finite": 0,
            "ci95": {"low": None, "high": None, "n_finite": 0, "includes_0": None},
            "p_value": None,
        }
    point = float(np.mean(finite))
    rng = np.random.default_rng(seed)
    samples = np.empty(n_boot, dtype=float)
    for draw in range(n_boot):
        take = finite[rng.integers(0, n, size=n)]
        samples[draw] = float(np.mean(take))
    low, high = np.percentile(samples, [2.5, 97.5])
    recentered = samples - point
    count = int(np.sum(np.abs(recentered) >= abs(point)))
    p_value = (1 + count) / (n_boot + 1)
    return {
        "point": point,
        "n_finite": n,
        "n_bootstrap": int(n_boot),
        "seed": int(seed),
        "ci95": {
            "low": float(low),
            "high": float(high),
            "n_finite": int(n_boot),
            "includes_0": bool(low <= 0.0 <= high),
        },
        "p_value": float(p_value),
        "n_as_extreme": count,
        "procedure": (
            "Resample the finite kinase-level values with replacement. "
            "The interval is the 2.5 and 97.5 percentiles. "
            "The two-sided p-value recenters each replicate by the observed mean "
            "and uses (1 + count of absolute recentered replicates at least as large "
            "as the absolute observed mean) / (n_bootstrap + 1)."
        ),
    }
