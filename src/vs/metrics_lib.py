from __future__ import annotations

import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score


def enrichment_factor(y_true_binary: np.ndarray, y_score: np.ndarray, fraction: float) -> float:
    """EF@fraction: (hits in top fraction) / (expected hits by chance)."""
    y_true_binary = np.asarray(y_true_binary).astype(int)
    y_score = np.asarray(y_score)
    n = len(y_true_binary)
    if n == 0:
        return float("nan")
    n_actives = int(y_true_binary.sum())
    if n_actives == 0:
        return float("nan")
    k = max(1, int(np.ceil(fraction * n)))
    top_idx = np.argsort(-y_score)[:k]
    hits = int(y_true_binary[top_idx].sum())
    expected = n_actives * (k / n)
    if expected == 0:
        return float("nan")
    return float(hits / expected)


def compute_metrics(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    binder_threshold: float,
    ef_fractions: list[float],
) -> dict[str, float]:
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    out: dict[str, float] = {}
    if len(y_true) < 2 or np.allclose(y_true, y_true[0]):
        out["spearman_rho"] = float("nan")
    else:
        rho, _ = spearmanr(y_true, y_pred)
        out["spearman_rho"] = float(rho)
    y_bin = (y_true >= binder_threshold).astype(int)
    if y_bin.min() == y_bin.max():
        out["auroc"] = float("nan")
    else:
        out["auroc"] = float(roc_auc_score(y_bin, y_pred))
    for frac in ef_fractions:
        key = f"ef_at_{int(frac * 100)}pct"
        out[key] = enrichment_factor(y_bin, y_pred, frac)
    out["n"] = float(len(y_true))
    out["n_binders"] = float(y_bin.sum())
    out["binder_rate"] = float(y_bin.mean()) if len(y_bin) else float("nan")
    return out
