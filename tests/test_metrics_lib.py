"""Hand-checked Spearman and enrichment on a tiny input.

Imports only src/vs/metrics_lib.py. No model fit and no chemistry stack.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from vs.metrics_lib import compute_metrics, enrichment_factor


def test_spearman_on_a_four_point_swap():
    # Ranks of [1, 2, 3, 4] versus [1, 2, 4, 3].
    # Squared rank gaps are 0, 0, 1, 1. rho = 1 - 6*2 / (4*(16-1)) = 0.8.
    y_true = np.array([1.0, 2.0, 3.0, 4.0])
    y_pred = np.array([1.0, 2.0, 4.0, 3.0])
    metrics = compute_metrics(y_true, y_pred, binder_threshold=3.0, ef_fractions=[0.5])
    assert math.isclose(metrics["spearman_rho"], 0.8, abs_tol=1e-12)

    reversed_pred = np.array([4.0, 3.0, 2.0, 1.0])
    reversed_metrics = compute_metrics(
        y_true, reversed_pred, binder_threshold=3.0, ef_fractions=[0.5]
    )
    assert math.isclose(reversed_metrics["spearman_rho"], -1.0, abs_tol=1e-12)


def test_enrichment_at_k_uses_ceil_of_the_fraction():
    # n = 10, one active, sitting at rank 3.
    # fraction 0.25 -> k = ceil(2.5) = 3, so the active is inside the slice.
    # expected hits = 1 * 3/10 = 0.3, so EF = 1/0.3.
    # A floor cutoff (k = 2) would have missed it.
    y_true = np.array([0, 0, 1, 0, 0, 0, 0, 0, 0, 0])
    y_score = np.array([0.9, 0.8, 0.7, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
    ef = enrichment_factor(y_true, y_score, 0.25)
    assert math.isclose(ef, 1.0 / 0.3, abs_tol=1e-12)

    # Same labels, active at the bottom of the list. Top 3 are inactive, so EF is 0.
    missed = np.array([0.9, 0.8, 0.0, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1])
    assert enrichment_factor(y_true, missed, 0.25) == 0.0
