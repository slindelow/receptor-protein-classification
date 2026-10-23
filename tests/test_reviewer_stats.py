"""Checks for the v1.7 ANOVA, censoring imputation, and bins."""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from scipy.stats import norm

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from vs.reviewer_stats import (
    balanced_two_way_shares,
    bin_index,
    bootstrap_mean,
    left_censored_latent_mean,
)
from vs.splits import _assign_by_group


def test_balanced_anova_recovers_main_effects_and_a_zero_residual():
    # Rows are ligands. Columns are kinases.
    # Ligand effects relative to the grand mean are -1 and +1.
    # Kinase effects are -2 and +2. No interaction.
    matrix = np.array([[0.0, 4.0], [2.0, 6.0]])
    ligand = np.array(["L0", "L0", "L1", "L1"])
    kinase = np.array(["K0", "K1", "K0", "K1"])
    shares = balanced_two_way_shares(matrix.ravel(), ligand, kinase)
    assert math.isclose(shares["share_ligand"], 0.2, abs_tol=1e-12)
    assert math.isclose(shares["share_kinase"], 0.8, abs_tol=1e-12)
    assert math.isclose(shares["share_residual"], 0.0, abs_tol=1e-12)
    assert math.isclose(shares["share_sum"], 1.0, abs_tol=1e-12)


def test_left_censored_mean_matches_the_standard_normal_mills_ratio():
    expected = -norm.pdf(0.0) / norm.cdf(0.0)
    got = float(left_censored_latent_mean(np.array([0.0]), 1.0, 0.0)[0])
    assert math.isclose(got, expected, abs_tol=1e-10)
    # A latent mean far below the cap stays far below the cap.
    below = float(left_censored_latent_mean(np.array([-5.0]), 1.0, 5.0)[0])
    assert below < -4.0


def test_bin_edges_are_right_open_until_the_last_bin():
    edges = [0.0, 0.30, 0.40, 0.50, 0.70, 1.0]
    assert bin_index(0.0, edges) == 0
    assert bin_index(0.299, edges) == 0
    assert bin_index(0.30, edges) == 1
    assert bin_index(0.40, edges) == 2
    assert bin_index(0.70, edges) == 4
    assert bin_index(1.0, edges) == 4


def test_group_assignment_does_not_split_a_group():
    groups = ["a", "a", "b", "c", "c", "c", "d"]
    assign = _assign_by_group(groups, 0.4, 0.2, 42)
    for name in set(groups):
        labels = {assign[i] for i, group in enumerate(groups) if group == name}
        assert labels <= {"train", "val", "test"}
        assert len(labels) == 1


def test_bootstrap_mean_of_a_constant_is_that_constant():
    out = bootstrap_mean(np.array([0.2, 0.2, 0.2, 0.2]), n_boot=50, seed=42)
    assert math.isclose(out["point"], 0.2, abs_tol=1e-12)
    assert math.isclose(out["ci95"]["low"], 0.2, abs_tol=1e-12)
    assert math.isclose(out["ci95"]["high"], 0.2, abs_tol=1e-12)
    assert out["ci95"]["includes_0"] is False
