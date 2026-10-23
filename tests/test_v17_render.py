"""Renderer checks. Numbers below are fixtures, not study results."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from vs.v17_render import assert_clean, render_readme


def _delta(point, low, high, includes, p_value=0.05):
    return {"point": point, "low": low, "high": high, "includes_0": includes, "p_value": p_value}


def _contrast(point=0.01):
    return {
        "n_rows": 100,
        "n_kinases": 10,
        "spearman_a": 0.2,
        "spearman_b": 0.19,
        "delta_spearman": _delta(point, -0.01, 0.03, True, 0.2),
        "ef_a": 1.0,
        "ef_b": 1.0,
        "delta_ef": _delta(0.0, -0.2, 0.2, True, 0.5),
        "interval_computed": True,
    }


def _within():
    return {
        "point": 0.0,
        "n_finite": 10,
        "n_kinases": 10,
        "ci95": {"low": -0.1, "high": 0.1, "includes_0": True},
        "p_value": 0.4,
        "interval_computed": True,
    }


def _split(threshold):
    block = {
        "threshold": threshold,
        "n_train_pairs": 1000,
        "n_test_pairs": 100,
        "n_train_kinases": 20,
        "n_test_kinases": 10,
        "n_train_ligands": 68,
        "n_test_ligands": 68,
        "n_train_clusters": 15,
        "n_test_clusters": 8,
        "n_clusters": 30,
        "largest_cluster": 4,
        "n_test_kinases_with_train_identity_at_least_threshold": 0 if threshold is not None else None,
        "n_test_kinases_with_train_fident_at_least_threshold": 2 if threshold is not None else None,
        "median_max_id_short_both80": 0.1,
        "max_of_max_id_short_both80": 0.2,
        "median_max_fident_any": 0.3,
        "max_of_max_fident_any": 0.55,
        "full_vs_ligand": _contrast(0.01),
        "uncensored_vs_ligand": _contrast(0.0),
        "within_full_vs_ligand": _within(),
        "within_uncensored_vs_ligand": _within(),
        "full_vs_additive": _contrast(0.0),
        "additive_vs_ligand": _contrast(0.01),
        "within_full_vs_additive": _within(),
        "bins_id_short": [
            {"label": "[0.00, 0.30)", "n_kinases": 10, "n_pairs": 100, "contrast": _contrast(0.01)}
        ],
        "bins_fident": [
            {"label": "[0.40, 0.50)", "n_kinases": 3, "n_pairs": 20, "contrast": _contrast(0.02)}
        ],
        "identity_association_id_short": {
            "n_kinases": 10,
            "spearman": 0.1,
            "low": -0.2,
            "high": 0.4,
            "includes_0": True,
            "p_value": 0.3,
        },
        "identity_association_fident": {
            "n_kinases": 10,
            "spearman": 0.1,
            "low": -0.2,
            "high": 0.4,
            "includes_0": True,
            "p_value": 0.3,
        },
        "oracle_spearman_uses_test_labels": 0.4,
        "ligand_only_max_within_drug_range": 0.0,
        "n_test_ligand_rows_missing_from_train": 0,
        "additive_alpha": 10.0,
        "additive_train_kinase_r2": 0.2,
        "tobit_vs_tobit_ligand": _contrast(0.01),
        "tobit_uncensored_vs_tobit_ligand": _contrast(0.0),
        "drop_censored_uncensored": _contrast(0.0),
        "drop_censored_full_labels": _contrast(0.0),
    }
    if threshold == 0.4:
        block["no_early_stop"] = _contrast(0.01)
    return block


def _payload():
    cold = _split(None)
    return {
        "freeze_id": "fixture",
        "data_freeze_id": "fixture-data",
        "seed": 42,
        "n_bootstrap": 1000,
        "config": "configs/v1.7_reviewer.yaml",
        "config_sha256": "abc",
        "published_cold_protein": {
            "recompute_matches": True,
            "recompute_abs_diff_point": 0.0,
            "recompute_abs_diff_low": 0.0,
            "spearman_protein": 0.55,
            "spearman_ligand_only": 0.5,
            "delta_spearman": _delta(0.05, 0.02, 0.07, False, 0.001),
            "delta_ef": _delta(1.0, -0.2, 3.0, True, 0.09),
            "ef_protein": 1.0,
            "ef_ligand_only": 1.0,
            "n_rows": 5168,
            "n_kinases": 76,
        },
        "censoring": {
            "n_censored": 18000,
            "n_pairs": 25000,
            "fraction_censored": 0.72,
            "n_kd_above_cap": 0,
            "n_uncensored": 7000,
        },
        "composition_uniqueness": {
            "n_unique_aac_dpc_rounded_1e8": 379,
            "n_kinases": 379,
            "n_unique_sequences": 379,
        },
        "variance_all_cells": {"share_ligand": 0.5, "share_kinase": 0.2, "share_residual": 0.3},
        "variance_uncensored_r2": {
            "r2_ligand": 0.4,
            "r2_kinase": 0.1,
            "r2_ligand_plus_kinase": 0.5,
            "share_outside_additive": 0.5,
        },
        "mmseqs": {
            "n_self_hits": 379,
            "release": "17-b804f",
            "cluster_mode": 2,
            "seq_id_mode": 1,
            "coverage": 0.8,
            "cov_mode": 0,
            "alignment_file": "artifacts/v1.7/mmseqs_pairs.tsv",
            "alignment_sha256": "abc",
            "identity_for_bins": "Shorter-sequence identity is nident over the shorter length.",
        },
        "fits": [{"n_iter": 200}, {"n_iter": 50}],
        "splits": {
            "cold_protein": cold,
            "seqid_c40": _split(0.4),
            "seqid_c50": _split(0.5),
            "seqid_c70": _split(0.7),
        },
    }


def test_readme_renders_the_fixture_interval_and_stays_clean():
    text = render_readme(_payload())
    assert_clean(text, "fixture")
    assert "0.0100 (95% interval -0.0100 to 0.0300, includes 0; two-sided p = 0.2000)" in text
    assert "0.0500 (95% interval 0.0200 to 0.0700, excludes 0; two-sided p = 0.0010)" in text
    assert "post hoc" in text
    assert "DESIGN.md" in text
    assert "stays open" in text


def test_banned_token_is_rejected():
    try:
        assert_clean("interval 0.1 \u2013 0.2", "fixture")
    except SystemExit:
        return
    raise AssertionError("en dash was accepted")
