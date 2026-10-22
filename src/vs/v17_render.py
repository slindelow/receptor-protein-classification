"""Render the v1.7 README opening, report section, and metrics table from JSON.

The prose quotes formatted values from the metrics payload. It does not invent
a second estimate. Banned characters are rejected before any file is written.
"""
from __future__ import annotations

from pathlib import Path


BANNED = ("\u2013", "\u2014", "/workspace", "/home/", "Mila", "significant", "Significant")

SPLIT_ORDER = ("cold_protein", "seqid_c40", "seqid_c50", "seqid_c70")


def fmt(value, digits: int = 4) -> str:
    if value is None:
        return "not computed"
    return f"{float(value):.{digits}f}"


def fmt_p(value) -> str:
    if value is None:
        return "not computed"
    return f"{float(value):.4f}"


def need(obj: dict, key: str):
    if key not in obj:
        raise SystemExit(f"metrics payload is missing {key}")
    return obj[key]


def zero_phrase(flag) -> str:
    if flag is True:
        return "includes 0"
    if flag is False:
        return "excludes 0"
    return "interval not computed"


def delta_phrase(delta: dict | None) -> str:
    if not delta or delta.get("point") is None:
        return "not computed"
    point = fmt(delta["point"])
    if delta.get("low") is None or delta.get("high") is None:
        return f"{point} (interval not computed)"
    p_bit = ""
    if delta.get("p_value") is not None:
        p_bit = f"; two-sided p = {fmt_p(delta['p_value'])}"
    return (
        f"{point} (95% interval {fmt(delta['low'])} to {fmt(delta['high'])}, "
        f"{zero_phrase(delta.get('includes_0'))}{p_bit})"
    )


def within_phrase(block: dict | None) -> str:
    if not block or block.get("point") is None:
        return "not computed"
    base = f"{fmt(block['point'])} (finite kinases {block.get('n_finite')})"
    ci = block.get("ci95") or {}
    if not block.get("interval_computed") or ci.get("low") is None:
        return base + " (interval not computed)"
    p_bit = ""
    if block.get("p_value") is not None:
        p_bit = f"; two-sided p = {fmt_p(block['p_value'])}"
    return (
        f"{base} (95% interval {fmt(ci['low'])} to {fmt(ci['high'])}, "
        f"{zero_phrase(ci.get('includes_0'))}{p_bit})"
    )


def split_label(name: str, block: dict) -> str:
    if name == "cold_protein":
        return "Cold protein, this file's pair-id order"
    threshold = block.get("threshold")
    if threshold is None:
        return name
    return f"Sequence cluster {float(threshold):.2f}"


def point_relation(new_point, old_point, new_label: str, old_label: str) -> str:
    if new_point is None or old_point is None:
        return "A point comparison is not available."
    new_text = fmt(new_point)
    old_text = fmt(old_point)
    if new_text == old_text:
        return f"The point estimate on {new_label} matches {old_label} at four decimals ({new_text})."
    if float(new_point) < float(old_point):
        direction = "smaller"
    else:
        direction = "larger"
    return (
        f"The point estimate on {new_label} is {new_text}, which is {direction} than "
        f"{old_text} on {old_label}."
    )


def interval_agreement(left: dict | None, right: dict | None, left_name: str, right_name: str) -> str:
    if not left or not right or left.get("includes_0") is None or right.get("includes_0") is None:
        return f"An interval comparison of {left_name} and {right_name} is not available."
    if bool(left["includes_0"]) == bool(right["includes_0"]):
        word = "include 0" if left["includes_0"] else "exclude 0"
        return f"The {left_name} interval and the {right_name} interval both {word}."
    return (
        f"The {left_name} interval {zero_phrase(left['includes_0'])}. "
        f"The {right_name} interval {zero_phrase(right['includes_0'])}."
    )


def assert_clean(text: str, label: str) -> None:
    for token in BANNED:
        if token in text:
            raise SystemExit(f"{label} contains banned token {token!r}")


def _contrast_cells(contrast: dict | None) -> tuple[str, str, str, str, str]:
    if not contrast:
        return ("not computed",) * 5
    return (
        str(contrast.get("n_rows")),
        str(contrast.get("n_kinases")),
        fmt(contrast.get("spearman_a")),
        fmt(contrast.get("spearman_b")),
        delta_phrase(contrast.get("delta_spearman")),
    )


def _md_table(headers: list[str], rows: list[list[str]]) -> str:
    head = "| " + " | ".join(headers) + " |"
    rule = "|" + "|".join(["---"] * len(headers)) + "|"
    body = ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join([head, rule, *body])


def _require_split_keys(name: str, block: dict) -> None:
    required = [
        "full_vs_ligand",
        "uncensored_vs_ligand",
        "within_full_vs_ligand",
        "within_uncensored_vs_ligand",
        "full_vs_additive",
        "additive_vs_ligand",
        "within_full_vs_additive",
        "bins_id_short",
        "bins_fident",
        "identity_association_id_short",
        "identity_association_fident",
        "oracle_spearman_uses_test_labels",
        "ligand_only_max_within_drug_range",
        "n_test_ligands",
        "n_train_ligands",
        "n_test_ligand_rows_missing_from_train",
        "additive_alpha",
        "additive_train_kinase_r2",
    ]
    for key in required:
        if key not in block:
            raise SystemExit(f"{name} is missing {key}")
    if name in ("cold_protein", "seqid_c40"):
        for key in (
            "tobit_vs_tobit_ligand",
            "tobit_uncensored_vs_tobit_ligand",
            "drop_censored_uncensored",
            "drop_censored_full_labels",
        ):
            if key not in block:
                raise SystemExit(f"{name} is missing {key}")
    if name == "seqid_c40" and "no_early_stop" not in block:
        raise SystemExit("seqid_c40 is missing no_early_stop")


def _leak_sentence(name: str, block: dict) -> str:
    n_test = block.get("n_test_kinases")
    if block.get("threshold") is None:
        return (
            "The original cold-protein split holds out target ids. "
            f"It places {n_test} kinases in the test set and leaves sequence neighbors free to sit in training. "
            "That homology path stays open on the original split."
        )
    n_short = block.get("n_test_kinases_with_train_identity_at_least_threshold")
    n_fident = block.get("n_test_kinases_with_train_fident_at_least_threshold")
    threshold = float(block["threshold"])
    short_status = (
        "No test kinase meets that bar, so this definition of the path is closed on this split."
        if n_short == 0
        else "This definition of the path stays open."
    )
    fident_status = (
        "No test kinase meets that bar."
        if n_fident == 0
        else "This unrestricted definition stays open."
    )
    return (
        f"On the {threshold:.2f} cluster split, {n_short} of {n_test} test kinases have "
        f"shorter-sequence identity at least {threshold:.2f} to a training kinase, "
        f"counting only hits with query coverage and target coverage both at least 0.8. {short_status} "
        f"Under unrestricted alignment identity (identical residues divided by alignment length, no coverage filter), "
        f"{n_fident} of {n_test} test kinases meet the same threshold. {fident_status} "
        f"The maximum shorter-sequence identity to the training set is {fmt(block.get('max_of_max_id_short_both80'))} "
        f"(median {fmt(block.get('median_max_id_short_both80'))}). "
        f"The maximum unrestricted alignment identity is {fmt(block.get('max_of_max_fident_any'))} "
        f"(median {fmt(block.get('median_max_fident_any'))}). "
        "The reference pool is the training kinases. Validation kinases are outside the fit and outside that pool."
    )


def _bin_rows(bins: list[dict]) -> list[list[str]]:
    rows = []
    for record in bins:
        contrast = record.get("contrast")
        if contrast is None:
            delta = "fewer than 2 pairs"
            ef = "fewer than 2 pairs"
        else:
            delta = delta_phrase(contrast.get("delta_spearman"))
            ef = delta_phrase(contrast.get("delta_ef"))
        rows.append(
            [
                str(record.get("label")),
                str(record.get("n_kinases")),
                str(record.get("n_pairs")),
                delta,
                ef,
            ]
        )
    return rows


def _bin_zero_sentence(bins: list[dict], definition: str) -> str:
    excluded = []
    included = []
    bare = []
    for record in bins:
        contrast = record.get("contrast") or {}
        delta = contrast.get("delta_spearman") or {}
        label = str(record.get("label"))
        if delta.get("point") is None:
            continue
        if delta.get("includes_0") is True:
            included.append(label)
        elif delta.get("includes_0") is False:
            excluded.append(label)
        else:
            n_kinases = record.get("n_kinases")
            unit = "kinase" if n_kinases == 1 else "kinases"
            bare.append(f"{label} point {fmt(delta.get('point'))} with {n_kinases} {unit}")
    parts = [f"For {definition}:"]
    if excluded:
        parts.append("Intervals that exclude 0: " + ", ".join(excluded) + ".")
    if included:
        parts.append("Intervals that include 0: " + ", ".join(included) + ".")
    if bare:
        parts.append("Bins with a point and no interval, because fewer than 8 kinases were present: " + "; ".join(bare) + ".")
    if len(parts) == 1:
        parts.append("No bin had a Spearman difference.")
    return " ".join(parts)


def _association_sentence(block: dict, key: str, definition: str) -> str:
    assoc = block.get(key) or {}
    return (
        f"Spearman between {definition} and the within-kinase Spearman difference "
        f"(protein features minus ligand-only) is {fmt(assoc.get('spearman'))} "
        f"on {assoc.get('n_kinases')} kinases "
        f"(95% interval {fmt(assoc.get('low'))} to {fmt(assoc.get('high'))}, "
        f"{zero_phrase(assoc.get('includes_0'))}"
        + (f"; two-sided p = {fmt_p(assoc.get('p_value'))}" if assoc.get("p_value") is not None else "")
        + ")."
    )


def _ef_tie_sentence(payload: dict) -> str:
    published = need(payload, "published_cold_protein")
    cold = need(need(payload, "splits"), "cold_protein")["full_vs_ligand"]
    if fmt(published.get("ef_protein")) == fmt(cold.get("ef_a")) and fmt(published.get("ef_ligand_only")) == fmt(cold.get("ef_b")):
        return (
            "Published and pair-id-sorted enrichment point estimates match at four decimals. "
            f"The published enrichment difference is {delta_phrase(published.get('delta_ef'))}."
        )
    return (
        f"Published enrichment at 1% is {fmt(published.get('ef_protein'))} for the protein-feature scores and "
        f"{fmt(published.get('ef_ligand_only'))} for ligand-only. "
        f"After those same scores are sorted by pair_id, protein-feature enrichment is {fmt(cold.get('ef_a'))} "
        f"and ligand-only enrichment is {fmt(cold.get('ef_b'))}. "
        "Ligand-only scores are constant within a drug, so the identities of the rows inside the top 1% can change when tied rows change order. "
        f"The published enrichment difference is {delta_phrase(published.get('delta_ef'))}. "
        f"The pair-id-sorted difference is {delta_phrase(cold.get('delta_ef'))}."
    )


def _exclude_zero_lines(payload: dict) -> list[str]:
    """Contrasts whose interval excludes 0. The primary 40% pooled row is included when it does."""
    lines = []
    for name in SPLIT_ORDER:
        block = need(need(payload, "splits"), name)
        label = split_label(name, block)
        named = [
            ("pooled Spearman", block["full_vs_ligand"].get("delta_spearman"), False),
            ("uncensored Spearman", block["uncensored_vs_ligand"].get("delta_spearman"), False),
            ("enrichment at 1%", block["full_vs_ligand"].get("delta_ef"), False),
            ("within-kinase Spearman", block.get("within_full_vs_ligand"), True),
            ("within-kinase Spearman on uncensored rows", block.get("within_uncensored_vs_ligand"), True),
            ("full model minus additive", block["full_vs_additive"].get("delta_spearman"), False),
            ("additive minus ligand-only", block["additive_vs_ligand"].get("delta_spearman"), False),
        ]
        for key, title in (
            ("tobit_vs_tobit_ligand", "Tobit, full labels"),
            ("tobit_uncensored_vs_tobit_ligand", "Tobit, uncensored rows"),
            ("drop_censored_uncensored", "drop-censored training, uncensored test"),
            ("drop_censored_full_labels", "drop-censored training, full test labels"),
            ("no_early_stop", "early stopping off"),
        ):
            if key in block:
                named.append((title, block[key].get("delta_spearman"), False))
        for title, delta, is_within in named:
            if is_within:
                flag = ((delta or {}).get("ci95") or {}).get("includes_0")
                phrase = within_phrase(delta)
            else:
                flag = (delta or {}).get("includes_0")
                phrase = delta_phrase(delta)
            if flag is False:
                lines.append(f"{label}, {title}: {phrase}")
        for column, bins in (("shorter-sequence identity", block.get("bins_id_short") or []), ("alignment identity", block.get("bins_fident") or [])):
            for record in bins:
                contrast = record.get("contrast") or {}
                delta = contrast.get("delta_spearman") or {}
                if delta.get("includes_0") is False:
                    lines.append(
                        f"{label}, {column} bin {record.get('label')} ({record.get('n_kinases')} kinases): "
                        f"{delta_phrase(delta)}"
                    )
        for key, definition in (
            ("identity_association_id_short", "max shorter-sequence identity to the training set"),
            ("identity_association_fident", "max alignment identity to the training set"),
        ):
            assoc = block.get(key) or {}
            if assoc.get("includes_0") is False:
                lines.append(f"{label}: {_association_sentence(block, key, definition)}")
    return lines


def _change_paragraph(payload: dict) -> str:
    primary = need(need(payload, "splits"), "seqid_c40")
    pooled = primary["full_vs_ligand"]["delta_spearman"]
    within = primary["within_full_vs_ligand"]
    uncensored = primary["uncensored_vs_ligand"]["delta_spearman"]
    enrichment = primary["full_vs_ligand"]["delta_ef"]
    additive = primary["full_vs_additive"]["delta_spearman"]
    sentences = [
        "The answer sentence is the 40% pooled Spearman difference and whether its interval includes 0.",
        "That sentence would change if the same kinase-cluster bootstrap on this same split produced the other include-0 result.",
        "A different feature set, a different model, or a split that was not written in configs/v1.7_reviewer.yaml before the fits would be a different question.",
    ]
    within_flag = (within.get("ci95") or {}).get("includes_0")
    if within_flag is True:
        sentences.append(
            "The within-kinase mean on this split has an interval that includes 0. "
            "The pooled sentence does not say that protein features change ranks inside a kinase."
        )
    elif within_flag is False:
        sentences.append(
            "The within-kinase sentence would change if its interval included 0. "
            f"The current within-kinase difference is {within_phrase(within)}."
        )
    if uncensored.get("includes_0") is True:
        sentences.append(
            "The uncensored-row interval includes 0. The pooled sentence uses the full label vector, including the Kd cap. "
            "It does not extend to the uncensored pairs."
        )
    elif uncensored.get("includes_0") is False:
        sentences.append(
            "The uncensored-row sentence would change if its interval included 0. "
            f"The current uncensored difference is {delta_phrase(uncensored)}."
        )
    if enrichment.get("includes_0") is True:
        sentences.append(
            "The enrichment-at-1% interval includes 0. The answer sentence is about Spearman. "
            f"The enrichment difference is {delta_phrase(enrichment)}."
        )
    elif enrichment.get("includes_0") is False:
        sentences.append(
            "The enrichment sentence would change if its interval included 0. "
            f"The current enrichment difference is {delta_phrase(enrichment)}."
        )
    within_excludes = within_flag is False
    if additive.get("includes_0") is True and within_flag is True and pooled.get("includes_0") is True:
        sentences.append(
            "On this split the pooled difference, the within-kinase difference, and the full-minus-additive "
            "difference all have intervals that include 0."
        )
    elif additive.get("includes_0") is True and within_flag is True:
        sentences.append(
            "Full model minus additive score has an interval that includes 0, and the within-kinase "
            "protein-minus-ligand interval includes 0. The pooled interval is the separate statement above."
        )
    elif additive.get("includes_0") is True and within_excludes:
        sentences.append(
            "Full model minus additive score has an interval that includes 0. The within-kinase "
            "protein-minus-ligand interval excludes 0. A pure kinase shift cannot change within-kinase ranks, "
            "so that within-kinase interval is the result a kinase-level shift does not explain."
        )
    elif additive.get("includes_0") is False:
        sentences.append(
            "The full-model minus additive-score sentence would change if its interval included 0. "
            f"The current difference is {delta_phrase(additive)}."
        )
    old_point = need(payload, "published_cold_protein")["delta_spearman"]["point"]
    if pooled.get("includes_0") is True and pooled.get("point") is not None and float(pooled["point"]) < float(old_point):
        sentences.append(
            "The original cold-protein interval remains the published interval. "
            "A 40% interval that includes 0 means that published difference does not carry over to this sequence-cluster split as a difference separated from 0."
        )
    return " ".join(sentences)


def render_readme(payload: dict) -> str:
    published = need(payload, "published_cold_protein")
    splits = need(payload, "splits")
    primary = need(splits, "seqid_c40")
    _require_split_keys("seqid_c40", primary)
    censoring = need(payload, "censoring")
    composition = need(payload, "composition_uniqueness")
    variance = need(payload, "variance_all_cells")
    variance_u = need(payload, "variance_uncensored_r2")
    cold = need(splits, "cold_protein")
    _require_split_keys("cold_protein", cold)

    pooled = primary["full_vs_ligand"]
    old = published["delta_spearman"]
    lines = []
    lines.append("## Research question")
    lines.append("")
    lines.append(
        "On DAVIS, does a histogram gradient booster that sees Morgan fingerprints plus amino-acid "
        "and dipeptide composition rank pKd higher, by Spearman correlation, than the same booster "
        "on fingerprints alone, when the test kinases are held out of the fit? The reviewer checks "
        "are sequence-cluster holdouts, the Kd = 10000 nM cap, and an additive ligand-plus-kinase score."
    )
    lines.append("")
    lines.append("## Answer")
    lines.append("")
    lines.append(
        f"On the 40% sequence-cluster holdout, the Spearman difference (protein features minus ligand-only) "
        f"is {delta_phrase(pooled['delta_spearman'])}."
    )
    lines.append("")
    lines.append(
        f"The original cold-protein difference, quoted from artifacts/metrics_uncertainty.json, is "
        f"{delta_phrase(old)}. "
        + point_relation(
            pooled["delta_spearman"]["point"],
            old["point"],
            "the 40% split",
            "the original cold-protein split",
        )
        + " "
        + interval_agreement(pooled["delta_spearman"], old, "40% split", "original cold-protein split")
    )
    lines.append("")
    lines.append(
        f"Protein-feature Spearman on the 40% split is {fmt(pooled['spearman_a'])}. "
        f"Ligand-only Spearman on the same rows is {fmt(pooled['spearman_b'])}. "
        f"The enrichment-at-1% difference is {delta_phrase(pooled['delta_ef'])}. "
        f"The mean within-kinase Spearman difference is {within_phrase(primary['within_full_vs_ligand'])}. "
        f"That within-kinase interval resamples kinase-level differences. It is a different bootstrap from the pooled kinase-cluster interval."
    )
    lines.append("")
    sens = []
    for name in ("seqid_c50", "seqid_c70"):
        block = need(splits, name)
        _require_split_keys(name, block)
        sens.append(
            f"{split_label(name, block)} difference {delta_phrase(block['full_vs_ligand']['delta_spearman'])}"
        )
    lines.append(
        "The 0.40 threshold is the primary threshold in configs/v1.7_reviewer.yaml. "
        "The same file specifies 0.50 and 0.70 before the fits. "
        + ". ".join(sens)
        + "."
    )
    lines.append("")
    lines.append(
        f"{censoring['n_censored']} of {censoring['n_pairs']} DAVIS pairs have Kd at the 10000 nM cap "
        f"(fraction {fmt(censoring['fraction_censored'])}). "
        f"Pairs with Kd above the cap: {censoring['n_kd_above_cap']}. "
        f"Uncensored pairs: {censoring['n_uncensored']}. "
        "Kd is right-censored at 10000 nM. pKd is left-censored at 5. "
        f"On uncensored rows of the 40% split, the Spearman difference is "
        f"{delta_phrase(primary['uncensored_vs_ligand']['delta_spearman'])}."
    )
    lines.append("")
    lines.append("## Where a leak could enter")
    lines.append("")
    lines.append(
        f"Kinase identity. Amino-acid plus dipeptide vectors, rounded to 1e-8, are unique for "
        f"{composition['n_unique_aac_dpc_rounded_1e8']} of {composition['n_kinases']} kinases. "
        f"Sequences are unique for {composition['n_unique_sequences']} of {composition['n_kinases']}. "
        "A unique composition vector can identify a kinase that has training rows. "
        "Held-out kinases have no training rows, so an exact lookup of their own labels is closed. "
        "A near-duplicate composition shared with a training kinase is still an open path. "
        "This revision does not hold out composition neighbors. Sequence identity is the neighbor check that was run."
    )
    lines.append("")
    lines.append(
        "Sequence homology. "
        + _leak_sentence("cold_protein", cold)
        + f" Max shorter-sequence identity from a test kinase to the training set is {fmt(cold.get('max_of_max_id_short_both80'))} "
        + f"(median {fmt(cold.get('median_max_id_short_both80'))}). "
        + _association_sentence(cold, "identity_association_id_short", "max shorter-sequence identity to the training set on the original cold-protein test")
        + " "
        + _bin_zero_sentence(cold["bins_id_short"], "the original cold-protein test, shorter-sequence identity")
    )
    lines.append("")
    lines.append(_leak_sentence("seqid_c40", primary))
    lines.append("")
    lines.append(
        "The clusterer is MMseqs2 easy-cluster, cluster-mode 2, seq-id-mode 1, coverage 0.8, cov-mode 0. "
        "The leak counts above come from a separate all-versus-all search with backtrace, using nident "
        "rather than the printed fident column. "
        f"Self-hits: {need(payload, 'mmseqs').get('n_self_hits')}, each with nident equal to the sequence length. "
        "Whole clusters stay together. Assignment uses src/vs/splits.py `_assign_by_group` with seed 42, "
        "test fraction 0.2, and validation fraction 0.1. The function shuffles group ids and fills the test "
        "budget, then the validation budget. A comment in that function mentions a large-first sort. The loop "
        "does not sort by size. This revision uses that function as it stands."
    )
    lines.append("")
    lines.append(_association_sentence(primary, "identity_association_id_short", "max shorter-sequence identity to the training set"))
    lines.append("")
    lines.append(
        f"Ligand scaffolds. The 40% test set has {primary['n_test_ligands']} ligands and the training fit has "
        f"{primary['n_train_ligands']}. "
        f"Test rows whose ligand is absent from training: {primary['n_test_ligand_rows_missing_from_train']}. "
        "These kinase holdouts keep ligands that the ligand-only model has already seen. "
        "The ligand-only score is the control for ligand-driven rank on those ligands. "
        "A split that holds out scaffolds and kinases together was not run. That path stays open. "
        "The original scaffold split remains in the original analysis and is a ligand holdout on its own."
    )
    lines.append("")
    noes = primary["no_early_stop"]
    lines.append(
        "Early-stopping validation. HistGBM early stopping uses an internal 10% of the training rows "
        "(validation_fraction 0.1, n_iter_no_change 15, max_iter 200, random_state 42). "
        "The exported validation split is written and is not that slice. The slice is not built from sequence clusters. "
        "Test kinases are outside the training rows, so their labels are outside the slice. "
        "Direct use of test labels for early stopping is closed. "
        "Dependence on a random training slice, and the mismatch with the exported validation clusters, stays open. "
        f"With early stopping on, the 40% Spearman difference is {delta_phrase(pooled['delta_spearman'])}. "
        f"With early stopping off and the same 200-tree cap, it is {delta_phrase(noes['delta_spearman'])}. "
        + interval_agreement(pooled["delta_spearman"], noes["delta_spearman"], "early-stopping", "no-early-stopping")
        + f" Ligand-only predictions on the 40% split have maximum within-drug range "
        f"{fmt(primary['ligand_only_max_within_drug_range'])}."
    )
    lines.append("")
    lines.append("## What was specified before the fits")
    lines.append("")
    lines.append(
        "Documented in DESIGN.md on 2026-10-02: Spearman and enrichment as ranking estimands, a scaffold split "
        "and a cold-protein split, a ligand-only control, and an ECFP plus gradient-boosting baseline. "
        "The binder cut pKd >= 7 is in configs/default.yaml and in artifacts/DECISIONS.md. "
        "DESIGN.md does not name amino-acid composition or dipeptide composition."
    )
    lines.append("")
    lines.append(
        "Documented with the results, in artifacts/DECISIONS.md: amino-acid composition as the first protein feature, "
        "then a v1.1 ablation. AAC+DPC is the mode that had the best cheap cold-protein Spearman in that ablation. "
        "scripts/uncertainty_bootstrap.py labels the AAC+DPC contrast pre-specified, and the same selection note "
        "says AAC+DPC was the best of the modes scored on that split. This revision treats the AAC+DPC choice as "
        "post hoc relative to DESIGN.md. The 0.05 Spearman-or-enrichment gate is in scripts/run_protein_ablation.py "
        "and in the decision log. It is not in DESIGN.md. The kinase-cluster bootstrap was added in "
        "scripts/uncertainty_bootstrap.py. It is not in the 2026-10-02 design lock. "
        "The published interval is not adjusted for the feature search that selected AAC+DPC."
    )
    lines.append("")
    lines.append(
        "Documented in configs/v1.7_reviewer.yaml before the fits of this revision: thresholds 0.40, 0.50, and 0.70, "
        "with 0.40 primary; identity bin edges 0.00, 0.30, 0.40, 0.50, 0.70, and 1.00; the censoring rule; "
        "four Tobit fits and a drop-censored refit on the original cold split and the 40% split; the additive "
        "baseline; within-kinase Spearman; and early stopping turned off on the 40% split. "
        "Those choices are pre-specified for this revision. They are not part of the 2026-10-02 lock. "
        "Seed 42 and 1000 bootstrap draws are in that file. "
        f"The saved-score recompute of the published Spearman difference and its lower bound matches within 1e-9: "
        f"{published['recompute_matches']} "
        f"(absolute point difference {published['recompute_abs_diff_point']}, "
        f"absolute lower-bound difference {published['recompute_abs_diff_low']})."
    )
    lines.append("")
    lines.append("## Effect sizes")
    lines.append("")
    lines.append(
        "Intervals are kinase-cluster percentile intervals unless a row says the interval was not computed. "
        "The original row is the published file, in the saved score file's row order. "
        "The cold-protein row from this revision is a bootstrap after sorting by pair_id, so its interval can differ "
        "even when the point matches. Four-decimal values are formatted from artifacts/metrics_v1.7.json."
    )
    lines.append("")
    headers = [
        "Split",
        "Pairs",
        "Kinases",
        "Spearman, protein features",
        "Spearman, ligand-only",
        "Delta Spearman",
        "Delta enrichment at 1%",
    ]
    rows = [[
        "Original cold protein, published",
        str(published["n_rows"]),
        str(published["n_kinases"]),
        fmt(published["spearman_protein"]),
        fmt(published["spearman_ligand_only"]),
        delta_phrase(published["delta_spearman"]),
        delta_phrase(published["delta_ef"]),
    ]]
    for name in SPLIT_ORDER:
        block = need(splits, name)
        contrast = block["full_vs_ligand"]
        cells = _contrast_cells(contrast)
        rows.append([split_label(name, block), *cells, delta_phrase(contrast.get("delta_ef"))])
    lines.append(_md_table(headers, rows))
    lines.append("")
    lines.append(_ef_tie_sentence(payload))
    lines.append("")
    lines.append("Uncensored rows, same splits:")
    lines.append("")
    urows = []
    for name in SPLIT_ORDER:
        block = need(splits, name)
        contrast = block["uncensored_vs_ligand"]
        cells = _contrast_cells(contrast)
        urows.append([split_label(name, block), *cells, delta_phrase(contrast.get("delta_ef"))])
    lines.append(_md_table(headers, urows))
    lines.append("")
    lines.append(
        f"Mean within-kinase Spearman difference on the 40% split, full labels: "
        f"{within_phrase(primary['within_full_vs_ligand'])}. "
        f"On uncensored rows: {within_phrase(primary['within_uncensored_vs_ligand'])}."
    )
    lines.append("")
    excluded = _exclude_zero_lines(payload)
    if not excluded:
        lines.append("None of the scanned contrasts have an interval that excludes 0.")
    else:
        lines.append(
            "Contrasts whose interval excludes 0 are listed here with the estimate. "
            "The primary 40% pooled Spearman interval is the answer sentence. "
            "A row in this list does not replace that sentence."
        )
        lines.append("")
        for line in excluded:
            lines.append(f"- {line}")
    lines.append("")
    lines.append("## Alternative explanations that were checked")
    lines.append("")
    lines.append(
        f"Ligand memorization. On the 40% split the ligand-only Spearman is {fmt(pooled['spearman_b'])} "
        f"and the protein-feature Spearman is {fmt(pooled['spearman_a'])}. "
        f"The difference is {delta_phrase(pooled['delta_spearman'])}. "
        "The ligand-only model is constant within a drug on this fit, up to the within-drug range above."
    )
    lines.append("")
    lines.append(
        "Homolog leakage. The original cold-protein split leaves it open. The 40% cluster split and the "
        "identity bins are the check. Bin rows are in the report section. "
        + _bin_zero_sentence(primary["bins_id_short"], "nident / shorter length, both coverages at least 0.8")
    )
    lines.append("")
    lines.append(
        f"Censoring mass point. The cap covers fraction {fmt(censoring['fraction_censored'])} of pairs. "
        f"The uncensored 40% difference is {delta_phrase(primary['uncensored_vs_ligand']['delta_spearman'])}. "
        f"The secondary Tobit model, four Gaussian left-censoring updates at pKd 5, differs from its ligand-only "
        f"Tobit fit by {delta_phrase(primary['tobit_vs_tobit_ligand']['delta_spearman'])} on the full 40% test "
        f"and by {delta_phrase(primary['tobit_uncensored_vs_tobit_ligand']['delta_spearman'])} on uncensored test rows. "
        f"A model fit only on uncensored training rows differs by "
        f"{delta_phrase(primary['drop_censored_uncensored']['delta_spearman'])} on uncensored test rows "
        f"and by {delta_phrase(primary['drop_censored_full_labels']['delta_spearman'])} on the full test labels. "
        "The drop-censored fit is secondary. Where its interval excludes 0 and the primary 40% interval includes 0, "
        "the drop-censored fit does not replace the primary interval."
    )
    lines.append("")
    lines.append(
        f"Kinase main effect versus an interaction. On the complete DAVIS matrix the balanced two-way shares are "
        f"ligand {fmt(variance['share_ligand'])}, kinase {fmt(variance['share_kinase'])}, "
        f"and residual {fmt(variance['share_residual'])}. "
        "The residual share is interaction plus measurement error. One observation per cell cannot separate those two. "
        f"On uncensored rows, least-squares R^2 is ligand {fmt(variance_u['r2_ligand'])}, "
        f"kinase {fmt(variance_u['r2_kinase'])}, and ligand plus kinase {fmt(variance_u['r2_ligand_plus_kinase'])}. "
        f"The share outside that additive fit is {fmt(variance_u['share_outside_additive'])}. "
        f"On the original cold-protein scores, additive minus ligand-only is {delta_phrase(cold['additive_vs_ligand']['delta_spearman'])}. "
        f"Full model minus additive is {delta_phrase(cold['full_vs_additive']['delta_spearman'])}. "
        f"Mean within-kinase protein-minus-ligand is {within_phrase(cold['within_full_vs_ligand'])}. "
        f"On the 40% split, additive minus ligand-only is {delta_phrase(primary['additive_vs_ligand']['delta_spearman'])}. "
        f"Full model minus additive is {delta_phrase(primary['full_vs_additive']['delta_spearman'])}. "
        f"Mean within-kinase full-minus-additive is {within_phrase(primary['within_full_vs_additive'])}. "
        f"An oracle additive score that uses each test kinase's own labels has Spearman "
        f"{fmt(primary['oracle_spearman_uses_test_labels'])}. That number uses test labels."
    )
    lines.append("")
    lines.append(
        "Post hoc feature selection. AAC+DPC was selected on the original cold-protein split after the v1.1 modes "
        "were scored. The sequence-cluster thresholds were not used to pick the feature. "
        "The intervals in this revision are not multiplicity-adjusted for that earlier search."
    )
    lines.append("")
    lines.append(
        "Early-stopping slice. The no-early-stopping fit on the 40% split is the check, reported above. "
        "KIBA was not refit. Chemprop was not refit. A composition-neighbor holdout was not run."
    )
    lines.append("")
    lines.append("## What would change the answer")
    lines.append("")
    lines.append(_change_paragraph(payload))
    lines.append("")
    lines.append("## Limits of this revision")
    lines.append("")
    lines.append(
        "DAVIS is one complete matrix of 68 ligands by 379 kinases. "
        "The additive scaler is fit on all training kinases before the ridge penalty is chosen. Test kinases are outside that fit. "
        f"The penalty on the 40% training kinases is {fmt(primary['additive_alpha'])}, "
        f"with training-kinase R^2 {fmt(primary['additive_train_kinase_r2'])}. "
        "The Tobit fit assumes Gaussian residuals and is labeled secondary. "
        "Tree counts are in artifacts/metrics_v1.7.json. "
        + _tree_sentence(payload.get("fits") or [])
        + " The original analysis, including its published interval, is unchanged."
    )
    return "\n".join(lines).rstrip() + "\n"


def _tree_sentence(fits: list) -> str:
    iters = [row.get("n_iter") for row in fits if row.get("n_iter") is not None]
    if not iters:
        return "No fit log was recorded."
    reached = sum(int(value) == 200 for value in iters)
    return f"{reached} of {len(iters)} recorded fits reached 200 trees."


def render_report(payload: dict) -> str:
    readme = render_readme(payload).rstrip()
    splits = need(payload, "splits")
    parts = [
        "## Reviewer revision (2026-10-22)",
        "",
        "This section is rendered from artifacts/metrics_v1.7.json. The sections above it are the original analysis.",
        "",
        readme,
        "",
        "### Split census",
        "",
    ]
    census_rows = []
    for name in SPLIT_ORDER:
        block = need(splits, name)
        census_rows.append(
            [
                split_label(name, block),
                str(block.get("n_train_pairs")),
                str(block.get("n_test_pairs")),
                str(block.get("n_train_kinases")),
                str(block.get("n_test_kinases")),
                str(block.get("n_train_clusters")),
                str(block.get("n_test_clusters")),
                str(block.get("n_clusters")),
                str(block.get("largest_cluster")),
                str(block.get("n_test_kinases_with_train_identity_at_least_threshold")),
                str(block.get("n_test_kinases_with_train_fident_at_least_threshold")),
            ]
        )
    parts.append(
        _md_table(
            [
                "Split",
                "Train pairs",
                "Test pairs",
                "Train kinases",
                "Test kinases",
                "Train clusters",
                "Test clusters",
                "Clusters in train or test",
                "Largest cluster",
                "Test kinases at threshold, short id",
                "Test kinases at threshold, alignment id",
            ],
            census_rows,
        )
    )
    parts.append("")
    parts.append("### Identity bins on the original cold-protein test")
    parts.append("")
    parts.append("These bins use the saved scores. The reference pool is the training kinases.")
    parts.append("")
    cold_block = need(splits, "cold_protein")
    parts.append(
        _md_table(
            ["Bin", "Kinases", "Pairs", "Delta Spearman", "Delta enrichment at 1%"],
            _bin_rows(cold_block["bins_id_short"]),
        )
    )
    parts.append("")
    parts.append(_bin_zero_sentence(cold_block["bins_id_short"], "original cold-protein shorter-sequence identity"))
    parts.append("")
    parts.append("### Identity bins on the 40% split")
    parts.append("")
    parts.append("Shorter-sequence identity, both coverages at least 0.8, maximum against the training set.")
    parts.append("")
    primary = need(splits, "seqid_c40")
    parts.append(
        _md_table(
            ["Bin", "Kinases", "Pairs", "Delta Spearman", "Delta enrichment at 1%"],
            _bin_rows(primary["bins_id_short"]),
        )
    )
    parts.append("")
    parts.append("Unrestricted alignment identity, nident / alignment length, no coverage filter.")
    parts.append("")
    parts.append(
        _md_table(
            ["Bin", "Kinases", "Pairs", "Delta Spearman", "Delta enrichment at 1%"],
            _bin_rows(primary["bins_fident"]),
        )
    )
    parts.append("")
    parts.append(_bin_zero_sentence(primary["bins_fident"], "unrestricted alignment identity"))
    parts.append("")
    parts.append(_association_sentence(primary, "identity_association_fident", "max unrestricted alignment identity to the training set"))
    parts.append("")
    parts.append("### Additive, Tobit, and drop-censored fits")
    parts.append("")
    extra_rows = []
    for name in SPLIT_ORDER:
        block = need(splits, name)
        label = split_label(name, block)
        extra_rows.append([label, "Full minus additive", delta_phrase(block["full_vs_additive"]["delta_spearman"]), within_phrase(block["within_full_vs_additive"])])
        extra_rows.append([label, "Additive minus ligand-only", delta_phrase(block["additive_vs_ligand"]["delta_spearman"]), ""])
        extra_rows.append([
            label,
            "Oracle additive Spearman, uses test labels",
            fmt(block["oracle_spearman_uses_test_labels"]),
            "",
        ])
        if "tobit_vs_tobit_ligand" in block:
            extra_rows.append([label, "Tobit, full labels", delta_phrase(block["tobit_vs_tobit_ligand"]["delta_spearman"]), ""])
            extra_rows.append([label, "Tobit, uncensored rows", delta_phrase(block["tobit_uncensored_vs_tobit_ligand"]["delta_spearman"]), ""])
            extra_rows.append([label, "Drop-censored training, uncensored test", delta_phrase(block["drop_censored_uncensored"]["delta_spearman"]), ""])
            extra_rows.append([label, "Drop-censored training, full test labels", delta_phrase(block["drop_censored_full_labels"]["delta_spearman"]), ""])
        if "no_early_stop" in block:
            extra_rows.append([label, "Early stopping off", delta_phrase(block["no_early_stop"]["delta_spearman"]), ""])
    parts.append(_md_table(["Split", "Contrast", "Spearman difference or Spearman", "Within-kinase note"], extra_rows))
    parts.append("")
    cold = need(splits, "cold_protein")
    published = need(payload, "published_cold_protein")
    parts.append(
        "Cold-protein protein-feature scores in this file are the saved scores in "
        "artifacts/scores_histgbm_cold_protein.csv. "
        f"Published delta Spearman {delta_phrase(published['delta_spearman'])}. "
        f"Pair-id-sorted recompute in this file {delta_phrase(cold['full_vs_ligand']['delta_spearman'])}. "
        "The original analysis keeps the published interval."
    )
    parts.append("")
    mmseqs = need(payload, "mmseqs")
    parts.append("### Alignment and model notes")
    parts.append("")
    parts.append(str(mmseqs.get("identity_for_bins")))
    parts.append("")
    parts.append(
        f"MMseqs2 release {mmseqs.get('release')}, version string contains the pinned commit, "
        f"cluster mode {mmseqs.get('cluster_mode')}, seq-id mode {mmseqs.get('seq_id_mode')}, "
        f"coverage {mmseqs.get('coverage')}, cov-mode {mmseqs.get('cov_mode')}. "
        f"Alignment file {mmseqs.get('alignment_file')} sha256 {mmseqs.get('alignment_sha256')}."
    )
    parts.append("")
    parts.append(
        "The additive kinase effect for an unseen kinase is a ridge regression on standardized amino-acid "
        "plus dipeptide composition, fit to training-kinase mean residuals after a training ligand effect. "
        "The scaler is fit on every training kinase before RidgeCV chooses the penalty. "
        "An unseen ligand would contribute ligand effect 0. "
        f"Missing ligand rows on the 40% test: {primary['n_test_ligand_rows_missing_from_train']}."
    )
    parts.append("")
    parts.append(_tree_sentence(payload.get("fits") or []))
    return "\n".join(parts).rstrip() + "\n"


def render_metrics_md(payload: dict) -> str:
    body = render_report(payload).rstrip()
    header = [
        "# v1.7 metrics",
        "",
        f"Freeze {payload.get('freeze_id')}. Data freeze {payload.get('data_freeze_id')}. "
        f"Seed {payload.get('seed')}. Bootstrap draws {payload.get('n_bootstrap')}. "
        f"Config {payload.get('config')} sha256 {payload.get('config_sha256')}.",
        "",
        "Formatted from artifacts/metrics_v1.7.json. Unrounded values live in that file.",
        "",
    ]
    return "\n".join(header) + body


def replace_span(path: Path, begin: str, end: str, body: str) -> None:
    text = path.read_text()
    start = text.find(begin)
    stop = text.find(end)
    if start < 0 or stop < 0 or stop < start:
        raise SystemExit(f"markers missing in {path.name}")
    if begin in body or end in body:
        raise SystemExit("rendered body contains a marker")
    updated = text[: start + len(begin)] + "\n" + body.rstrip() + "\n" + text[stop:]
    path.write_text(updated)


def write_documents(payload: dict, root: Path) -> None:
    readme = render_readme(payload)
    report = render_report(payload)
    metrics = render_metrics_md(payload)
    assert_clean(readme, "README opening")
    assert_clean(report, "report section")
    assert_clean(metrics, "metrics markdown")
    art = root / "artifacts"
    v17 = art / "v1.7"
    v17.mkdir(parents=True, exist_ok=True)
    (art / "metrics_v1.7.md").write_text(metrics)
    (v17 / "readme_opening.md").write_text(readme)
    (v17 / "report_section.md").write_text(report)
    replace_span(root / "README.md", "<!-- V17_README_BEGIN -->", "<!-- V17_README_END -->", readme)
    replace_span(
        root / "RESEARCH_REPORT.md",
        "<!-- V17_RESULTS_BEGIN -->",
        "<!-- V17_RESULTS_END -->",
        report,
    )
