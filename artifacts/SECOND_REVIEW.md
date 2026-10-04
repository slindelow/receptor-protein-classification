# Second scientific review (fresh eyes)

Reviewer role: secondary check of written claims against exported artifacts and metric code. Models were not rerun. Numbers were not invented.

Date: 2026-10-03 (America/Toronto). Project: `.`.

---

## 1. Verdict

Headline ranking numbers in the abstract and Results match the exported metric tables at the paper’s four-decimal rounding, including cold-protein Spearman 0.5501 versus 0.5002 and EF at 1% 10.2458 versus 8.4528 from `artifacts/metrics_v1.1.md`. The paper’s enrichment definition matches `src/vs/metrics_lib.py` (top-slice binder rate over overall binder rate; top-1% library cutoff k = ceil(0.01 × 8068) = 81). Cold-protein and ligand-only comparisons use the same test sizes and binder counts on the saved splits. One clear factual slip: Methods says primary-split training and validation counts “were not exported,” but `artifacts/metrics.json` and `artifacts/metrics_v1.1.json` already carry `n_train` (17813 scaffold, 18020 cold-protein). Several results are point estimates without uncertainty, with sketch gates, one kinase (SLK) chosen from 30 eligible targets, and heavy ligand reuse, so they are not valid to read as statistical significance. For a public repo the training-size sentence, reproduce coverage beyond v1/v1.1/v1.2/v1.4, and `.gitignore` for large raw files need fixing before a stranger can trust and rebuild the freeze.

---

## 2. Confirmed

Claims that match exports (file cited). Rounding is to four decimals unless noted.

### Abstract / sample sizes / manifests

| Claim | Export |
|---|---|
| DAVIS 25,772 pairs, 68 ligands, 379 targets; pKd = 9 − log10(Kd_nM); binder pKd ≥ 7 | `data/MANIFEST.json`; `configs/default.yaml` |
| KIBA 118,254 pairs, 2,068 ligands, 229 targets; binder score ≥ 12.1; TDC 403 then DeepDTA mirror | `data/MANIFEST_kiba.json`; `artifacts/metrics_kiba.md` |
| Library 8,068 molecules, 7,582 scaffolds, HIV subsample seed 42 + 68 DAVIS spike-in | `data/screen_library/MANIFEST_screen.json`; `artifacts/metrics_screen.md` |
| Scaffold test 5,306 pairs / 302 binders; cold-protein test 5,168 / 388 | `artifacts/metrics.md`, `artifacts/metrics_v1.1.md` |
| Cold-protein AAC+DPC Spearman 0.5501 vs ligand-only 0.5002; EF@1% 10.2458 vs 8.4528 | `artifacts/metrics_v1.1.md` (also abstract) |
| Selectivity “positive”; ligand-only difference constant | `artifacts/metrics_selectivity.md` (sketch gate YES; ligand-only Δ constant) |
| LCK warm library recovery above chance; LCK not held out; SLK cold recovery mixed, ligand-only wins broader ranks | `artifacts/metrics_screen.md`; `artifacts/metrics_screen_cold.md` |

### Protein ablation table (`RESEARCH_REPORT.md` Results, scaffold and cold-protein)

All ten rows (ligand-only, AAC, DPC, AAC+DPC, ESM-2 × two splits) match `artifacts/metrics_v1.1.md` / `.json` at four decimals, including feat dims 2048 / 2068 / 2448 / 2468 / 2368. AAC rows match the earlier `artifacts/metrics.md` fit. Auto gate text matches: ΔSpearman exported as 0.0500, `rho_gate=False`, `ef_gate=True`.

### Graph encoder

Chemprop rows match `artifacts/metrics_v1.2.md`. Cold-protein Chemprop+protein vs Chemprop ligand-only: rho 0.4895 vs 0.4997 (Δ −0.0102); EF@1% 8.4528 vs 8.7090. HistGBM baselines cited as not retrained match v1.1.

### KIBA

All eight cells match `artifacts/metrics_kiba.md`. Cold-protein ΔSpearman 0.0609; EF@1% does not improve (3.8872 vs 3.9041).

### Calibration

All ECE / Brier / AUROC / n_eval rows match `artifacts/metrics_calibration.md`. Protocol split (HistGBM isotonic on validation; Chemprop isotonic on 30% of test) matches the export and `scripts/run_calibration.py`. AUROC unchanged by the map is consistent with rank preservation.

### Selectivity (SRC / LCK)

| Claim | Export |
|---|---|
| Train excludes SRC and LCK, 25,636 pairs | 25,772 − 136 = 25,636; `artifacts/metrics_selectivity.json` protocol |
| Best-window identity 0.4008; AAC+DPC cosine 0.9621 | `metrics_selectivity.json` `pair` |
| 68 dual ligands; 12 with \|ΔpKd\| ≥ 1, all prefer LCK | same |
| Spearman(Δ) 0.3694; mean \|Δpred\| 0.1938; EF@5/10/20% 2.83 / 3.24 / 2.43 | `metrics_selectivity.md` (two-decimal EF as exported) |
| Ligand-only Δ constant → Spearman undefined, EF 0 | same |

### Library screens

| Claim | Export |
|---|---|
| LCK fit n_train = 20,466 (scaffold train+val); seq len 509; 17 binders | `metrics_screen.md` / `.json` |
| LCK mean rank 378.4118; random mean 4034.5; top-50/100/500 = 7/7/15; enrichment top-500 ≈ 14.2376×; ranks 1310 and 2412 among binders; rank-1 pred 9.1475 | same |
| SLK hold-out 68 pairs; n_train 25,704; seq len 1,235; nearest NEK1 cosine 0.9702; rule and runners-up GAK/LOK/DDR2/VEGFR2/FLT1 | `metrics_screen_cold.md` / `.json` |
| Protein: 8/17 top 1%, 10/17 top 500, mean rank 1311.2353 | same |
| Ligand-only: 7/17, 12/17, mean rank 469.1765; verdict mixed; deltas +0.0588 / −0.1176 / −842.0588 | same |
| SLK binder true pKd 7.5850 at protein rank 7735 / ligand-only 2367; pred 5.1183; HIV neighbor at protein rank 2 not a known binder | `metrics_screen_cold.md` detail and top-10 |

### Metric definitions vs code

- Spearman: `scipy.stats.spearmanr` in `src/vs/metrics_lib.py` `compute_metrics`.
- EF: `enrichment_factor` = hits in top k / (n_actives × k/n) with `k = max(1, ceil(fraction × n))`. Equivalent to the paper’s “binder rate in the top slice / binder rate in the test list.”
- AUROC: `sklearn.metrics.roc_auc_score` on binder labels vs continuous score.
- ECE: 10 equal-width bins in `scripts/run_calibration.py` `expected_calibration_error`.
- Selectivity EF: same `enrichment_factor` on \|Δpred\| ranks (`scripts/run_selectivity_v1.4.py`), fractions 0.05 / 0.10 / 0.20 from `configs/v1.4_selectivity.yaml`.
- Same rows for protein vs ligand-only on a split: shared split CSVs; identical `n_test` and `n_binders` in exports.

### Notebook drift (`artifacts/RESEARCH_REPORT_NOTEBOOK.md`)

Structure and voice differ (versioned first-person notebook vs formal paper). Headline cold-protein and screen numbers match. The notebook is more precise about training sizes: it says they were not printed in `metrics.md`, not that they were never exported. No numerical drift of abstract headlines between notebook and paper relative to the exports.

---

## 3. Mismatches

### Paper vs export

1. **Primary-split training sizes “not exported”**  
   - Paper (`## Methods` / `### Splits and what was held out`): “Exact training and validation counts for those two primary splits were not exported.”  
   - Export: `artifacts/metrics.json` (and `metrics_v1.1.json`) include `n_train`: scaffold **17813**, cold-protein **18020**. Implied validation ≈ 2653 and 2584 if total pairs = 25772.  
   - `artifacts/metrics.md` tables omit `n_train`, which is likely what was meant; the paper sentence as written is false.

2. **Exported ΔSpearman 0.0500 vs gate failure**  
   - Paper and `metrics_v1.1.md` report Δ = **0.0500** and “strict Spearman part … not met.”  
   - Unrounded JSON: 0.550125… − 0.500158… = **0.049968…**, so `Δ ≥ 0.05` fails.  
   - Not a table typo, but a reader who only sees 0.0500 will not see why the Spearman gate failed.

3. **Selectivity EF decimal places**  
   - Paper / `metrics_selectivity.md`: **2.83, 3.24, 2.43**.  
   - JSON: **2.8333…, 3.2381…, 2.4286…**.  
   - Consistent with the selectivity markdown export, inconsistent with four-decimal EF elsewhere in the paper.

No other headline abstract or Results table cell failed a four-decimal check against `metrics_v1.1`, `metrics_v1.2`, `metrics_kiba`, `metrics_calibration`, `metrics_screen`, or `metrics_screen_cold`.

### Paper vs code

No definition mismatch found for Spearman, EF, AUROC, or ECE. Library top 1% = rank ≤ 81 matches code and `metrics_screen_cold.json` (`top1_k`: 81).

---

## 4. Statistics that are not valid to claim

These are limitations of inference, not export typos. The paper often softens them; they still cannot support strong language.

1. **No confidence intervals or tests** on Spearman, EF, AUROC, recovery fractions, or selectivity EF. Discussion notes missing CI on enrichment; Results still use “beat,” “pass,” “positive,” and “far above chance” without a test. “Significant” as a stats word does not appear; the risk is everyday English that reads like inference.

2. **Pre-stated gates are not hypothesis tests.** Cold-protein pass is OR of ΔSpearman ≥ 0.05 and “clear EF lift.” Selectivity “sketch” gate is ρ > 0.1 and EF@10/20% > 1.2. Those are checklist thresholds, not p-values.

3. **Selectivity n is tiny and one-sided.** Twelve selective ligands out of 68; EF@5% uses k = ceil(0.05 × 68) = **4**, with expected hits ≈ 0.71. EF@5% ≈ 2.83 is a handful of ranks. Zero SRC-preferring ligands at the cut, so directionality was not tested. One kinase pair only.

4. **SLK is one draw from 30 eligible kinases** under a binder-count rule tied to LCK. Runners-up were not scored. A single cold target cannot support a general claim that protein features help (or fail) on cold screens.

5. **Ligand reuse / not cold-ligand.** DAVIS has 68 ligands reused across kinases. Cold-protein and cold-target (SLK) still leave those chemotypes labeled elsewhere. Ligand-only cold-protein Spearman ≈ 0.50 already shows the leak. Abstract and Discussion state this; any claim of “sequence understanding” would be invalid (the paper mostly avoids that claim).

6. **Multiple protein modes** (AAC, DPC, AAC+DPC, ESM-2) compared; AAC+DPC chosen as winner for later experiments. No multiplicity adjustment. Later v1.2–v1.6 inherit that selection.

7. **Warm LCK screen** trains with LCK pairs present (scaffold train+val). Recovery is a sanity check, not prospective targeting. Ligand-only recovery for LCK was **not** exported (`metrics_screen.json` has no `recovery_ligand_only`), so protein vs ligand-only recovery on LCK cannot be claimed from exports.

8. **Chemprop calibration ECE after isotonic on a test slice** is not an unbiased test ECE (paper states this). Do not cite those post-isotonic ECE values as held-out calibration quality.

9. **Graph model bit-stability not claimed**; documented tolerance 0.05. Fine as a negative result for this freeze, not as a definitive Chemprop-vs-HistGBM benchmark.

10. **KIBA binder cut and label** are not pKd; transfer on Spearman only. Valid only as a secondary composite-score check.

---

## 5. Presentation problems

| Issue | Where |
|---|---|
| Methods says primary train/val counts were not exported, while JSON has them | `RESEARCH_REPORT.md` → `### Splits and what was held out` |
| Four-decimal metrics without uncertainty read like high precision | Results tables under protein features, graph, KIBA, calibration, screens |
| Selectivity EF at two decimals next to four-decimal Spearman | `### Selectivity on one kinase pair` table |
| “Clear lift,” “was positive,” “far above chance,” “beat” without a test | Abstract; Results protein features; selectivity; LCK recovery |
| Abstract leads with warm LCK recovery before the cold SLK caveat | `## Abstract` |
| ΔSpearman printed as 0.0500 beside “Spearman gate not met” | Results protein features paragraph on AAC+DPC |
| Calibration table is wide (protocol + six metrics); Chemprop “ECE after” looks excellent unless the protocol column is read | `### Calibration` |
| Notebook is a different document (first person, version headings); a reader opening both may think two papers disagree on style, not on numbers | `artifacts/RESEARCH_REPORT_NOTEBOOK.md` vs `RESEARCH_REPORT.md` |
| README reproduce list stops at v1.4; no one-command path for KIBA, calibration, v1.5, v1.6 | `README.md` → `## Reproduce` |
| `.gitignore` ignores `artifacts/models_v1.1/` and some raw DAVIS files, but not `kiba_pairs.csv` (~98 MB), `moleculenet_hiv.csv`, or other `artifacts/models_*` dirs | `.gitignore` |

Tables are otherwise labeled (split, model, metrics, n test, n binders). Column headers are legible. Limits and non-claims in Discussion are clearer than the Abstract’s density.

---

## 6. What must be fixed before a public repo (ordered)

1. **Correct the training-size sentence** in Methods: either report `n_train` (and validation) from `metrics.json` / `metrics_v1.1.json`, or say they are absent from the markdown tables but present in JSON. Do not leave “were not exported.”

2. **Ship reproduce entry points for every exported freeze the paper cites:** KIBA (`run_kiba_v1.3.py`), calibration (`run_calibration.py`), screen v1.5, cold screen v1.6, or one orchestrating script plus README steps. Today only `reproduce.sh`, `reproduce_v1.1.sh`, `reproduce_v1.2.sh`, `reproduce_v1.4.sh` exist.

3. **Harden `.gitignore` and download docs for publication size and hashes:** ignore or Git-LFS `data/raw/kiba_pairs.csv`, HIV raw, and large model dirs consistently; keep manifests and split ID lists; document fetch-by-URL + SHA256 verify from `MANIFEST.json`, `MANIFEST_kiba.json`, `MANIFEST_screen.json`. KIBA split rebuild from frozen pairs is already noted; make the command explicit.

4. **Tone down inferential wording** in Abstract/Results (“beat,” “positive,” “far above chance”) or attach “point estimate; no CI.” Keep the honest Discussion limits; align the Abstract with them (especially warm LCK vs cold SLK).

5. **Annotate selectivity** as a one-pair sketch with n_selective = 12 and EF@5% over top-4; do not let Abstract “positive” stand without that constraint. Prefer four-decimal EF to match other tables, or state “rounded as in `metrics_selectivity.md`.”

6. **Clarify the 0.0500 / gate-fail** in one sentence (unrounded Δ < 0.05; rounded display 0.0500) so readers do not think the export contradicts the gate.

7. **Optional but useful:** add `n_train` / `n_val` columns to `metrics.md` / `metrics_v1.1.md`; export LCK ligand-only recovery if any warm-screen protein-vs-ligand comparison is desired; add a short “how to read the notebook vs the paper” note so `RESEARCH_REPORT_NOTEBOOK.md` is not mistaken for a second conflicting manuscript.

---

End of review. `RESEARCH_REPORT.md` was not edited.
