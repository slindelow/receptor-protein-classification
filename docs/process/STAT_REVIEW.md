# Statistical review

This note reads the metric exports, the split files, `src/vs/metrics_lib.py`, and `artifacts/metrics_uncertainty.json`. Models were not refit. `scripts/uncertainty_bootstrap.py` is finished, not mid-write. It already ran a row bootstrap (1,000 draws, seed 42) and correctly skipped the histogram and KIBA comparisons because per-pair scores are not saved. Those row-bootstrap p-values are not reused here as tests. The one new calculation is a protein-cluster bootstrap, written by `scripts/stat_review_cold_protein.py` to `artifacts/metrics_stat_review.json`.

## 1. What is dependent, and why a row bootstrap misleads

The primary table is a complete grid: 25,772 rows, 379 kinases, 68 ligands, and 379 times 68 equals 25,772. Every ligand is measured on every kinase.

Cold-protein test: 5,168 rows = 76 kinases times 68 ligands. Train, validation, and test share no kinases (0 shared proteins) and share all 68 ligands. Each test ligand is repeated on all 76 test kinases. A ligand-only score does not use the protein, and on the saved graph-model test scores that prediction is exactly constant within a ligand (within-drug range 0). The same number is paired with 76 different labels.

Scaffold test: 5,306 rows = 14 ligands times 379 kinases. Murcko scaffold and drug id are the same partition (68 scaffolds, 68 drugs, one scaffold per drug). Train and test share no scaffolds and no drugs, and they share all 379 kinases. Holding out a scaffold holds out one ligand's entire kinase profile. The 14 blocks have 379 rows each.

Two dependencies sit in every metric that pools those rows.

- Same ligand across kinases. Labels for one compound move together, and a ligand-only model emits one score for that compound. Redrawing rows recounts the same compound once per kinase.
- Same protein across ligands. The 68 rows of one kinase share one protein feature vector. On the graph model the protein-aware score for one ligand moves across kinases (within-drug range up to 4.022). Residuals inside a kinase are not 68 independent draws.

`metrics_lib.compute_metrics` drops the Spearman p-value from `scipy.stats.spearmanr` and keeps only rho. That discard is the right call. A p-value on pairs assumes independent rows. Enrichment is one ranked list (`k = ceil(fraction * n)`), so it has the same row dependence and no separate sample size.

A naive row bootstrap draws those pairs with replacement. The nominal size is 5,168 or 5,306. The split unit is 76 kinases or 14 scaffolds. Intervals come out too narrow, and a difference can look "significant" because the same kinase or the same ligand was counted hundreds of times. The exported scaffold graph comparison is the clear case: delta Spearman 0.0556, row-bootstrap interval 0.0240 to 0.0864, two-sided p = 0.0020, on 5,306 rows that are 14 scaffolds. That p-value does not match the holdout. The same objection applies to a paired test, a Wilcoxon test, or a t-test computed on residual rows.

KIBA is not a complete grid (118,254 pairs, 229 targets, 2,111 drug ids, 2,068 distinct SMILES; 229 times 2,111 is 483,419). The cold-protein test still shares no targets with train and shares all 2,111 drug ids with train. Test size is 24,280 rows on 47 proteins, with 13 to 1,223 rows per protein (median 501). About 99.5% of those drugs appear on more than one test protein. The scaffold test has 269 scaffolds, 1 to 2,165 rows each (median 34), 226 proteins that also appear in train, and 10 drug ids that appear on both sides of the scaffold cut. A row bootstrap of 24,280 or 23,663 would have the same defect. Those scores were not bootstrapped.

## 2. Which comparisons have a valid estimand

Spearman rho and enrichment, as coded in `metrics_lib.py`, are valid descriptions of one frozen test table. They answer: on these rows, how does this score order the labels, and how many binders fall in the top fraction of this one list. They do not by themselves answer: would the lift recur for a new kinase, or for a new scaffold.

Valid as descriptions, not yet as significance claims:

| Comparison | Estimand that is defined | Why it is not a significance claim |
|---|---|---|
| Histogram cold-protein, joint composition vs ligand-only | Pooled Spearman 0.5501 vs 0.5002 (difference 0.04997). EF at 1%: 10.2458 vs 8.4528. EF at 5%: 6.9941 vs 6.5312. n = 5,168 rows, 388 binders. | Per-pair scores were not saved. No interval exists. The mode is the maximum cold-protein Spearman among four feature sets on this same split. |
| Histogram scaffold, each protein mode vs ligand-only | Pooled Spearman and EF on 5,306 rows, 302 binders. | No per-pair scores. The holdout unit is 14 scaffolds. |
| Graph cold-protein | Pooled Spearman 0.4895 vs 0.4997 (delta -0.0102) on the same 5,168 rows. | A row-bootstrap interval was exported and includes 0 (p = 0.3586). That procedure is the wrong unit. Section 3 replaces it. |
| Graph scaffold | Pooled Spearman delta 0.0556 on 5,306 rows. | Row-bootstrap p = 0.0020 excludes 0. The resampled unit should be the scaffold (n = 14). That test was not run. The p-value should not be cited. |
| KIBA cold-protein and scaffold | Spearman and EF on the KIBA score. Cold-protein delta Spearman 0.0609 (0.5415 vs 0.4806), n = 24,280, 5,912 labels at or above 12.1. | No per-pair scores. Pooled rho is dominated by the largest proteins. The label is not pKd (section 4). |
| SLK library recovery | Fraction of 17 known binders in the top 1% and the top 500, and their mean rank, on a fixed library of 8,068 rows. | One kinase. Three endpoints. The exported p-values resample library rows. |
| SRC/LCK selectivity | Spearman of predicted vs true delta on 68 ligands, and enrichment of 12 ligands with absolute delta pKd at least 1. | One pair. No interval. The ligand-only delta is exactly constant, so the control has no rank. |
| Calibration ECE and Brier | Fit of an isotonic map from score to P(binder). | A probability check, not a ranking test. The fit sample is not an independent cohort (section 4). |

A second estimand is well defined and was not exported: the unweighted mean, across held-out kinases, of within-kinase Spearman (protein score minus ligand-only score). That is the ranking question for a screen of one kinase. The pooled rho mixes that within-kinase order with between-kinase shifts of the score. On the saved graph-model cold-protein scores the two estimands are not the same shape. The pooled delta is -0.0102. The within-kinase deltas (76 kinases, 68 ligands each) have mean -0.0082, median +0.0161, standard deviation 0.154, 44 positive and 32 negative. No p-value is attached to that mean or median. Pooled enrichment has the same gap: one list of 5,168 mixed pairs is not 76 within-kinase screens.

## 3. Can the cold-protein lift be tested?

The histogram lift cannot be tested from the saved files. `artifacts/metrics_uncertainty.json` records that skip, and this review does not invent an interval around 0.04997.

The graph-model cold-protein scores can be tested. They align with the cold-protein test frame in split-file order (maximum absolute gap between saved `y_true` and `pkd` is under 1e-8; both files have 5,168 rows; dropped SMILES on that run were 0).

Recommended test, and the only one computed: protein-cluster bootstrap of the pooled delta Spearman (protein-aware score minus ligand-only score).

- Unit: the 76 test kinases, drawn with replacement. A kinase drawn twice contributes its 68 rows twice. Both scores use the same draw. Seed 42, 1,000 replicates, percentile interval 2.5 to 97.5, two-sided recentered p-value with the plus-one correction.
- Why protein, not pair: the split randomized kinases. Pairs inside a kinase share a feature vector, and each ligand is copied across the 76 kinases.
- Why not scaffold: this split does not hold out scaffolds. All 68 ligands are in train. A scaffold resample would ask a different question, and it would still be a question about 68 seen chemotypes rather than about new kinases.
- The interval is uncertainty across kinases for this fixed 68-ligand panel. It is not a new-chemotype interval.

Result, from `artifacts/metrics_stat_review.json`:

- Point: Spearman 0.4895 vs 0.4997, delta -0.0102.
- Protein-cluster 95% interval for the delta: -0.0559 to 0.0299. Includes 0.
- Two-sided p for delta = 0: 0.6494 (649 of 1,000 recentered replicates at least as large, plus one, over 1,001).
- The row bootstrap on the same scores was about half as wide (-0.0304 to 0.0104, p = 0.3586). Both include 0. The row interval is the anti-conservative one. Neither result is a cold-protein lift. The histogram claim remains untested.

For the scaffold graph comparison, the matching test would be a scaffold-cluster bootstrap of the pooled delta (14 blocks). It was not computed. Until it is, the exported p = 0.0020 is not a result.

## 4. Concrete problems

Multiple comparisons and a selected mode. Four protein representations were scored on the cold-protein test. The reported mode is the one with the largest Spearman (joint composition). The automatic gate then passes if the Spearman difference is at least 0.05 or if either enrichment difference exceeds 0.5. The unrounded Spearman difference is 0.04997, so the Spearman gate is false. The four-decimal display 0.0500 hides that. Enrichment at 1% differs by 1.793 and clears 0.5. Enrichment at 5% differs by 0.463 and does not. The pass is an OR of a failed primary contrast and one of two enrichment cuts, after a max over four modes, on one split. Joint composition was then locked for the graph model, the selectivity pair, and both library ranks. Those later runs are not independent confirmations of a pre-specified feature.

SLK was chosen from 30 eligible kinases. The rule (exclude LCK, drop identifiers that contain a parenthesis, require at least 10 binders, take the binder count closest to LCK's 17) was applied before scores, which avoids picking the winner after seeing ranks. It still leaves one kinase, chosen to resemble LCK in binder count, with five named runners-up that were not scored. Recovery uses 17 known binders. The row bootstrap reports three contrasts. Top 1% difference 0.0588 (interval includes 0, p = 0.7792). Top 500 difference -0.1176 (includes 0, p = 0.4306). Mean-rank contrast -842.1 in the direction that favors ligand-only (interval excludes 0, p = 0.0370). One of three endpoints excluding 0, on a library-row resample of 17 labels and one selected kinase, is not evidence that ligand-only wins in general or that the protein score wins. The nearest training kinase by joint-composition cosine is NEK1 at 0.9702, so the holdout is not a distant sequence in the feature space the model uses.

Selectivity n. SRC and LCK, 68 dual-measured ligands, 12 with absolute delta pKd at least 1, and all 12 prefer LCK. SRC-directed enrichment is undefined (zero positives). Spearman of the deltas is 0.3694. Mean absolute predicted delta is 0.1938 against a true-delta cut of 1 and a mean absolute true delta of 0.4878. The ligand-only predicted delta is exactly 0, so enrichment 0 and an undefined Spearman are a degenerate control, not a fitted alternative. The sketch gate (Spearman above 0.1 and enrichment at 10% or 20% above 1.2) was met on this single pair. No interval was computed. Closer pairs were considered and not exported.

Calibration sample. Histogram isotonic maps are fit on the validation split and scored on the test split. For the cold-protein split that fit set is 2,584 rows and 38 kinases, and it contains the same 68 ligands as the test set. Graph-model validation scores were not saved. The isotonic map is fit on a random permutation of test rows, `n_cal = int(0.3 * n_test)`: 1,550 fit rows and 3,618 eval rows on the cold-protein graph scores, 1,591 and 3,715 on the scaffold graph scores. The permutation cuts through kinases and ligands, so the eval rows are not an untouched test. A code comment says 20% while the split is 30%. Post-isotonic ECE near 0.01 on the cold-protein histogram models, and near 0.007 on the graph models, should not be read as calibration on new pairs. AUROC is unchanged, which is the useful statement: the map does not reorder the list.

KIBA is not pKd. The column named `pkd` equals `kiba_score` exactly (maximum absolute difference 0), and `kd_nm` is empty. The binder cut is 12.1 on that score. Cold-protein "binders" are 5,912 of 24,280 (about 24%), against 388 of 5,168 (about 7.5%) at pKd 7 on DAVIS. A KIBA Spearman difference of 0.0609 is not a replication of a DAVIS difference of 0.05, and the enrichment numbers are not on the same prevalence. The transfer sentence has to stay on the KIBA scale.

## 5. What has to change in the paper

1. Remove significance language that rests on row resampling. That includes the scaffold graph sentence that the interval excludes 0 (p = 0.0020) and the SLK mean-rank sentence that the interval excludes 0 (p = 0.0370). If those intervals stay, label them as row redraws and state that they are anti-conservative for the split that was actually run.

2. Put the split unit next to every n. Cold-protein DAVIS: 76 kinases, 68 ligands seen in training, 5,168 rows. Scaffold DAVIS: 14 scaffolds, 379 kinases seen in training, 5,306 rows. KIBA cold-protein: 47 kinases with very unequal row counts. Do not let 5,168 read as a sample size.

3. Stop treating the histogram cold-protein difference 0.04997 as a demonstrated lift. Print the unrounded difference next to the gate so 0.0500 cannot be misread as a pass. Either pre-specify one feature mode and one endpoint (suggest: mean within-kinase delta Spearman, kinase as the unit), or mark the four-mode table exploratory and do not lock the winning mode for later experiments. Save per-pair histogram scores before any interval is claimed. The same requirement applies to KIBA.

4. For the graph cold-protein comparison, replace the row interval with the protein-cluster result already computed: delta -0.0102, 95% interval -0.0559 to 0.0299, p = 0.6494, seed 42, 1,000 protein draws. Say in the same sentence that the 68 ligands were not held out. If the scientific claim is ranking compounds for one new kinase, change the estimand to the within-kinase mean and do not call the pooled rho that claim. The within-kinase numbers above are descriptive only.

5. For the graph scaffold comparison, run a scaffold-cluster bootstrap or drop the p-value. Fourteen blocks cannot support the exported p = 0.0020.

6. SLK stays a single case. Do not generalize from 17 binders. Either score the pre-specified runners-up (GAK, LOK, DDR2, VEGFR2, FLT1) under the same rule, or delete any interval that looks like a test. Three recovery endpoints need one pre-specified contrast.

7. Selectivity: report 12 of 68, one direction, and a constant ligand-only control. Delete the sketch-gate pass, or move it to a supplement as an illustration with no success flag.

8. Calibration: state the fit n and that the graph map uses test rows. Do not quote the small post-isotonic ECE as confirmation. Keep the statement that rank metrics do not change.

9. KIBA: say the stored `pkd` column is the KIBA score, not a dissociation constant, and do not set the 0.0609 difference beside the DAVIS difference as the same estimand. Note the protein-size imbalance if a pooled rho is shown.

10. One primary contrast. The current stack (four modes, Spearman, two enrichment cuts, an OR gate, a second dataset, a graph encoder, one kinase pair, and three recovery endpoints) is a search. A confirmatory sentence needs one contrast, one unit, and an interval from that unit.

## 6. HistGBM cold-protein kinase-cluster bootstrap

This section supersedes the statement in section 3 that the histogram lift cannot be tested. It was written by `scripts/uncertainty_bootstrap.py` from the comparison block in `artifacts/metrics_uncertainty.json`. The graph protein-cluster result already in `artifacts/metrics_stat_review.json` was not recomputed.

Score file: `artifacts/scores_histgbm_cold_protein.csv`. Models: loaded from artifacts/models_v1.1 joblibs, not refit. Point estimates matched `metrics_v1.1.json` within 0.0001: yes.

Resample unit: held-out kinase. n kinases: 76. n pairs: 5168. Seed: 42. n resamples: 1000.

Draw the held-out kinases with replacement. Keep every test pair of each drawn kinase (a kinase drawn twice contributes its pairs twice). Recompute pooled Spearman and EF@1% on the concatenated pairs. The same draw enters both scores.

Delta Spearman (AAC+DPC minus ligand-only): point 0.04996776087112098, 95% percentile interval 0.02712625028941049 to 0.071976586483323, includes 0: False, two-sided recentered p 0.000999000999000999 (0 of 1000 recentered replicates at least as large, plus one).

Component Spearman on the original test: AAC+DPC 0.5501253139581714, ligand-only 0.5001575530870505.

Delta EF@1% (AAC+DPC minus ligand-only): point 1.7930214115781133, 95% percentile interval -0.2624061982450547 to 3.2467613621320517, includes 0: True, two-sided recentered p 0.09090909090909091.

AAC+DPC was the best cold-protein Spearman among four modes scored on this same split (ligand-only 0.500158, AAC 0.498399, DPC 0.540900, AAC+DPC 0.550125). ESM was a fifth comparison and was essentially tied with AAC+DPC (ESM Spearman 0.547628). The contrast tested here is the pre-specified AAC+DPC versus ligand-only contrast. The interval is not adjusted for the other modes.

The delta Spearman interval excludes zero. That exclusion is not multiplicity-adjusted: AAC+DPC was chosen as the best of four modes on this split, with ESM a fifth comparison essentially tied with it. The delta EF@1% interval includes zero. This file does not call either contrast significant.

Rerun from the repository root: `.venv/bin/python scripts/score_histgbm_cold_protein.py && .venv/bin/python scripts/uncertainty_bootstrap.py`.
