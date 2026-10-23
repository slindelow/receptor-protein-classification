# Ranking Public Kinase-Ligand Affinities With and Without Protein Features

Sofia Lindelow. Original analysis 2026-10-03. Reviewer revision 2026-10-22.

<!-- V17_README_BEGIN -->
## Research question

On DAVIS, does a histogram gradient booster that sees Morgan fingerprints plus amino-acid and dipeptide composition rank pKd higher, by Spearman correlation, than the same booster on fingerprints alone, when the test kinases are held out of the fit? The reviewer checks are sequence-cluster holdouts, the Kd = 10000 nM cap, and an additive ligand-plus-kinase score.

## Answer

On the 40% sequence-cluster holdout, the Spearman difference (protein features minus ligand-only) is -0.0024 (95% interval -0.0277 to 0.0236, includes 0; two-sided p = 0.8671).

The original cold-protein difference, quoted from artifacts/metrics_uncertainty.json, is 0.0500 (95% interval 0.0271 to 0.0720, excludes 0; two-sided p = 0.0010). The point estimate on the 40% split is -0.0024, which is smaller than 0.0500 on the original cold-protein split. The 40% split interval includes 0. The original cold-protein split interval excludes 0.

Protein-feature Spearman on the 40% split is 0.4987. Ligand-only Spearman on the same rows is 0.5011. The enrichment-at-1% difference is -1.7485 (95% interval -2.6863 to 1.2587, includes 0; two-sided p = 0.3287). The mean within-kinase Spearman difference is -0.0058 (finite kinases 77) (95% interval -0.0170 to 0.0063, includes 0; two-sided p = 0.3526). That within-kinase interval resamples kinase-level differences. It is a different bootstrap from the pooled kinase-cluster interval.

The 0.40 threshold is the primary threshold in configs/v1.7_reviewer.yaml. The same file specifies 0.50 and 0.70 before the fits. Sequence cluster 0.50 difference 0.0128 (95% interval -0.0147 to 0.0427, includes 0; two-sided p = 0.3756). Sequence cluster 0.70 difference 0.0090 (95% interval -0.0202 to 0.0391, includes 0; two-sided p = 0.5584).

18343 of 25772 DAVIS pairs have Kd at the 10000 nM cap (fraction 0.7117). Pairs with Kd above the cap: 0. Uncensored pairs: 7429. Kd is right-censored at 10000 nM. pKd is left-censored at 5. On uncensored rows of the 40% split, the Spearman difference is -0.0039 (95% interval -0.0404 to 0.0379, includes 0; two-sided p = 0.8492).

## Where a leak could enter

Kinase identity. Amino-acid plus dipeptide vectors, rounded to 1e-8, are unique for 379 of 379 kinases. Sequences are unique for 379 of 379. A unique composition vector can identify a kinase that has training rows. Held-out kinases have no training rows, so an exact lookup of their own labels is closed. A near-duplicate composition shared with a training kinase is still an open path. This revision does not hold out composition neighbors. Sequence identity is the neighbor check that was run.

Sequence homology. The original cold-protein split holds out target ids. It places 76 kinases in the test set and leaves sequence neighbors free to sit in training. That homology path stays open on the original split. Max shorter-sequence identity from a test kinase to the training set is 0.9087 (median 0.1568). Spearman between max shorter-sequence identity to the training set on the original cold-protein test and the within-kinase Spearman difference (protein features minus ligand-only) is 0.2734 on 76 kinases (95% interval 0.0557 to 0.4841, excludes 0; two-sided p = 0.0130). For the original cold-protein test, shorter-sequence identity: Intervals that exclude 0: [0.50, 0.70), [0.70, 1.00]. Intervals that include 0: [0.00, 0.30). Bins with a point and no interval, because fewer than 8 kinases were present: [0.30, 0.40) point -0.0453 with 4 kinases; [0.40, 0.50) point 0.0577 with 7 kinases.

On the 0.40 cluster split, 2 of 77 test kinases have shorter-sequence identity at least 0.40 to a training kinase, counting only hits with query coverage and target coverage both at least 0.8. This definition of the path stays open. Under unrestricted alignment identity (identical residues divided by alignment length, no coverage filter), 77 of 77 test kinases meet the same threshold. This unrestricted definition stays open. The maximum shorter-sequence identity to the training set is 0.6384 (median 0.0000). The maximum unrestricted alignment identity is 0.8000 (median 0.5238). The reference pool is the training kinases. Validation kinases are outside the fit and outside that pool.

The clusterer is MMseqs2 easy-cluster, cluster-mode 2, seq-id-mode 1, coverage 0.8, cov-mode 0. The leak counts above come from a separate all-versus-all search with backtrace, using nident rather than the printed fident column. Self-hits: 379, each with nident equal to the sequence length. Whole clusters stay together. Assignment uses src/vs/splits.py `_assign_by_group` with seed 42, test fraction 0.2, and validation fraction 0.1. The function shuffles group ids and fills the test budget, then the validation budget. A comment in that function mentions a large-first sort. The loop does not sort by size. This revision uses that function as it stands.

Spearman between max shorter-sequence identity to the training set and the within-kinase Spearman difference (protein features minus ligand-only) is 0.1523 on 77 kinases (95% interval -0.1066 to 0.3878, includes 0; two-sided p = 0.2567).

Ligand scaffolds. The 40% test set has 68 ligands and the training fit has 68. Test rows whose ligand is absent from training: 0. These kinase holdouts keep ligands that the ligand-only model has already seen. The ligand-only score is the control for ligand-driven rank on those ligands. A split that holds out scaffolds and kinases together was not run. That path stays open. The original scaffold split remains in the original analysis and is a ligand holdout on its own.

Early-stopping validation. HistGBM early stopping uses an internal 10% of the training rows (validation_fraction 0.1, n_iter_no_change 15, max_iter 200, random_state 42). The exported validation split is written and is not that slice. The slice is not built from sequence clusters. Test kinases are outside the training rows, so their labels are outside the slice. Direct use of test labels for early stopping is closed. Dependence on a random training slice, and the mismatch with the exported validation clusters, stays open. With early stopping on, the 40% Spearman difference is -0.0024 (95% interval -0.0277 to 0.0236, includes 0; two-sided p = 0.8671). With early stopping off and the same 200-tree cap, it is -0.0110 (95% interval -0.0344 to 0.0138, includes 0; two-sided p = 0.3976). The early-stopping interval and the no-early-stopping interval both include 0. Ligand-only predictions on the 40% split have maximum within-drug range 0.0000.

## What was specified before the fits

Documented in DESIGN.md on 2026-10-02: Spearman and enrichment as ranking estimands, a scaffold split and a cold-protein split, a ligand-only control, and an ECFP plus gradient-boosting baseline. The binder cut pKd >= 7 is in configs/default.yaml and in artifacts/DECISIONS.md. DESIGN.md does not name amino-acid composition or dipeptide composition.

Documented with the results, in artifacts/DECISIONS.md: amino-acid composition as the first protein feature, then a v1.1 ablation. AAC+DPC is the mode that had the best cheap cold-protein Spearman in that ablation. scripts/uncertainty_bootstrap.py labels the AAC+DPC contrast pre-specified, and the same selection note says AAC+DPC was the best of the modes scored on that split. This revision treats the AAC+DPC choice as post hoc relative to DESIGN.md. The 0.05 Spearman-or-enrichment gate is in scripts/run_protein_ablation.py and in the decision log. It is not in DESIGN.md. The kinase-cluster bootstrap was added in scripts/uncertainty_bootstrap.py. It is not in the 2026-10-02 design lock. The published interval is not adjusted for the feature search that selected AAC+DPC.

Documented in configs/v1.7_reviewer.yaml before the fits of this revision: thresholds 0.40, 0.50, and 0.70, with 0.40 primary; identity bin edges 0.00, 0.30, 0.40, 0.50, 0.70, and 1.00; the censoring rule; four Tobit fits and a drop-censored refit on the original cold split and the 40% split; the additive baseline; within-kinase Spearman; and early stopping turned off on the 40% split. Those choices are pre-specified for this revision. They are not part of the 2026-10-02 lock. Seed 42 and 1000 bootstrap draws are in that file. The saved-score recompute of the published Spearman difference and its lower bound matches within 1e-9: True (absolute point difference 0.0, absolute lower-bound difference 0.0).

## Effect sizes

Intervals are kinase-cluster percentile intervals unless a row says the interval was not computed. The original row is the published file, in the saved score file's row order. The cold-protein row from this revision is a bootstrap after sorting by pair_id, so its interval can differ even when the point matches. Four-decimal values are formatted from artifacts/metrics_v1.7.json.

| Split | Pairs | Kinases | Spearman, protein features | Spearman, ligand-only | Delta Spearman | Delta enrichment at 1% |
|---|---|---|---|---|---|---|
| Original cold protein, published | 5168 | 76 | 0.5501 | 0.5002 | 0.0500 (95% interval 0.0271 to 0.0720, excludes 0; two-sided p = 0.0010) | 1.7930 (95% interval -0.2624 to 3.2468, includes 0; two-sided p = 0.0909) |
| Cold protein, this file's pair-id order | 5168 | 76 | 0.5501 | 0.5002 | 0.0500 (95% interval 0.0274 to 0.0719, excludes 0; two-sided p = 0.0010) | 2.3053 (95% interval -0.4827 to 3.1638, includes 0; two-sided p = 0.0919) |
| Sequence cluster 0.40 | 5236 | 77 | 0.4987 | 0.5011 | -0.0024 (95% interval -0.0277 to 0.0236, includes 0; two-sided p = 0.8671) | -1.7485 (95% interval -2.6863 to 1.2587, includes 0; two-sided p = 0.3287) |
| Sequence cluster 0.50 | 5168 | 76 | 0.4820 | 0.4692 | 0.0128 (95% interval -0.0147 to 0.0427, includes 0; two-sided p = 0.3756) | 0.2761 (95% interval -1.5829 to 2.2438, includes 0; two-sided p = 0.8472) |
| Sequence cluster 0.70 | 5168 | 76 | 0.5239 | 0.5149 | 0.0090 (95% interval -0.0202 to 0.0391, includes 0; two-sided p = 0.5584) | 0.5951 (95% interval -1.8182 to 2.5340, includes 0; two-sided p = 0.6613) |

Published enrichment at 1% is 10.2458 for the protein-feature scores and 8.4528 for ligand-only. After those same scores are sorted by pair_id, protein-feature enrichment is 10.2458 and ligand-only enrichment is 7.9405. Ligand-only scores are constant within a drug, so the identities of the rows inside the top 1% can change when tied rows change order. The published enrichment difference is 1.7930 (95% interval -0.2624 to 3.2468, includes 0; two-sided p = 0.0909). The pair-id-sorted difference is 2.3053 (95% interval -0.4827 to 3.1638, includes 0; two-sided p = 0.0919).

Uncensored rows, same splits:

| Split | Pairs | Kinases | Spearman, protein features | Spearman, ligand-only | Delta Spearman | Delta enrichment at 1% |
|---|---|---|---|---|---|---|
| Cold protein, this file's pair-id order | 1508 | 76 | 0.3960 | 0.3185 | 0.0775 (95% interval 0.0457 to 0.1127, excludes 0; two-sided p = 0.0010) | 0.7287 (95% interval 0.0000 to 1.9499, includes 0; two-sided p = 0.1918) |
| Sequence cluster 0.40 | 1398 | 77 | 0.3314 | 0.3352 | -0.0039 (95% interval -0.0404 to 0.0379, includes 0; two-sided p = 0.8492) | -0.5891 (95% interval -1.8057 to 0.7239, includes 0; two-sided p = 0.4276) |
| Sequence cluster 0.50 | 1494 | 76 | 0.3727 | 0.3288 | 0.0439 (95% interval 0.0069 to 0.0801, excludes 0; two-sided p = 0.0200) | 0.5533 (95% interval -0.6446 to 1.9700, includes 0; two-sided p = 0.4745) |
| Sequence cluster 0.70 | 1446 | 76 | 0.3580 | 0.3332 | 0.0249 (95% interval -0.0204 to 0.0694, includes 0; two-sided p = 0.2667) | 0.5772 (95% interval -0.6035 to 2.2026, includes 0; two-sided p = 0.4486) |

Mean within-kinase Spearman difference on the 40% split, full labels: -0.0058 (finite kinases 77) (95% interval -0.0170 to 0.0063, includes 0; two-sided p = 0.3526). On uncensored rows: 0.0315 (finite kinases 77) (95% interval -0.0090 to 0.0830, includes 0; two-sided p = 0.1978).

Contrasts whose interval excludes 0 are listed here with the estimate. The primary 40% pooled Spearman interval is the answer sentence. A row in this list does not replace that sentence.

- Cold protein, this file's pair-id order, pooled Spearman: 0.0500 (95% interval 0.0274 to 0.0719, excludes 0; two-sided p = 0.0010)
- Cold protein, this file's pair-id order, uncensored Spearman: 0.0775 (95% interval 0.0457 to 0.1127, excludes 0; two-sided p = 0.0010)
- Cold protein, this file's pair-id order, within-kinase Spearman: 0.0271 (finite kinases 76) (95% interval 0.0131 to 0.0411, excludes 0; two-sided p = 0.0020)
- Cold protein, this file's pair-id order, within-kinase Spearman on uncensored rows: 0.0480 (finite kinases 75) (95% interval 0.0225 to 0.0757, excludes 0; two-sided p = 0.0020)
- Cold protein, this file's pair-id order, full model minus additive: 0.0512 (95% interval 0.0249 to 0.0771, excludes 0; two-sided p = 0.0010)
- Cold protein, this file's pair-id order, Tobit, full labels: 0.0546 (95% interval 0.0355 to 0.0736, excludes 0; two-sided p = 0.0010)
- Cold protein, this file's pair-id order, Tobit, uncensored rows: 0.0774 (95% interval 0.0484 to 0.1090, excludes 0; two-sided p = 0.0010)
- Cold protein, this file's pair-id order, drop-censored training, uncensored test: 0.0599 (95% interval 0.0336 to 0.0869, excludes 0; two-sided p = 0.0010)
- Cold protein, this file's pair-id order, drop-censored training, full test labels: 0.0948 (95% interval 0.0685 to 0.1194, excludes 0; two-sided p = 0.0010)
- Cold protein, this file's pair-id order, shorter-sequence identity bin [0.50, 0.70) (14 kinases): 0.0759 (95% interval 0.0016 to 0.1527, excludes 0; two-sided p = 0.0490)
- Cold protein, this file's pair-id order, shorter-sequence identity bin [0.70, 1.00] (13 kinases): 0.1004 (95% interval 0.0526 to 0.1466, excludes 0; two-sided p = 0.0010)
- Cold protein, this file's pair-id order, alignment identity bin [0.70, 1.00] (20 kinases): 0.0898 (95% interval 0.0497 to 0.1275, excludes 0; two-sided p = 0.0010)
- Cold protein, this file's pair-id order: Spearman between max shorter-sequence identity to the training set and the within-kinase Spearman difference (protein features minus ligand-only) is 0.2734 on 76 kinases (95% interval 0.0557 to 0.4841, excludes 0; two-sided p = 0.0130).
- Cold protein, this file's pair-id order: Spearman between max alignment identity to the training set and the within-kinase Spearman difference (protein features minus ligand-only) is 0.2713 on 76 kinases (95% interval 0.0567 to 0.4523, excludes 0; two-sided p = 0.0080).
- Sequence cluster 0.40, drop-censored training, uncensored test: 0.0466 (95% interval 0.0071 to 0.0838, excludes 0; two-sided p = 0.0130)
- Sequence cluster 0.40, drop-censored training, full test labels: 0.0932 (95% interval 0.0660 to 0.1199, excludes 0; two-sided p = 0.0010)
- Sequence cluster 0.50, uncensored Spearman: 0.0439 (95% interval 0.0069 to 0.0801, excludes 0; two-sided p = 0.0200)
- Sequence cluster 0.50, within-kinase Spearman on uncensored rows: 0.0311 (finite kinases 76) (95% interval 0.0083 to 0.0524, excludes 0; two-sided p = 0.0070)
- Sequence cluster 0.50, full model minus additive: 0.0291 (95% interval 0.0051 to 0.0509, excludes 0; two-sided p = 0.0170)
- Sequence cluster 0.50, alignment identity bin [0.70, 1.00] (10 kinases): 0.0656 (95% interval 0.0229 to 0.0988, excludes 0; two-sided p = 0.0010)
- Sequence cluster 0.50: Spearman between max alignment identity to the training set and the within-kinase Spearman difference (protein features minus ligand-only) is 0.2861 on 76 kinases (95% interval 0.0477 to 0.4961, excludes 0; two-sided p = 0.0120).
- Sequence cluster 0.70, within-kinase Spearman on uncensored rows: 0.0587 (finite kinases 76) (95% interval 0.0182 to 0.1106, excludes 0; two-sided p = 0.0200)
- Sequence cluster 0.70: Spearman between max shorter-sequence identity to the training set and the within-kinase Spearman difference (protein features minus ligand-only) is 0.3898 on 76 kinases (95% interval 0.1639 to 0.5730, excludes 0; two-sided p = 0.0020).

## Alternative explanations that were checked

Ligand memorization. On the 40% split the ligand-only Spearman is 0.5011 and the protein-feature Spearman is 0.4987. The difference is -0.0024 (95% interval -0.0277 to 0.0236, includes 0; two-sided p = 0.8671). The ligand-only model is constant within a drug on this fit, up to the within-drug range above.

Homolog leakage. The original cold-protein split leaves it open. The 40% cluster split and the identity bins are the check. Bin rows are in the report section. For nident / shorter length, both coverages at least 0.8: Intervals that include 0: [0.00, 0.30), [0.30, 0.40). Bins with a point and no interval, because fewer than 8 kinases were present: [0.40, 0.50) point 0.1556 with 1 kinase; [0.50, 0.70) point 0.1003 with 1 kinase.

Censoring mass point. The cap covers fraction 0.7117 of pairs. The uncensored 40% difference is -0.0039 (95% interval -0.0404 to 0.0379, includes 0; two-sided p = 0.8492). The secondary Tobit model, four Gaussian left-censoring updates at pKd 5, differs from its ligand-only Tobit fit by 0.0110 (95% interval -0.0096 to 0.0341, includes 0; two-sided p = 0.3207) on the full 40% test and by 0.0011 (95% interval -0.0348 to 0.0405, includes 0; two-sided p = 0.9690) on uncensored test rows. A model fit only on uncensored training rows differs by 0.0466 (95% interval 0.0071 to 0.0838, excludes 0; two-sided p = 0.0130) on uncensored test rows and by 0.0932 (95% interval 0.0660 to 0.1199, excludes 0; two-sided p = 0.0010) on the full test labels. The drop-censored fit is secondary. Where its interval excludes 0 and the primary 40% interval includes 0, the drop-censored fit does not replace the primary interval.

Kinase main effect versus an interaction. On the complete DAVIS matrix the balanced two-way shares are ligand 0.2610, kinase 0.1158, and residual 0.6231. The residual share is interaction plus measurement error. One observation per cell cannot separate those two. On uncensored rows, least-squares R^2 is ligand 0.1562, kinase 0.1172, and ligand plus kinase 0.3067. The share outside that additive fit is 0.6933. On the original cold-protein scores, additive minus ligand-only is -0.0012 (95% interval -0.0235 to 0.0199, includes 0; two-sided p = 0.9211). Full model minus additive is 0.0512 (95% interval 0.0249 to 0.0771, excludes 0; two-sided p = 0.0010). Mean within-kinase protein-minus-ligand is 0.0271 (finite kinases 76) (95% interval 0.0131 to 0.0411, excludes 0; two-sided p = 0.0020). On the 40% split, additive minus ligand-only is 0.0118 (95% interval -0.0037 to 0.0280, includes 0; two-sided p = 0.1598). Full model minus additive is -0.0142 (95% interval -0.0323 to 0.0029, includes 0; two-sided p = 0.1399). Mean within-kinase full-minus-additive is -0.0057 (finite kinases 77) (95% interval -0.0168 to 0.0066, includes 0; two-sided p = 0.3546). An oracle additive score that uses each test kinase's own labels has Spearman 0.5625. That number uses test labels.

Post hoc feature selection. AAC+DPC was selected on the original cold-protein split after the v1.1 modes were scored. The sequence-cluster thresholds were not used to pick the feature. The intervals in this revision are not multiplicity-adjusted for that earlier search.

Early-stopping slice. The no-early-stopping fit on the 40% split is the check, reported above. KIBA was not refit. Chemprop was not refit. A composition-neighbor holdout was not run.

## What would change the answer

The answer sentence is the 40% pooled Spearman difference and whether its interval includes 0. That sentence would change if the same kinase-cluster bootstrap on this same split produced the other include-0 result. A different feature set, a different model, or a split that was not written in configs/v1.7_reviewer.yaml before the fits would be a different question. The within-kinase mean on this split has an interval that includes 0. The pooled sentence does not say that protein features change ranks inside a kinase. The uncensored-row interval includes 0. The pooled sentence uses the full label vector, including the Kd cap. It does not extend to the uncensored pairs. The enrichment-at-1% interval includes 0. The answer sentence is about Spearman. The enrichment difference is -1.7485 (95% interval -2.6863 to 1.2587, includes 0; two-sided p = 0.3287). On this split the pooled difference, the within-kinase difference, and the full-minus-additive difference all have intervals that include 0. The original cold-protein interval remains the published interval. A 40% interval that includes 0 means that published difference does not carry over to this sequence-cluster split as a difference separated from 0.

## Limits of this revision

DAVIS is one complete matrix of 68 ligands by 379 kinases. The additive scaler is fit on all training kinases before the ridge penalty is chosen. Test kinases are outside that fit. The penalty on the 40% training kinases is 1000.0000, with training-kinase R^2 0.4394. The Tobit fit assumes Gaussian residuals and is labeled secondary. Tree counts are in artifacts/metrics_v1.7.json. 7 of 16 recorded fits reached 200 trees. The original analysis, including its published interval, is unchanged.
<!-- V17_README_END -->

## Known limitations

1. Held-out kinases are split by kinase identity, not by sequence similarity, so close homologs can sit in both train and test.
2. DAVIS labels are censored at Kd = 10000 nM (pKd 5), and this is not yet modeled.
3. There is no additive ligand-plus-kinase baseline yet.
4. The top-1% enrichment gain is not settled, since its interval includes 0.

The Discussion in [RESEARCH_REPORT.md](RESEARCH_REPORT.md) has the detail. The cold-protein split groups rows by `target_id` (`src/vs/splits.py`). Labels are fit as exact `pKd = 9 - log10(Kd_nM)` values (`src/vs/data.py`).

The full write-up, with every number taken from the exported tables, is [RESEARCH_REPORT.md](RESEARCH_REPORT.md).

![Cold-protein and scaffold-split comparisons. Protein features raise rank correlation when scaffolds are held out. On held-out kinases the Spearman difference is smaller, and enrichment at 1% is not a settled gain.](figures/kinase-binding-figures.png)

## Repo map

- [README.md](README.md): this page, including how to check the saved scores.
- [RESEARCH_REPORT.md](RESEARCH_REPORT.md): the full write-up. Reported numbers there are the exported tables.
- [figures/kinase-binding-figures.png](figures/kinase-binding-figures.png): cold-protein and scaffold-split figure.
- [artifacts/metrics_uncertainty.md](artifacts/metrics_uncertainty.md): bootstrap intervals from the saved scores.
- [artifacts/metrics_v1.7.md](artifacts/metrics_v1.7.md): sequence-identity holdouts, censored labels, and the additive baseline.

## Data

- DAVIS pairs, the primary labels, downloaded from Therapeutics Data Commons via Harvard Dataverse file 5219748: https://dataverse.harvard.edu/api/access/datafile/5219748. Papers already cited in the report: Davis et al. 2011 (https://doi.org/10.1038/nbt.1990) and Huang et al. 2021, Therapeutics Data Commons.
- KIBA, the secondary table. The TDC Dataverse file returned HTTP 403, so the frozen table came from the DeepDTA public mirror: https://github.com/hkmztrk/DeepDTA/tree/master/data/kiba. Papers: Tang et al. 2014 (https://doi.org/10.1021/ci400709d) and Öztürk et al. 2018 (https://doi.org/10.1093/bioinformatics/bty593). The three KIBA files over 100 MB (`data/raw/kiba_pairs.csv`, `data/splits_kiba/cold_protein_split.csv`, and `data/splits_kiba/scaffold_split.csv`) are not in the repo. The smaller DeepDTA text files that are in the repo (`data/raw/kiba_deepdta/ligands_can.txt`, `data/raw/kiba_deepdta/proteins.txt`, and `data/raw/kiba_deepdta/kiba_binding_affinity_v2.txt`) came from that same mirror.
- The screening library is a Murcko subsample of MoleculeNet HIV plus the 68 DAVIS ligands. HIV activity labels were not used. Source file: https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/HIV.csv

## Data licenses

The MIT text in [LICENSE](LICENSE) covers this repository's code. It does not cover the data files.

- DAVIS, via Therapeutics Data Commons, Harvard Dataverse file 5219748. The TDC drug-target page lists the DAVIS dataset license as "Not Specified" and links [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/): https://tdcommons.ai/multi_pred_tasks/dti/#davis. Papers already cited in the report: Davis et al. 2011 and Huang et al. 2021.
- KIBA, from the DeepDTA public mirror, because the TDC file returned HTTP 403: https://github.com/hkmztrk/DeepDTA/tree/master/data/kiba. That repository has no LICENSE file. Papers already cited: Tang et al. 2014 and Öztürk et al. 2018. TDC's own KIBA entry, which this repo did not download, also lists the license as "Not Specified" and links CC BY 4.0: https://tdcommons.ai/multi_pred_tasks/dti/#kiba.
- MoleculeNet HIV, the CSV at https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/HIV.csv. Wu et al. 2018 (https://doi.org/10.1039/C7SC02664A), already cited in [docs/process/RESEARCH_PACK.md](docs/process/RESEARCH_PACK.md), describe this table as the NCI Drug Therapeutics Program AIDS Antiviral Screen. The MoleculeNet code repository is MIT (https://github.com/deepchem/moleculenet) and does not contain this CSV. This note does not assign the CSV a license of its own. HIV activity labels were not used.

## Original analysis

The reviewer revision is the opening section. The points in this section are the original splits and the original exported tables.

The working scorer is a histogram gradient booster on Morgan fingerprints, plus amino-acid composition and dipeptide composition. On the saved DAVIS splits, those protein features beat a ligand-only control when scaffolds are held out, and they add a smaller lift when whole proteins are held out (cold-protein Spearman 0.5501 vs 0.5002). A Chemprop graph model did not keep that cold-protein lift. The same features transferred to KIBA on Spearman, not on enrichment at 1%. A one-pair selectivity check (SRC vs LCK) was positive, and ligand-only could not rank it. A warm screen of LCK looked strong because that kinase was in training. A cold screen of SLK stayed above chance, and ligand-only won the broader ranks.

## Reproduce

```bash
uv venv --python 3.12 .venv && source .venv/bin/activate
uv pip install -e .
./reproduce.sh
```

Sequence-identity holdouts, censored-label metrics, and the additive baseline:

```bash
./reproduce_v1.7.sh
```

`reproduce_v1.7.sh` installs pandas 2.2.3, numpy 1.26.4, scikit-learn 1.5.2, rdkit 2024.3.5, scipy 1.13.1, joblib 1.4.2, and pyyaml 6.0.2 if the virtual environment does not already have them. It then runs `scripts/run_reviewer_v1.7.py`, which downloads MMseqs2 release 17-b804f (pinned sha256 in `configs/v1.7_reviewer.yaml`) and the DAVIS table from Harvard Dataverse file 5219748 when those files are absent, checks the DAVIS sha256 against `data/MANIFEST.json`, and writes `artifacts/metrics_v1.7.json` plus `artifacts/metrics_v1.7.md`. The script sets `OMP_NUM_THREADS`, `MKL_NUM_THREADS`, and `OPENBLAS_NUM_THREADS` to 4. It does not rewrite `artifacts/metrics_v1.1.json`, `artifacts/metrics_uncertainty.json`, or `artifacts/scores_histgbm_cold_protein.csv`.

`reproduce.sh` checks the v1 metrics within 1e-4. Later shell entry points that exist are `reproduce_v1.1.sh`, `reproduce_v1.2.sh`, and `reproduce_v1.4.sh`. There is no `reproduce_v1.3.sh`, `reproduce_v1.5.sh`, or `reproduce_v1.6.sh`. The matching Python scripts in the tree are `scripts/run_kiba_v1.3.py`, `scripts/run_calibration.py`, `scripts/run_screen_v1.5.py`, and `scripts/run_screen_v1.6.py`. Raw DAVIS tables, the KIBA pair table, the ESM cache, and the KIBA split CSVs are on disk in this working copy and are listed in `.gitignore`, so a public commit should not include them. Manifests under `data/` record the source files and hashes. Rebuilding the KIBA splits needs the frozen pair table, which is also gitignored because it is over GitHub's 100 MB file limit.

## What is not claimed

No new binder. No prospective screen. No external concordance check. See the Discussion in [RESEARCH_REPORT.md](RESEARCH_REPORT.md).

## Verify

From the project root, with the existing `.venv`:

```bash
.venv/bin/python scripts/uncertainty_bootstrap.py
```

That command does not refit any model. It reads saved scores only and rewrites `artifacts/metrics_uncertainty.json` and `artifacts/metrics_uncertainty.md`. It also refreshes section 6 of `docs/process/STAT_REVIEW.md`.

The cold-protein histogram comparison is one of those saved-score checks. `artifacts/scores_histgbm_cold_protein.csv` stores per-pair ligand-only and amino-acid-plus-dipeptide scores for the frozen DAVIS cold-protein test (5,168 pairs, 76 kinases). [RESEARCH_REPORT.md](RESEARCH_REPORT.md) and [artifacts/metrics_uncertainty.md](artifacts/metrics_uncertainty.md) report the Spearman difference 0.0500 with a 95% kinase-cluster bootstrap interval of 0.0271 to 0.0720 (two-sided p = 0.001 in the report; the metrics table prints 0.0010). That interval excludes 0. The same sources report the enrichment-at-1% difference 1.7930 with a 95% interval of -0.2624 to 3.2468 (two-sided p = 0.091 in the report; the metrics table prints 0.0909). That interval includes 0. The script resamples held-out kinases from the CSV. It does not impute the interval from the aggregates in `artifacts/metrics_v1.1.json`.

The interval can also be recomputed for the graph model, because `artifacts/models_v1.2/preds_cold_protein_protein.npz`, `preds_cold_protein_ligand_only.npz`, `preds_scaffold_protein.npz`, and `preds_scaffold_ligand_only.npz` store per-pair scores. It can also be recomputed for the SLK recovery comparison, because `artifacts/screen_ranked_cold.csv` stores both scores and the known-binder flag.

It cannot be recomputed for the histogram amino-acid-plus-dipeptide versus ligand-only comparison on the DAVIS scaffold split. Those per-pair scores were not saved. The same is true for KIBA: `artifacts/metrics_kiba.json` has aggregates, not per-pair scores. The script skips those comparisons and does not impute an interval. The KIBA cold-protein difference of 0.0609 stays a point estimate.

`./reproduce_v1.2.sh` retrains the graph model and can replace the npz files above. Run the read-only checks first if those saved scores should stay as they are.

A reviewer recomputes the ranking math in this order:

1. `src/vs/metrics_lib.py` is the Spearman rho and enrichment definition used by the training scripts.
2. `./reproduce.sh` runs `scripts/run_train_eval.py` and `scripts/check_metrics_tolerance.py`, refits the first DAVIS models, and writes `artifacts/metrics.json`.
3. `./reproduce_v1.1.sh` runs `scripts/run_protein_ablation.py` and refits the histogram protein-feature comparison. It writes `artifacts/metrics_v1.1.json`. Cold-protein per-pair scores from the saved fit are already in `artifacts/scores_histgbm_cold_protein.csv`. `scripts/score_histgbm_cold_protein.py` does not refit. It loads the joblibs in `artifacts/models_v1.1/`, scores the frozen cold-protein test, and rewrites that CSV only if Spearman rho and enrichment at 1% match the exported point estimates within 1e-4. It does not rewrite `artifacts/metrics_v1.1.json`. Those joblib files are committed. `.gitignore` does not exclude `artifacts/models_v1.1/`. Regenerating the CSV also needs `data/raw/davis_pairs.csv`, which `.gitignore` does exclude. The bootstrap in step 7 reads the saved CSV and does not need that raw table.

```bash
.venv/bin/python scripts/score_histgbm_cold_protein.py
```

4. Graph-model Spearman rho and enrichment, without refitting, from the four saved npz files. Run this from the project root after `pip install -e .` so `vs` imports. `./reproduce_v1.2.sh` (`scripts/run_chemprop_v1.2.py`) is the refit path for the same metrics and is not required to check the saved scores.

```bash
.venv/bin/python - << 'PY'
import numpy as np
from vs.metrics_lib import compute_metrics
files = [
    "artifacts/models_v1.2/preds_cold_protein_protein.npz",
    "artifacts/models_v1.2/preds_cold_protein_ligand_only.npz",
    "artifacts/models_v1.2/preds_scaffold_protein.npz",
    "artifacts/models_v1.2/preds_scaffold_ligand_only.npz",
]
for path in files:
    z = np.load(path)
    m = compute_metrics(z["y_true"], z["y_pred"], 7.0, [0.01, 0.05])
    print(path, {k: round(m[k], 4) for k in ("spearman_rho", "ef_at_1pct", "ef_at_5pct", "auroc")})
PY
```

5. `scripts/run_kiba_v1.3.py` refits KIBA and is the only way to recompute that Spearman rho and enrichment. Per-pair scores are not saved, so there is no interval.
6. `scripts/run_screen_v1.6.py` refits the SLK screen. The saved table is `artifacts/screen_ranked_cold.csv`. Recovery on that table, including the new intervals, is step 7. The LCK table `artifacts/screen_ranked.csv` has scores but not a known-binder flag, and no ligand-only recovery table was exported, so that warm-screen comparison is not bootstrapped.
7. `.venv/bin/python scripts/uncertainty_bootstrap.py` writes the bootstrap intervals from the saved histogram cold-protein scores, the saved graph scores, and the saved SLK table.
8. `./reproduce_v1.7.sh` refits the sequence-identity models and rewrites `artifacts/metrics_v1.7.json`. A run that only rebuilds the prose from that json is `.venv/bin/python scripts/run_reviewer_v1.7.py --stage report`.
