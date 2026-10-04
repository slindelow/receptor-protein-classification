# Ranking kinase ligands from public affinity labels

Sofia Lindelow
2026-10-03 (America/Toronto)

Public write-up of the local experiments. Numbers are copied from exported artifacts under `artifacts/` and from `data/MANIFEST.json`, `data/MANIFEST_kiba.json`, `data/screen_library/MANIFEST_screen.json`, `DESIGN.md`, `LAB_NOTEBOOK.md`, `artifacts/DECISIONS.md`, and `configs/`. Models were not rerun for this report. If a figure was not in an export, it is not here.

## 1. Abstract

I asked whether a small model can rank protein-ligand affinity well enough to sort a public compound library for one kinase, and whether that rank uses the protein or only the ligand. Labels are experimental DAVIS Kd values converted to pKd, with a binder cut at pKd >= 7. The working scorer is a histogram gradient booster on Morgan fingerprints, with amino-acid composition plus dipeptide composition as the protein features that won the cold-protein check. On the saved DAVIS splits, those protein features beat a ligand-only control when scaffolds are held out, and they add a smaller lift when whole proteins are held out (cold-protein Spearman 0.5501 vs 0.5002; EF@1% 10.2458 vs 8.4528). A Chemprop graph encoder did not keep that cold-protein lift. The same composition features did transfer to a KIBA secondary freeze on Spearman, not on EF@1%. A one-pair selectivity sketch on SRC versus LCK was positive and ligand-only was a null. Screening an 8068-molecule library for LCK recovered known binders far above chance, but LCK was not held out of training. Holding out SLK reversed the story: recovery stayed above chance, and the ligand-only model won the broader ranks. These are predicted scores on public labels. They are not wet-lab hits.

## 2. Question and estimand

The product question, locked in `DESIGN.md` on 2026-10-02 (ET), is: given a target sequence, rank a fixed public library. The scientific object underneath is a score for one (protein, molecule) pair.

What is predicted is continuous affinity, not a class. On DAVIS the label is pKd = 9 - log10(Kd in nM). The model is `HistGradientBoostingRegressor` (or, in v1.2 only, a Chemprop regressor). Binary labels exist only to score enrichment and AUROC. A binder is pKd >= 7, which is Kd <= 100 nM. That cut is locked in `configs/default.yaml` and was not changed in later configs.

Ranking is the primary estimand. On held-out experimental pairs I report Spearman rho of predicted score versus pKd, and enrichment factor at 1% and 5% of the ranked test list (EF@1%, EF@5%). AUROC at the binder cut is secondary discrimination. Calibration asks a different question: after mapping the score to P(binder), does the probability match the observed binder rate (ECE, 10 bins, and Brier). Calibration does not change the ranking metrics.

The control is a ligand-only model with protein features removed, on the same rows. If it matches the full model, the rank is ligand-driven.

Two holdouts define the claim, both aiming at about 20% test and about 10% validation, seed 42:

- Scaffold: Murcko scaffolds held out.
- Cold-protein: entire `target_id`s held out.

Success is not a high AUROC on a random split. I do not claim a new binder unless an experiment measured it. NOVA or other external submission lists were optional concordance checks in the design. They were not run.

Later experiments keep that estimand and add three narrower ones:

- v1.3: same ranking metrics on KIBA score, not pKd, with binder defined as KIBA score >= 12.1.
- v1.4: Spearman of predicted delta versus true delta pKd on one kinase pair, plus enrichment of ligands with |delta pKd| >= 1.
- v1.5 and v1.6: recovery of known DAVIS binders inside a ranked library (fraction in the top 1% and in the top 500, and mean rank), versus a ligand-only rank and versus a random rank. v1.6 is the cold-target version of that recovery.

## 3. Data

### DAVIS (primary)

Freeze id: `sofia-vs-v1-2026-10`. Source: Therapeutics Data Commons multi_pred DTI DAVIS, Harvard Dataverse file `5219748`, downloaded 2026-10-02 (America/Toronto). Citation recorded in the manifest: Huang et al., Therapeutics Data Commons, NeurIPS 2021; Davis et al., Nat Biotechnol 2011. Raw file `data/raw/davis_tdc.tab`. Sanitized pairs `data/raw/davis_pairs.csv`.

From `data/MANIFEST.json`:

| Item | Value |
|---|---|
| Pairs | 25772 |
| Unique ligands | 68 |
| Unique targets | 379 |
| pKd | 9 - log10(Kd_nM); Kd in the TDC file is nM |
| Raw SHA256 | `6d4c6809dcb7c5da2b91a32d594d6935b75484940bde4d18055eb5e1059262f4` |
| Pairs SHA256 | `6b6281de1815c6a4131299091e64d91284b9fb151022a00c0f4bb9e2e687098e` |

PyTDC 0.4.1 was not the runtime reader. It depends on `rdkit-pypi`, which had no cp312 wheel. DAVIS was fetched from the same Dataverse URL TDC uses. `rdkit==2024.3.5` is pinned.

DAVIS has 68 ligands, each measured across the kinase panel. That is why a ligand-only model can look strong on a cold protein: the chemotypes are not new.

Split id lists were written once to `data/splits/` and later experiments loaded them. They were not regenerated. Exact train and validation counts for the v1 scaffold and cold-protein splits were not printed in `artifacts/metrics.md`. What that file does export is test size: scaffold n_test = 5306 (302 binders); cold-protein n_test = 5168 (388 binders). Later screen fits state their own train sizes (below).

### KIBA (secondary, v1.3 only)

Freeze id: `sofia-vs-kiba-v1.3`. The TDC Dataverse file returned HTTP 403 from this box. Pairs were assembled from the DeepDTA public mirror (`data/raw/kiba_deepdta/` to `data/raw/kiba_pairs.csv`). Documented in `data/MANIFEST_kiba.json`. Citation recorded there: Tang et al. JCIM 2014; Ozturk et al. DeepDTA Bioinformatics 2018.

| Item | Value |
|---|---|
| Pairs | 118254 |
| Unique ligands | 2068 |
| Unique targets | 229 |
| Label | KIBA score, stored in a `pkd` column so the same code path could run. Not pKd. |
| Binder cut | KIBA score >= 12.1 |

DAVIS split files were not touched. KIBA splits are under `data/splits_kiba/`, same fraction and seed protocol as `configs/default.yaml`.

### Screen library (v1.5, reused in v1.6)

Not a training label set. From `data/screen_library/MANIFEST_screen.json` and `artifacts/metrics_screen.md`:

| Item | Value |
|---|---|
| Source | MoleculeNet HIV, DeepChem S3 `HIV.csv` |
| URL | `https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/HIV.csv` |
| Raw SHA256 | `9ffa7fe57dc86c342627ee1d5255e937e2ab812393c73c4d16c697022f6e1d22` |
| Built | 2026-10-03T08:17:53-04:00 |
| HIV molecules kept | 8000 (diverse Murcko subsample, seed 42) |
| DAVIS ligands spiked in | 68 |
| Library n | 8068 |
| Unique Murcko scaffolds | 7582 |
| Library SHA256 | `5721c043c87b744e4f288d31319e7743bf72e7a064de658708b0a16c61df0604` |
| Path | `data/screen_library/library_v1.5.csv` |

HIV activity labels were not used. v1.6 reused this file and did not redownload it.

An earlier demo library (`data/demo/library.csv`, n=400) pads DAVIS ligands with public Delaney/ESOL SMILES. ESOL is not training labels. The demo library is the CLI toy. The 8068 library is the v1.5/v1.6 product.

## 4. Methods

### Features

Ligand features in every HistGBM run: Morgan/ECFP, radius 2, 2048 bits.

Protein features, concatenated to the fingerprint:

| Mode | What it is | Dim exported in v1.1 (`feat_dim` includes ECFP) |
|---|---|---|
| `aac` | Amino-acid composition | 2068 (20 + 2048) |
| `dpc` | Dipeptide composition | 2448 |
| `aac_dpc` | Both | 2468 |
| `esm2` | Frozen ESM-2 `esm2_t6_8M_UR50D`, mean pool, dim 320, sequences truncated at 1022 residues | 2368 |
| ligand-only | ECFP only | 2048 |

MVP (v1) used `aac` only, because composition can be computed for an unseen sequence, unlike a target-id one-hot. v1.1 compared the four modes above. From v1.4 onward the deployed protein features are `aac_dpc`.

v1.2 replaces the ECFP ligand encoder with Chemprop bond message passing (D-MPNN). Protein features, when used, are `aac_dpc` passed as Chemprop `x_d`. ESM-2 was not the Chemprop extra. `aac_dpc` was chosen because it was the best cheap cold-protein mode in v1.1.

### Models

HistGBM hyperparameters, locked in `configs/default.yaml` and copied into later configs: `max_iter` 200, `learning_rate` 0.08, `max_depth` 8, `min_samples_leaf` 20, `l2_regularization` 0.1, early stopping on, `validation_fraction` 0.1, `n_iter_no_change` 15, `random_state` 42.

Chemprop (`configs/v1.2_chemprop.yaml`): chemprop 2.2.1, `max_epochs` 15, batch 64, hidden 300, depth 3, dropout 0, warmup 2 epochs, init/max/final lr 0.0001 / 0.001 / 0.0001, early stopping patience 5, CPU, seed 42. Checkpoint reload used `weights_only=False`. Wall clock for the four trainings was about 18 minutes (decision log, 2026-10-02).

Pins recorded for the MVP stack: Python 3.12, pandas 2.2.3, numpy 1.26.4, scikit-learn 1.5.2, rdkit 2024.3.5. v1.1 added a CPU torch wheel and `fair-esm==2.0.0`. The decision log names the torch build as 2.14.1. I am not restating a version that was not written there.

### Splits and what each experiment trained on

| Experiment | Train rows | Held out of training |
|---|---|---|
| v1, v1.1, v1.2 | Scaffold split and, separately, cold-protein split. IDs from `data/splits/`. | Test scaffolds, or test proteins. Val is inside the split file. |
| Calibration | No new fit of the ranker. Isotonic map only. | HistGBM: isotonic fit on val, scored on test. Chemprop: val scores were not saved; isotonic fit on a 30% test holdout. |
| v1.3 | New scaffold and cold-protein splits on KIBA only. | DAVIS splits unchanged. |
| v1.4 | All DAVIS pairs except SRC and LCK. n_train = 25636 (decision log). | Both kinases together (cold-pair). |
| v1.5 | Scaffold train+val. n_train = 20466. | Scaffold test not used for this fit. LCK itself was not held out. |
| v1.6 | All DAVIS pairs except the 68 SLK pairs. n_train = 25704. | SLK only. The same 68 ligands still have labels on other kinases. |

Screening weights are a separate fit from the test-set models. v1 `artifacts/models/hgb_screen.joblib` is the MVP screen model (AAC, scaffold train+val). v1.5 and v1.6 have their own joblibs. I did not reuse v1.4 weights for the screens.

### Metrics

Pair ranking: Spearman rho, EF@1%, EF@5%, AUROC. Binder cut as above.

Calibration: ECE with 10 bins, Brier, before and after isotonic regression. AUROC is reported again and does not move, which is the check that calibration did not reorder scores.

Selectivity: delta = score(SRC) - score(LCK), same definition on true pKd. Selective ligand: |delta pKd| >= 1. EF at 5%, 10%, and 20% of the 68 ligands, ranked by |delta pred|.

Library recovery: fraction of known binders (pKd >= 7 on that kinase) inside the top 1% and inside the top 500, plus mean rank. Top 1% uses k = ceil(0.01 * 8068) = 81. Random expected counts and random mean rank are in the export. They are not a second model.

### Verifiability

Frozen inputs, SHA256, and download date are in the manifests. Seeds and thresholds are in YAML. Split ids are files, not a paragraph. Metrics markdown is written by the runner scripts. The instruction in those files is not to hand-edit them.

Reproduce scripts and tolerances that exist:

| Script | Tolerance |
|---|---|
| `reproduce.sh` | 1e-4 absolute on Spearman, AUROC, EF. Notebook: tolerance OK on 2026-10-02. |
| `reproduce_v1.1.sh` | 1e-4 vs `artifacts/.metrics_v1.1_reference.json`. Notebook: tolerance OK, 40 metrics, 2026-10-02 16:41 EDT. |
| `reproduce_v1.2.sh` | Tries 1e-4, then a documented fallback of 0.05 because Chemprop on CPU is not bit-stable. |
| `reproduce_v1.4.sh` | 1e-4 vs `artifacts/.metrics_v1.4_reference.json`. The decision log says the full retrain was not re-run in that session; the first export is the reference. |

No reproduce script for calibration, KIBA, v1.5, or v1.6 was listed in the project root at the time of this write-up. Those results rest on the exported markdown and the decision log.

## 5. Results

Preferred scorer after v1.2, and the one used for v1.4 through v1.6: HistGBM + ECFP + `aac_dpc`.

Tables below are the export tables. Rounded one-pagers in `artifacts/RESULTS_SUMMARY.md` match them at the precision that page chose. Where they differ in digits, this report keeps the metrics file.

### 5.1 v1: ECFP + HistGBM + AAC

Freeze `sofia-vs-v1-2026-10`. Generated 2026-10-02T14:52:41-04:00. Question: does AAC help versus ligand-only, on scaffold and on cold-protein?

Source: `artifacts/metrics.md`.

| Split | Model | Spearman ρ | EF@1% | EF@5% | AUROC | n_test | n_binders |
|---|---|---:|---:|---:|---:|---:|---:|
| scaffold | ECFP+AAC | 0.2124 | 6.5072 | 3.3025 | 0.6276 | 5306 | 302 |
| scaffold | ligand-only (ECFP) | -0.0107 | 0.0000 | 0.1982 | 0.4413 | 5306 | 302 |
| cold_protein | ECFP+AAC | 0.4984 | 8.7090 | 6.1712 | 0.8352 | 5168 | 388 |
| cold_protein | ligand-only (ECFP) | 0.5002 | 8.4528 | 6.5312 | 0.8280 | 5168 | 388 |

Scaffold: protein features are the whole signal. Ligand-only Spearman is about zero and EF@1% is zero. Cold-protein: AAC does not help. Rho is 0.4984 versus 0.5002. EF@1% is a bit higher (8.7090 vs 8.4528) and EF@5% is a bit lower (6.1712 vs 6.5312). AUROC is close (0.8352 vs 0.8280). Cold-protein ranking on this freeze is largely ligand-driven. That is the result the ablation was there to catch. I do not read cold-protein rho of 0.50 as evidence that AAC generalized.

### 5.2 v1.1: stronger protein features

Experiment freeze `sofia-vs-v1.1-protein`. Same data freeze. Generated 2026-10-02T16:35:33-04:00. Same splits, same HistGBM, same binder cut. ESM-2 ran on CPU.

Source: `artifacts/metrics_v1.1.md`.

| Split | Protein mode | Model | Spearman ρ | EF@1% | EF@5% | AUROC | n_test | n_binders | feat_dim |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| scaffold | n/a | ligand-only (ECFP) | -0.0107 | 0.0000 | 0.1982 | 0.4413 | 5306 | 302 | 2048 |
| scaffold | aac | ECFP+AAC | 0.2124 | 6.5072 | 3.3025 | 0.6276 | 5306 | 302 | 2068 |
| scaffold | dpc | ECFP+DPC | 0.2219 | 4.8804 | 3.5007 | 0.6341 | 5306 | 302 | 2448 |
| scaffold | aac_dpc | ECFP+AAC+DPC | 0.2397 | 5.2058 | 3.6328 | 0.6437 | 5306 | 302 | 2468 |
| scaffold | esm2 | ECFP+ESM2 | 0.2566 | 5.8565 | 4.0952 | 0.6702 | 5306 | 302 | 2368 |
| cold_protein | n/a | ligand-only (ECFP) | 0.5002 | 8.4528 | 6.5312 | 0.8280 | 5168 | 388 | 2048 |
| cold_protein | aac | ECFP+AAC | 0.4984 | 8.7090 | 6.1712 | 0.8352 | 5168 | 388 | 2068 |
| cold_protein | dpc | ECFP+DPC | 0.5409 | 9.9897 | 7.1484 | 0.8652 | 5168 | 388 | 2448 |
| cold_protein | aac_dpc | ECFP+AAC+DPC | 0.5501 | 10.2458 | 6.9941 | 0.8689 | 5168 | 388 | 2468 |
| cold_protein | esm2 | ECFP+ESM2 | 0.5476 | 10.2458 | 7.4055 | 0.8783 | 5168 | 388 | 2368 |

AAC rows match v1 exactly. Splits were not rebuilt.

Primary gate as exported: best cold-protein mode is `aac_dpc`, rho 0.5501 vs ligand-only 0.5002 (delta 0.0500). EF@1% 10.2458 vs 8.4528. EF@5% 6.9941 vs 6.5312. Pass if delta Spearman >= 0.05 or a clear EF lift: YES. The rho gate itself is recorded False. The EF gate is True. The unrounded delta that failed the strict float check was not given as a separate number.

ESM-2 is close on cold-protein (rho 0.5476, same EF@1% 10.2458, higher EF@5% 7.4055, higher AUROC 0.8783). DPC alone already carries most of the cold-protein lift. AAC alone still does not.

On scaffold, every protein mode beats ligand-only. ESM-2 has the highest scaffold rho (0.2566) and AUROC (0.6702). `aac_dpc` is the mode I kept because it was the best cheap cold-protein Spearman, tied with ESM-2 on EF@1%, and it does not need an embedding model at train time.

### 5.3 v1.2: Chemprop D-MPNN

Experiment freeze `sofia-vs-v1.2-chemprop`. Generated 2026-10-02T17:06:25-04:00. Same DAVIS freeze and the same split CSVs. HistGBM was not retrained. The HistGBM rows below are cited from v1.1.

Source: `artifacts/metrics_v1.2.md`.

| Split | Model | Spearman rho | EF@1% | EF@5% | AUROC | n_test | n_binders |
|---|---|---:|---:|---:|---:|---:|---:|
| scaffold | Chemprop ligand-only | 0.2621 | 0.3254 | 0.5284 | 0.6645 | 5306 | 302 |
| scaffold | Chemprop+aac_dpc | 0.3176 | 7.4833 | 4.5575 | 0.7373 | 5306 | 302 |
| cold_protein | Chemprop ligand-only | 0.4997 | 8.7090 | 6.6341 | 0.8279 | 5168 | 388 |
| cold_protein | Chemprop+aac_dpc | 0.4895 | 8.4528 | 5.9655 | 0.8541 | 5168 | 388 |

Cited HistGBM baselines:

| Split | Model | Spearman rho | EF@1% | EF@5% | AUROC |
|---|---|---:|---:|---:|---:|
| scaffold | HistGBM ligand-only (ECFP) | -0.0107 | 0.0000 | 0.1982 | 0.4413 |
| scaffold | HistGBM+AAC+DPC | 0.2397 | 5.2058 | 3.6328 | 0.6437 |
| cold_protein | HistGBM ligand-only (ECFP) | 0.5002 | 8.4528 | 6.5312 | 0.8280 |
| cold_protein | HistGBM+AAC+DPC | 0.5501 | 10.2458 | 6.9941 | 0.8689 |

Primary test: Chemprop+`aac_dpc` versus Chemprop ligand-only on cold-protein. Rho 0.4895 vs 0.4997 (delta -0.0102). EF@1% 8.4528 vs 8.7090 (delta -0.2561). EF@5% 5.9655 vs 6.6341 (delta -0.6686). Primary pass: NO.

Against HistGBM+`aac_dpc` on cold-protein: Chemprop+protein rho 0.4895 vs 0.5501 (delta -0.0606).

Chemprop ligand-only on cold-protein matches HistGBM ligand-only (0.4997 vs 0.5002). The graph encoder is not worse than ECFP on that split. It also does not use the protein extra. AUROC for Chemprop+`aac_dpc` on cold-protein is 0.8541, above Chemprop ligand-only 0.8279, while Spearman and EF go down. Discrimination and ranking disagreed. I kept the ranking metrics as primary, so this is a failed protein extra, not a win.

Scaffold is a different picture. Chemprop ligand-only already reaches rho 0.2621, where HistGBM ligand-only was -0.0107, but its EF@1% is only 0.3254. Adding `aac_dpc` raises scaffold EF@1% to 7.4833 and rho to 0.3176. Protein features still matter when the held-out unit is the scaffold.

### 5.4 Calibration

Not a new ranker. Export `artifacts/metrics_calibration.md`, generated 2026-10-02T17:25:51-04:00. ECE uses 10 bins. HistGBM isotonic regression is fit on the validation split. Chemprop validation scores were not saved, so that isotonic fit uses a 30% holdout of the test set. The Chemprop ECE after calibration is not an unbiased test number.

| Split | Model | Protocol | ECE unc | ECE iso | Brier unc | Brier iso | AUROC | n_eval |
|---|---|---|---:|---:|---:|---:|---:|---:|
| scaffold | hgb_ligand_only | isotonic_fit_on_val | 0.1405 | 0.0750 | 0.0934 | 0.0602 | 0.4413 | 5306 |
| scaffold | hgb_aac_dpc | isotonic_fit_on_val | 0.0780 | 0.0482 | 0.0638 | 0.0572 | 0.6437 | 5306 |
| scaffold | hgb_esm2 | isotonic_fit_on_val | 0.0899 | 0.0528 | 0.0644 | 0.0567 | 0.6702 | 5306 |
| cold_protein | hgb_ligand_only | isotonic_fit_on_val | 0.1036 | 0.0097 | 0.0744 | 0.0573 | 0.8280 | 5168 |
| cold_protein | hgb_aac_dpc | isotonic_fit_on_val | 0.0849 | 0.0135 | 0.0613 | 0.0540 | 0.8689 | 5168 |
| cold_protein | hgb_esm2 | isotonic_fit_on_val | 0.0939 | 0.0147 | 0.0628 | 0.0538 | 0.8783 | 5168 |
| scaffold | chemprop_ligand_only | isotonic_on_30pct_test_holdout | 0.4108 | 0.0067 | 0.2685 | 0.0536 | 0.6645 | 3715 |
| scaffold | chemprop_protein_aac_dpc | isotonic_on_30pct_test_holdout | 0.2576 | 0.0071 | 0.1258 | 0.0509 | 0.7373 | 3715 |
| cold_protein | chemprop_ligand_only | isotonic_on_30pct_test_holdout | 0.1116 | 0.0049 | 0.0769 | 0.0573 | 0.8279 | 3618 |
| cold_protein | chemprop_protein_aac_dpc | isotonic_on_30pct_test_holdout | 0.1953 | 0.0072 | 0.0971 | 0.0565 | 0.8541 | 3618 |

On cold-protein HistGBM, isotonic regression cuts ECE from about 0.08 to 0.10 down to about 0.01 to 0.015. The exact pairs are in the table (ligand-only 0.1036 to 0.0097; `aac_dpc` 0.0849 to 0.0135; ESM-2 0.0939 to 0.0147). Scaffold HistGBM stays worse after calibration (0.0482 to 0.0750). Brier moves with ECE. AUROC is unchanged, so the map did not reorder pairs. I would use the isotonic scores only if a later step needs a probability. The screen ranks raw predicted pKd.

### 5.5 v1.3: KIBA transfer

Experiment freeze `sofia-vs-kiba-v1.3`. Generated 2026-10-02T17:26:49-04:00. HistGBM, ECFP, with and without `aac_dpc`, MVP hyperparameters. Label is KIBA score. Binder cut is 12.1, not pKd 7.

Source: `artifacts/metrics_kiba.md`.

| Split | Model | Spearman rho | EF@1% | EF@5% | AUROC | n_test | n_binders |
|---|---|---:|---:|---:|---:|---:|---:|
| scaffold | HistGBM ligand-only | 0.3394 | 1.6314 | 2.2196 | 0.6671 | 23663 | 4529 |
| scaffold | HistGBM+AAC+DPC | 0.5584 | 4.9602 | 3.8171 | 0.7704 | 23663 | 4529 |
| cold_protein | HistGBM ligand-only | 0.4806 | 3.9041 | 3.3085 | 0.7485 | 24280 | 5912 |
| cold_protein | HistGBM+AAC+DPC | 0.5415 | 3.8872 | 3.4980 | 0.7688 | 24280 | 5912 |

Cold-protein delta Spearman = 0.0609. Pass (delta >= 0.05 or EF lift): YES. The pass is the Spearman gap. EF@1% is flat to slightly down (3.8872 vs 3.9041). EF@5% and AUROC move up (3.4980 vs 3.3085; 0.7688 vs 0.7485).

Scaffold on KIBA is a large protein lift: rho 0.3394 to 0.5584, EF@1% 1.6314 to 4.9602. Unlike DAVIS, ligand-only scaffold rho is not zero. KIBA has 2068 ligands, not 68, so a scaffold holdout still leaves other chemotypes to learn from. I do not treat KIBA rho as a pKd result.

### 5.6 v1.4: selectivity on SRC / LCK

Experiment freeze `sofia-vs-v1.4-selectivity`. Timestamp 2026-10-02 17:43 EDT. Data freeze unchanged.

Pair choice, recorded before the fit: SRC and LCK, Src family. Best-window sequence identity 0.4008. AAC+DPC cosine 0.9621. All 68 DAVIS ligands are measured on both. 12 ligands have |delta pKd| >= 1. The decision log says closer pairs (ERK1/ERK2, CDK2/CDK3) had almost no selectivity signal, and ABL1/ABL2 looked similar by composition but not by an N-terminal identity window. Those alternate pairs were not scored in the metrics export, so I am not adding numbers for them.

Protocol: cold-pair holdout. Train excludes both SRC and LCK (n_train 25636). Scorer: HistGBM + ECFP + `aac_dpc`. Delta = score(SRC) - score(LCK).

Source: `artifacts/metrics_selectivity.md`.

| Model | Spearman(Δpred, Δtrue) | mean|Δpred| | EF@5% | EF@10% | EF@20% |
|---|---:|---:|---:|---:|---:|
| HistGBM+aac_dpc | 0.3694 | 0.1938 | 2.83 | 3.24 | 2.43 |
| HistGBM ligand-only | n/a (constant Δ) | 0.0000 | 0.00 | 0.00 | 0.00 |

EF is enrichment of the 12 selective ligands among 68, ranked by |delta pred|. All 12 prefer LCK (delta pKd <= -1). None prefer SRC by at least 1. LCK-directed EF equals the absolute EF. SRC-directed EF is undefined because there are zero positives.

Ligand-only delta is exactly 0. Spearman is undefined. EF is 0. A model with no protein features cannot rank a selectivity gap. That control worked.

Sketch gate (rho > 0.1 and EF@10% or EF@20% > 1.2): YES. This is one pair, n=68, and every selective ligand points the same way. It is not a selectivity product.

### 5.7 v1.5: screen a library for LCK

Experiment freeze `sofia-vs-v1.5-screen`. Generated 2026-10-03T08:26:08-04:00. Question: can HistGBM+`aac_dpc` rank about 8k public molecules for one DAVIS kinase?

Weights: `artifacts/models_v1.5/hgb_screen_aac_dpc.joblib`, fit on scaffold train+val, n_train = 20466. Ligand-only control saved beside it. Target: LCK, sequence length 509, 68 labeled ligands, 17 binders at pKd >= 7. LCK is in the scaffold train. This is not a cold target.

Library n = 8068 as in section 3. Ranked file: `artifacts/screen_ranked.csv`.

Recovery of the 17 known LCK binders, all present in the library (`artifacts/metrics_screen.md`):

| Metric | Value |
|---|---|
| Mean rank | 378.4118 |
| Random expected mean rank | 4034.5000 |
| Median rank | 263.0000 |
| In top 50 | 7 (random expected 0.1054) |
| In top 100 | 7 (random expected 0.2107) |
| In top 500 | 15 (random expected 1.0535); enrichment 14.2376x |

The top of the list is mostly DAVIS spike-ins. Top 20 by predicted pKd:

| rank | pred_pkd | ligand_only_pred_pkd | source | smiles |
|---:|---:|---:|---|---|
| 1 | 9.1475 | 5.6711 | DAVIS_spike_in | `Cc1nc(Nc2ncc(C(=O)Nc3c(C)cccc3Cl)s2)cc(N2CCN(CCO)CC2)n1` |
| 2 | 8.3247 | 5.9842 | DAVIS_spike_in | `COc1cc2c(Oc3ccc(NC(=O)C4(C(=O)Nc5ccc(F)cc5)CC4)cc3F)ccnc2cc1OCCCN1CCOCC1` |
| 3 | 8.3230 | 5.6955 | DAVIS_spike_in | `CSc1cccc(Nc2ncc3cc(-c4c(Cl)cccc4Cl)c(=O)n(C)c3n2)c1` |
| 4 | 8.0111 | 5.6600 | MoleculeNet_HIV | `Cc1ccc(C(=O)Nc2ccc(S(=O)(=O)O)c3cc(S(=O)(=O)O)cc(S(=O)(=O)O)c23)cc1NC(=O)c1cccc(NC(=S)Nc2cccc(C(=O)Nc3cc(C(=O)Nc4ccc(S(=O)(=O)O)c5cc(S(=O)(=O)O)cc(S(=O)(=O)O)c45)ccc3C)c2)c1.[NaH]` |
| 5 | 7.9527 | 5.7380 | DAVIS_spike_in | `CCN(CCO)CCCOc1ccc2c(Nc3cc(CC(=O)Nc4cccc(F)c4)[nH]n3)ncnc2c1` |
| 6 | 7.8527 | 5.6665 | MoleculeNet_HIV | `N#Cc1cccc(NC(=O)c2ccc(C(=O)O)c(C(=O)c3ccc(C(=O)O)c(C(=O)Nc4cccc(C#N)c4)c3)c2)c1` |
| 7 | 7.8383 | 5.6446 | MoleculeNet_HIV | `COc1cccc(NC(=O)c2c(SC)c(C#N)c(=O)n(N)c2N)c1` |
| 8 | 7.8130 | 5.6655 | MoleculeNet_HIV | `O=C(CCC(=O)Nc1cccc(Cl)c1)CC(=O)c1ccc(F)cc1` |
| 9 | 7.8033 | 5.6370 | MoleculeNet_HIV | `COc1cccc(NC(=O)C(=O)C2C=CCS2(=O)=O)c1` |
| 10 | 7.7801 | 5.6476 | MoleculeNet_HIV | `COC(=O)NN=C(C(=O)Nc1cccc([N+](=O)[O-])c1)c1nc2ccc([N+](=O)[O-])cc2nc1O` |
| 11 | 7.7433 | 5.6174 | MoleculeNet_HIV | `O=C(Nc1cccc([N+](=O)[O-])c1)C(=NNC(=O)c1ccncc1)c1nc2ccc([N+](=O)[O-])cc2nc1O` |
| 12 | 7.7383 | 5.6701 | MoleculeNet_HIV | `Cc1cc2c(cc1C)[n+]([O-])c(C(=O)CC(=NNC(=S)NN)C(=O)Nc1cccc(C(F)(F)F)c1)c[n+]2[O-]` |
| 13 | 7.7366 | 5.6692 | MoleculeNet_HIV | `O=C(Nc1cccc(-c2nc3ccccc3[nH]2)c1)c1ccc(C(=O)Nc2cccc(-c3nc4ccccc4[nH]3)c2)cc1` |
| 14 | 7.7351 | 5.9222 | DAVIS_spike_in | `CCN1CCN(Cc2ccc(NC(=O)Nc3ccc(Oc4cc(NC)ncn4)cc3)cc2C(F)(F)F)CC1` |
| 15 | 7.7018 | 5.7473 | MoleculeNet_HIV | `Cc1ccc(C(=O)Nc2ccc(CP(=O)(O)O)cc2)cc1NC(=O)c1cccc(NC(=O)Nc2cccc(C(=O)Nc3cc(C(=O)Nc4ccc(CP(=O)(O)O)cc4)ccc3C)c2)c1.[NaH]` |
| 16 | 7.7012 | 5.6462 | MoleculeNet_HIV | `COC(=O)c1cc(C(=O)OC)c(CC(=NNC(=O)C[N+](C)(C)C)C(=O)Nc2cccc(C(F)(F)F)c2)nc1CC(=NNC(=O)C[N+](C)(C)C)C(=O)Nc1cccc(C(F)(F)F)c1.[Cl-]` |
| 17 | 7.6853 | 5.6526 | MoleculeNet_HIV | `O=C(NCCCC1C2CCC(C2)C1CNC(=O)Nc1cccc(Cl)c1)Nc1cccc(Cl)c1` |
| 18 | 7.6849 | 5.7239 | DAVIS_spike_in | `C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1` |
| 19 | 7.6502 | 5.6574 | MoleculeNet_HIV | `N#CCc1ccc(C(C#N)C(CCCC(=O)Nc2cccc([N+](=O)[O-])c2)=NNC(=O)NN)cc1` |
| 20 | 7.6297 | 5.6635 | MoleculeNet_HIV | `O=C(Nc1cccc(O)c1)C(=O)C1C(=O)NC(=S)NC1=O` |

Known binders by rank:

| rank | true_pkd | pred_pkd | drug_id | smiles |
|---:|---:|---:|---|---|
| 1 | 9.6990 | 9.1475 | 3062316 | `Cc1nc(Nc2ncc(C(=O)Nc3c(C)cccc3Cl)s2)cc(N2CCN(CCO)CC2)n1` |
| 2 | 8.2218 | 8.3247 | 42642645 | `COc1cc2c(Oc3ccc(NC(=O)C4(C(=O)Nc5ccc(F)cc5)CC4)cc3F)ccnc2cc1OCCCN1CCOCC1` |
| 3 | 8.9586 | 8.3230 | 447077 | `CSc1cccc(Nc2ncc3cc(-c4c(Cl)cccc4Cl)c(=O)n(C)c3n2)c1` |
| 14 | 7.9586 | 7.7351 | 11409972 | `CCN1CCN(Cc2ccc(NC(=O)Nc3ccc(Oc4cc(NC)ncn4)cc3)cc2C(F)(F)F)CC1` |
| 37 | 8.2076 | 7.5369 | 9809715 | `COC(=O)c1ccc2c(c1)NC(=O)C2=C(Nc1ccc(N(C)C(=O)CN2CCN(C)CC2)cc1)c1ccccc1` |
| 42 | 7.5229 | 7.5259 | 44259 | `CNC1CC2OC(C)(C1OC)n1c3ccccc3c3c4c(c5c6ccccc6n2c5c31)C(=O)NC4` |
| 47 | 7.1427 | 7.5057 | 16722836 | `Cc1cnc(Nc2ccc(OCCN3CCCC3)cc2)nc1Nc1cccc(S(=O)(=O)NC(C)(C)C)c1` |
| 252 | 7.5229 | 6.9911 | 11984591 | `COc1cc(Nc2ncc(F)c(Nc3ccc4c(n3)NC(=O)C(C)(C)O4)n2)cc(OC)c1OC.O=S(=O)(O)c1ccccc1` |
| 263 | 7.2147 | 6.9555 | 5494449 | `Cc1cc(Nc2cc(N3CCN(C)CC3)nc(Sc3ccc(NC(=O)C4CC4)cc3)n2)n[nH]1` |
| 267 | 9.2291 | 6.9518 | 5328940 | `COc1cc(Nc2c(C#N)cnc3cc(OCCCN4CCN(C)CC4)c(OC)cc23)c(Cl)cc1Cl` |
| 342 | 7.3279 | 6.7284 | 644241 | `Cc1cn(-c2cc(NC(=O)c3ccc(C)c(Nc4nccc(-c5cccnc5)n4)c3)cc(C(F)(F)F)c2)cn1` |
| 345 | 7.3979 | 6.7232 | 5291 | `Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1` |
| 352 | 7.5086 | 6.6994 | 10074640 | `Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nc(-c2cccnc2)cs1` |
| 367 | 7.1739 | 6.6632 | 9933475 | `COc1cc2c(Oc3ccc4[nH]c(C)cc4c3F)ncnc2cc1OCCCN1CCCC1` |
| 377 | 7.7696 | 6.6235 | 3081361 | `COc1cc2c(Nc3ccc(Br)cc3F)ncnc2cc1OCC1CCN(C)CC1` |
| 1310 | 7.3098 | 5.9382 | 16038120 | `COc1cc(N2CCC(N3CCN(C)CC3)CC2)ccc1Nc1ncc(Cl)c(Nc2ccccc2S(=O)(=O)C(C)C)n1` |
| 2412 | 7.5229 | 5.7634 | 11626560 | `CC(Oc1cc(-c2cnn(C3CCNCC3)c2)cnc1N)c1c(Cl)ccc(F)c1Cl` |

Two of the 17 binders land at ranks 1310 and 2412. High recovery in the top 500 is not uniform. Ligand-only predicted pKd on the top 20 sits in a narrow band around 5.6 to 6.0, while the full model spreads from 9.15 down. The export says the full model separates known actives more clearly at the top. I do not have a ligand-only recovery table for v1.5. That comparison was not exported. It is the point of v1.6.

This recovery is a sanity check. Scaffold train+val includes some LCK pairs. Scores are predicted pKd.

### 5.8 v1.6: cold-target screen on SLK

Experiment freeze `sofia-vs-v1.6-cold-screen`. Generated 2026-10-03T12:45:44-04:00. Question: if the scored kinase is fully held out, does HistGBM+`aac_dpc` still recover its known binders better than ligand-only?

Selection rule, locked before any model score: exclude LCK; drop `target_id`s that contain "("; require at least 10 binders at pKd >= 7; pick the binder count closest to LCK (17). Tie break: more binders, then `target_id`. Result: SLK, 17 binders. 30 kinases were eligible. Runners-up: GAK (16), LOK (19), DDR2 (15), VEGFR2 (15), FLT1 (20). SLK sequence length 1235. Nearest other kinase by AAC+DPC cosine: NEK1 at 0.9702. All 68 SLK pairs held out. n_train = 25704. Library reused, n = 8068. Top 1% cutoff is rank <= 81.

This is a cold target, not a cold ligand. The same 68 ligands still have labels on other kinases. Ligand-only can rank chemotypes that were potent elsewhere.

Source: `artifacts/metrics_screen_cold.md`.

| model | n binders in library | fraction top 1% | n top 1% | fraction top 500 | n top 500 | mean rank | random mean rank |
|---|---:|---:|---:|---:|---:|---:|---:|
| protein (ECFP+aac_dpc) | 17 | 0.4706 | 8 | 0.5882 | 10 | 1311.2353 | 4034.5000 |
| ligand-only | 17 | 0.4118 | 7 | 0.7059 | 12 | 469.1765 | 4034.5000 |

Random expected counts: top 1% about 0.1707 (fraction 0.0100); top 500 about 1.0535 (fraction 0.0620). Mean rank versus random (random/model, higher is better than chance): protein 3.0769, ligand-only 8.5991.

Does the protein model beat ligand-only? The export's verdict is mixed. Rule used in the export: yes if protein is strictly better on at least one of {top 1% fraction, top 500 fraction, mean rank} and worse on none; no if ligand-only is strictly better on at least one and protein on none; mixed if each wins at least one. Lower mean rank is better. Deltas (positive means protein better): top 1% fraction +0.0588, top 500 fraction -0.1176, mean-rank improvement -842.0588.

Both beat a random ranking. Neither result looks like v1.5, where the LCK mean rank was 378.4118 with the target partly in train. Once the kinase is unseen, known-binder recovery is largely ligand-driven. The protein model wins only the very top of the list (8/17 vs 7/17 inside rank 81).

Known SLK binders, sorted by protein rank:

| protein_rank | ligand_only_rank | true_pkd | pred_pkd | ligand_only_pred_pkd | drug_id | smiles |
|---:|---:|---:|---:|---:|---|---|
| 1 | 1 | 10.6198 | 7.9691 | 7.2437 | 44259 | `CNC1CC2OC(C)(C1OC)n1c3ccccc3c3c4c(c5c6ccccc6n2c5c31)C(=O)NC4` |
| 4 | 45 | 7.4559 | 7.0830 | 6.7909 | 126565 | `CC12OC(CC1(O)CO)n1c3ccccc3c3c4c(c5c6ccccc6n2c5c31)CNC4=O` |
| 46 | 55 | 7.2518 | 6.6076 | 6.0519 | 5329102 | `CCN(CC)CCNC(=O)c1c(C)[nH]c(C=C2C(=O)Nc3ccc(F)cc32)c1C` |
| 49 | 87 | 7.6383 | 6.5307 | 5.9204 | 11409972 | `CCN1CCN(Cc2ccc(NC(=O)Nc3ccc(Oc4cc(NC)ncn4)cc3)cc2C(F)(F)F)CC1` |
| 52 | 221 | 8.3279 | 6.3702 | 5.7678 | 5328940 | `COc1cc(Nc2c(C#N)cnc3cc(OCCCN4CCN(C)CC4)c(OC)cc23)c(Cl)cc1Cl` |
| 61 | 84 | 7.2924 | 6.2457 | 5.9451 | 9809715 | `COC(=O)c1ccc2c(c1)NC(=O)C2=C(Nc1ccc(N(C)C(=O)CN2CCN(C)CC2)cc1)c1ccccc1` |
| 62 | 48 | 7.4815 | 6.2338 | 6.1218 | 11984591 | `COc1cc(Nc2ncc(F)c(Nc3ccc4c(n3)NC(=O)C(C)(C)O4)n2)cc(OC)c1OC.O=S(=O)(O)c1ccccc1` |
| 73 | 75 | 7.8861 | 6.1520 | 6.0025 | 11427553 | `O=C(c1ccc(C=Cc2n[nH]c3ccccc23)cc1)N1CCNCC1` |
| 135 | 47 | 7.6778 | 5.9782 | 6.2726 | 16038120 | `COc1cc(N2CCC(N3CCN(C)CC3)CC2)ccc1Nc1ncc(Cl)c(Nc2ccccc2S(=O)(=O)C(C)C)n1` |
| 208 | 1007 | 7.1805 | 5.9064 | 5.3441 | 9915743 | `CCOc1cc2ncc(C#N)c(Nc3ccc(OCc4ccccn4)c(Cl)c3)c2cc1NC(=O)C=CCN(C)C` |
| 731 | 322 | 7.0757 | 5.5589 | 5.6877 | 25243800 | `CC(C)N1NC(=C2C=c3cc(O)ccc3=N2)c2c(N)ncnc21` |
| 810 | 76 | 8.5229 | 5.5342 | 5.9969 | 42642645 | `COc1cc2c(Oc3ccc(NC(=O)C4(C(=O)Nc5ccc(F)cc5)CC4)cc3F)ccnc2cc1OCCCN1CCOCC1` |
| 1361 | 499 | 7.0862 | 5.4310 | 5.5319 | 5494449 | `Cc1cc(Nc2cc(N3CCN(C)CC3)nc(Sc3ccc(NC(=O)C4CC4)cc3)n2)n[nH]1` |
| 2898 | 1003 | 7.0223 | 5.2976 | 5.3445 | 3081361 | `COc1cc2c(Nc3ccc(Br)cc3F)ncnc2cc1OCC1CCN(C)CC1` |
| 3654 | 513 | 7.7447 | 5.2656 | 5.5245 | 11626560 | `CC(Oc1cc(-c2cnn(C3CCNCC3)c2)cnc1N)c1c(Cl)ccc(F)c1Cl` |
| 4411 | 1526 | 7.6576 | 5.2402 | 5.2962 | 9933475 | `COc1cc2c(Oc3ccc4[nH]c(C)cc4c3F)ncnc2cc1OCCCN1CCCC1` |
| 7735 | 2367 | 7.5850 | 5.1183 | 5.2366 | 176870 | `C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1` |

One known binder (drug_id 176870, true pKd 7.5850) is ranked 7735 by the protein model and 2367 by ligand-only. Predicted pKd for that row is 5.1183. The protein model is not failing only in the middle of the list.

Top 10 by protein score (`artifacts/metrics_screen_cold.md`):

| rank | pred_pkd | ligand_only_rank | ligand_only_pred_pkd | source | known_binder | smiles |
|---:|---:|---:|---:|---|---|---|
| 1 | 7.9691 | 1 | 7.2437 | DAVIS_spike_in | True | `CNC1CC2OC(C)(C1OC)n1c3ccccc3c3c4c(c5c6ccccc6n2c5c31)C(=O)NC4` |
| 2 | 7.8507 | 2 | 7.2421 | MoleculeNet_HIV | False | `CNC1CC2OC(C)(C1OC)n1c3ccccc3c3c4c(c5c6ccccc6n2c5c31)C(=O)NC4O` |
| 3 | 7.2743 | 26 | 7.1078 | MoleculeNet_HIV | False | `CC1(C)C2CCC13CS(=O)(=O)N(C(=O)C1CC=CCC1)C3C2` |
| 4 | 7.0830 | 45 | 6.7909 | DAVIS_spike_in | True | `CC12OC(CC1(O)CO)n1c3ccccc3c3c4c(c5c6ccccc6n2c5c31)CNC4=O` |
| 5 | 7.0353 | 9 | 7.1167 | MoleculeNet_HIV | False | `CN(C)CCCNC(=O)c1cn2c(ccc3sc4ccccc4c(=O)c32)n1` |
| 6 | 6.9399 | 20 | 7.1129 | MoleculeNet_HIV | False | `CCC(C)(O)C(=O)OC1C(=O)OC2CC3C(C)=CC(=O)C(O)C3(C)C3C4(O)OCC23C1C(C)C4O` |
| 7 | 6.9384 | 13 | 7.1162 | MoleculeNet_HIV | False | `CCN1CCN(CC)[Pd-2]12[n+]1ccccc1-c1cccc[n+]12.[O-][Cl+3]([O-])([O-])O` |
| 8 | 6.9004 | 18 | 7.1140 | MoleculeNet_HIV | False | `CCN=C(NC(=O)c1ccc(F)cc1)SCSC(=NC(=O)c1ccc(F)cc1)NCC` |
| 9 | 6.9004 | 4 | 7.1244 | MoleculeNet_HIV | False | `CCN(N=C1CC2C=CCN(C2)C(=O)Cc2c1[nH]c1ccccc21)S(=O)(=O)c1ccc(C)cc1` |
| 10 | 6.8526 | 25 | 7.1096 | MoleculeNet_HIV | False | `Cc1ccc(NCCCN2CCN(CCCNc3ccc(C)c4sc5ccccc5c(=O)c34)CC2)c2c(=O)c3ccccc3sc12` |

Rank 2 is an HIV molecule one heavy atom away from the known binder at rank 1, by inspection of the two SMILES in the export (the HIV row adds an oxygen on the last ring system). I am not calling that a measured SLK binder. The table marks `known_binder` False. Full ranking: `artifacts/screen_ranked_cold.csv`.

## 6. What we can and cannot claim

We can claim, on these freezes and these saved splits:

- Scaffold holdout on DAVIS: a ligand-only ECFP HistGBM does not rank (rho -0.0107, EF@1% 0). Adding protein composition does.
- Cold-protein holdout on DAVIS: AAC does not beat ligand-only. `aac_dpc` and ESM-2 do, by about 0.05 Spearman and by EF@1% (10.2458 vs 8.4528). The strict rho gate was recorded False; the EF gate was True.
- Chemprop D-MPNN matches ECFP ligand-only on cold-protein and does not benefit from `aac_dpc` concat on that split. HistGBM+`aac_dpc` remains the stronger cold-protein ranker in the exports (rho 0.5501, EF@1% 10.2458).
- Isotonic regression, fit on validation, makes cold-protein HistGBM probabilities much closer to observed binder rates. It does not change rank order.
- On KIBA score, `aac_dpc` beats ligand-only on cold-protein Spearman by 0.0609. EF@1% does not improve.
- On SRC versus LCK, with both kinases held out, predicted deltas track true deltas (rho 0.3694) and enrich the 12 selective ligands. Ligand-only deltas are constant.
- A scaffold-trained HistGBM+`aac_dpc` ranks known LCK binders high in an 8068-molecule list (15/17 in the top 500, mean rank 378.4118). That target was not held out.
- With SLK fully held out, both models beat a random rank. The protein model does not beat ligand-only on recovery. The export calls that mixed: protein higher only in the top 1% (8/17 vs 7/17); ligand-only better in the top 500 (12/17 vs 10/17) and on mean rank (469.1765 vs 1311.2353).

We cannot claim:

- A wet-lab hit, a prospective screen, or a NOVA concordance. None of those were run.
- That cold-protein rho near 0.50 is sequence understanding. Ligand-only is already there, because the 68 ligands are shared.
- That v1.5 recovery would survive a cold target. v1.6 is the check, and it does not support that reading.
- A general selectivity model. One pair, 12 ligands, all LCK-preferring.
- That Chemprop probabilities are calibrated on unseen pairs. The low ECE uses a slice of the test set.
- That KIBA metrics are pKd metrics.
- Bit-stable Chemprop reruns. The reproduce note allows a 0.05 fallback.
- Performance on a library other than this HIV subsample plus the DAVIS spike-in, or on a target other than the ones named in each export.

## 7. Limitations

DAVIS is 68 ligands and 379 kinases. The panel is kinase-heavy. Affinities are noisy, and the studied targets are the ones people assay. A model can memorize ligand potency and look good as soon as the ligand has been seen on any kinase.

Cold-protein and cold-target are not cold-ligand. v1.6 states this directly. The ligand-only win on SLK mean rank is the expected failure mode, not a surprise appendix.

`aac_dpc` is a composition vector. SRC and LCK already have cosine 0.9621. SLK's nearest neighbor in that space, NEK1, is at 0.9702. Composition cannot separate close kinases. The selectivity sketch still moved, but the dynamic range of predicted delta is small (mean absolute predicted delta 0.1938 against a true cut of 1.0).

ESM-2 here is the 8M t6 model, mean-pooled, truncated at 1022 residues. SLK is length 1235, so a later ESM screen of SLK would truncate. I did not run ESM on the library.

Chemprop was a short CPU run (15 epochs). A longer run or a different way of fusing protein features might differ. That experiment was not exported, so I do not know. What was exported is a negative result for `x_d` concat on this freeze.

Calibration for Chemprop is optimistically fit. Scaffold ECE for HistGBM stays higher than cold-protein ECE after isotonic regression. I would not threshold scaffold probabilities without saying so.

v1.4 has no SRC-preferring ligand at the cut, so a screen that needed SRC-over-LCK selectivity was not tested. The gate was a sketch gate (rho > 0.1 and EF > 1.2), not a pre-registered clinical bar.

The screen library is a diversity subsample of an HIV set plus every DAVIS ligand. Spike-in makes recovery computable. It also puts the answer keys in the list. HIV labels were ignored; I am not claiming anything about HIV activity. No confidence interval on enrichment was exported.

KIBA is a different assay composite from a different mirror, because the TDC file was blocked. Good for a transfer check. Not a second DAVIS.

Reproduce coverage stops at v1, v1.1, v1.2 (loose tolerance), and v1.4 (reference not re-trained in the session that wrote the script). Calibration, KIBA, v1.5, and v1.6 do not have a matching reproduce script in the tree I used.

## 8. What a follow-up would test

I would not add another protein feature on the same DAVIS cold-protein split until the ligand leak is closed. The next test is a cold-ligand screen: hold out the ligand, not only the kinase, and ask whether `aac_dpc` still beats ligand-only on recovery. v1.6 already shows that a cold kinase is not enough.

Second, repeat the v1.6 protocol on more than SLK. The selection rule and the runners-up are written down (GAK, LOK, DDR2, VEGFR2, FLT1). Those kinases were not scored. One kinase can be a fluke, in either direction.

Third, selectivity needs a pair with both directions populated, and more than one pair. The SRC/LCK sketch can stay as a positive control. It cannot be the claim.

Fourth, if probabilities are going to gate a list, refit Chemprop isotonic regression on a real validation split and report ECE on the untouched test. Ranking screens can keep raw pKd. v1.5 and v1.6 already do.

Fifth, an external assay set, or the NOVA concordance check that `DESIGN.md` left optional, is still undone. Until that exists, the ranked CSV is a model sort of public SMILES. I will not call the top rows candidates for synthesis on the basis of these files.
