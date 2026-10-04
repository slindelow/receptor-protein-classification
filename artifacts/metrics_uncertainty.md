# Uncertainty from saved predictions

Numbers in this file were printed by `scripts/uncertainty_bootstrap.py`. Models were not refit.

Seed 42, 1000 resamples, percentile 95% interval. The two-sided p-value recenters each replicate by the observed delta and uses a plus-one correction. Graph Spearman blocks and the SLK recovery block resample rows. The HistGBM cold-protein block, when computed, resamples held-out kinases. Delta is the protein score minus the ligand-only score unless a row says otherwise.

## What was saved

- NPZ prediction files found: `artifacts/models_v1.2/preds_cold_protein_ligand_only.npz`, `artifacts/models_v1.2/preds_cold_protein_protein.npz`, `artifacts/models_v1.2/preds_scaffold_ligand_only.npz`, `artifacts/models_v1.2/preds_scaffold_protein.npz`.
- Per-pair HistGBM AAC+DPC and ligand-only scores for the DAVIS cold-protein test are in `artifacts/scores_histgbm_cold_protein.csv`. Models were loaded from artifacts/models_v1.1, not refit. Scaffold HistGBM scores and KIBA per-pair scores are still absent. The warm LCK screen was not bootstrapped.
- No per-pair KIBA scores were stored. metrics_kiba.json stores aggregate metrics only. KIBA cold-protein delta Spearman was not resampled.
- artifacts/models_v1.2/preds_*.npz are the Chemprop runs named in metrics_v1.2.json (graph model, joint composition versus ligand-only), not the histogram model.

## Computed comparisons

### DAVIS cold-protein, HistGBM amino-acid plus dipeptide versus ligand-only

Status: computed.

Source: `artifacts/scores_histgbm_cold_protein.csv`.

Resample unit: held-out kinase. n kinases: 76. n pairs: 5168. Seed: 42. n resamples: 1000. Models: loaded from artifacts/models_v1.1 joblibs, not refit. Point estimates matched published metrics_v1.1.json within 0.0001: yes.

Draw the held-out kinases with replacement. Keep every test pair of each drawn kinase (a kinase drawn twice contributes its pairs twice). Recompute pooled Spearman and EF@1% on the concatenated pairs. The same draw enters both scores.

| Quantity | Point | 95% CI low | 95% CI high | Includes 0 | Two-sided p |
|---|---:|---:|---:|---|---:|
| Spearman rho, AAC+DPC | 0.5501 | 0.5118 | 0.5848 | no |  |
| Spearman rho, ligand-only | 0.5002 | 0.4613 | 0.5386 | no |  |
| Delta Spearman (AAC+DPC minus ligand-only) | 0.0500 | 0.0271 | 0.0720 | no | 0.0010 |
| Delta EF@1% (AAC+DPC minus ligand-only) | 1.7930 | -0.2624 | 3.2468 | yes | 0.0909 |

Unrounded delta Spearman: 0.04996776087112098. Unrounded delta EF@1%: 1.7930214115781133. Finite delta-Spearman replicates: 1000.

AAC+DPC was the best cold-protein Spearman among four modes scored on this same split (ligand-only 0.500158, AAC 0.498399, DPC 0.540900, AAC+DPC 0.550125). ESM was a fifth comparison and was essentially tied with AAC+DPC (ESM Spearman 0.547628). The contrast tested here is the pre-specified AAC+DPC versus ligand-only contrast. The interval is not adjusted for the other modes.

The delta Spearman interval excludes zero. That exclusion is not multiplicity-adjusted: AAC+DPC was chosen as the best of four modes on this split, with ESM a fifth comparison essentially tied with it. The delta EF@1% interval includes zero. This file does not call either contrast significant.

### DAVIS scaffold, HistGBM amino-acid plus dipeptide versus ligand-only

Status: skipped.

Skipped. Per-pair scores for this comparison are not in artifacts/.

### KIBA cold-protein, HistGBM joint composition versus ligand-only

Status: skipped.

Skipped. Per-pair KIBA scores are not in artifacts/. The recovery-style or delta-Spearman interval was not imputed from aggregate metrics.

### DAVIS cold-protein, graph model plus joint composition versus graph ligand-only

Status: computed.

Source: `artifacts/models_v1.2/preds_cold_protein_protein.npz and artifacts/models_v1.2/preds_cold_protein_ligand_only.npz`.

| Quantity | Point | 95% CI low | 95% CI high | Includes 0 | Two-sided p |
|---|---:|---:|---:|---|---:|
| Spearman rho, protein | 0.4895 | 0.4683 | 0.5113 | no |  |
| Spearman rho, ligand-only | 0.4997 | 0.4794 | 0.5196 | no |  |
| Delta Spearman (protein minus ligand-only) | -0.0102 | -0.0304 | 0.0104 | yes | 0.3586 |

Rows: 5168. Finite delta replicates: 1000.

### DAVIS scaffold, graph model plus joint composition versus graph ligand-only

Status: computed.

Source: `artifacts/models_v1.2/preds_scaffold_protein.npz and artifacts/models_v1.2/preds_scaffold_ligand_only.npz`.

| Quantity | Point | 95% CI low | 95% CI high | Includes 0 | Two-sided p |
|---|---:|---:|---:|---|---:|
| Spearman rho, protein | 0.3176 | 0.2902 | 0.3412 | no |  |
| Spearman rho, ligand-only | 0.2621 | 0.2370 | 0.2860 | no |  |
| Delta Spearman (protein minus ligand-only) | 0.0556 | 0.0240 | 0.0864 | no | 0.0020 |

Rows: 5306. Finite delta replicates: 1000.

### SLK library recovery, HistGBM joint composition versus ligand-only

Status: computed.

Source: `artifacts/screen_ranked_cold.csv`.

| Quantity | Point | 95% CI low | 95% CI high | Includes 0 | Two-sided p |
|---|---:|---:|---:|---|---:|
| Delta fraction in top 1% (protein minus ligand-only) | 0.0588 | -0.2308 | 0.3158 | yes | 0.7792 |
| Delta fraction in top 500 (protein minus ligand-only) | -0.1176 | -0.3750 | 0.1334 | yes | 0.4306 |
| Mean-rank improvement (ligand-only minus protein; positive favors protein) | -842.0588 | -1718.5427 | -150.3346 | no | 0.0370 |

Continuous per-pair pKd is not stored for the library, so delta Spearman was not computed. The interval is on the recovery difference only.

Library rows: 8068. Known binders: 17. Replicates with no known binder: 0.
