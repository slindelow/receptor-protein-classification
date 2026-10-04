# Ranking Public Kinase-Ligand Affinities With and Without Protein Features

Sofia Lindelow. 2026-10-03.

A score for one protein-ligand pair was tested for whether it can order a public compound library for one kinase, and for whether that order uses the protein or only the ligand. On the saved splits, amino-acid plus dipeptide composition beats a ligand-only control when scaffolds are held out, and it adds a smaller lift when whole proteins are held out (cold-protein Spearman rho 0.5501 versus 0.5002; enrichment at 1% 10.2458 versus 8.4528). A graph model did not keep that cold-protein lift. The composition lift transferred to KIBA on Spearman rho, not on enrichment at 1%. A warm library rank for LCK recovered known binders far above chance, and a cold rank for SLK did not beat ligand-only on the broader ranks. The scores are predictions on public labels. They are not wet-lab hits.

The full write-up, with every number taken from the exported tables, is [RESEARCH_REPORT.md](RESEARCH_REPORT.md).

## What held up

The working scorer is a histogram gradient booster on Morgan fingerprints, plus amino-acid composition and dipeptide composition. On the saved DAVIS splits, those protein features beat a ligand-only control when scaffolds are held out, and they add a smaller lift when whole proteins are held out (cold-protein Spearman 0.5501 vs 0.5002). A Chemprop graph model did not keep that cold-protein lift. The same features transferred to KIBA on Spearman, not on enrichment at 1%. A one-pair selectivity check (SRC vs LCK) was positive, and ligand-only could not rank it. A warm screen of LCK looked strong because that kinase was in training. A cold screen of SLK stayed above chance, and ligand-only won the broader ranks.

## Reproduce

```bash
uv venv --python 3.12 .venv && source .venv/bin/activate
uv pip install -e .
./reproduce.sh
```

`reproduce.sh` checks the v1 metrics within 1e-4. Later shell entry points that exist are `reproduce_v1.1.sh`, `reproduce_v1.2.sh`, and `reproduce_v1.4.sh`. There is no `reproduce_v1.3.sh`, `reproduce_v1.5.sh`, or `reproduce_v1.6.sh`. The matching Python scripts in the tree are `scripts/run_kiba_v1.3.py`, `scripts/run_calibration.py`, `scripts/run_screen_v1.5.py`, and `scripts/run_screen_v1.6.py`. Raw DAVIS tables, the KIBA pair table, the ESM cache, and the KIBA split CSVs are on disk in this working copy and are listed in `.gitignore`, so a public commit should not include them. Manifests under `data/` record the source files and hashes. Rebuilding the KIBA splits needs the frozen pair table, which is also gitignored because it is over GitHub's 100 MB file limit.

## What is not claimed

No new binder. No prospective screen. No external concordance check. See the limits section of the report.

## Verify

From the project root, with the existing `.venv`:

```bash
.venv/bin/python scripts/uncertainty_bootstrap.py
```

That command does not refit any model. It reads saved scores only and rewrites `artifacts/metrics_uncertainty.json` and `artifacts/metrics_uncertainty.md`.

The interval can be recomputed for the graph model, because `artifacts/models_v1.2/preds_cold_protein_protein.npz`, `preds_cold_protein_ligand_only.npz`, `preds_scaffold_protein.npz`, and `preds_scaffold_ligand_only.npz` store per-pair scores. It can also be recomputed for the SLK recovery comparison, because `artifacts/screen_ranked_cold.csv` stores both scores and the known-binder flag.

It cannot be recomputed for the histogram amino-acid-plus-dipeptide versus ligand-only comparison on DAVIS, cold-protein or scaffold. Those per-pair scores were not saved. `artifacts/metrics_v1.1.json` has the aggregate Spearman rho and enrichment only. The same is true for KIBA: `artifacts/metrics_kiba.json` has aggregates, not per-pair scores. The script skips those comparisons and does not impute an interval. The cold-protein histogram difference of 0.0500 and the KIBA cold-protein difference of 0.0609 stay point estimates.

`./reproduce_v1.2.sh` retrains the graph model and can replace the npz files above. Run the read-only checks first if those saved scores should stay as they are.

A reviewer recomputes the ranking math in this order:

1. `src/vs/metrics_lib.py` is the Spearman rho and enrichment definition used by the training scripts.
2. `./reproduce.sh` runs `scripts/run_train_eval.py` and `scripts/check_metrics_tolerance.py`, refits the first DAVIS models, and writes `artifacts/metrics.json`.
3. `./reproduce_v1.1.sh` runs `scripts/run_protein_ablation.py` and refits the histogram protein-feature comparison. It writes `artifacts/metrics_v1.1.json`. Per-pair predictions for that fit were not saved. `scripts/score_histgbm_cold_protein.py` does not refit. It loads the saved v1.1 joblibs, scores the frozen cold-protein test, and writes `artifacts/scores_histgbm_cold_protein.csv` only if Spearman rho and enrichment at 1% match the exported point estimates within 1e-4. It does not rewrite `artifacts/metrics_v1.1.json`. It needs `data/raw/davis_pairs.csv` and `artifacts/models_v1.1/`, which `.gitignore` excludes.

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
7. `.venv/bin/python scripts/uncertainty_bootstrap.py` writes the bootstrap intervals from the saved graph scores and the saved SLK table.
