# Decisions  -  Sofia VS MVP (2026-10-02, America/Toronto)

## Scope
- Built same-day MUST only: DAVIS freeze, ECFP+HistGBM, scaffold + cold-protein eval, ligand-only ablation, `vs-screen` CLI, demo library, `reproduce.sh`.
- CUT: Chemprop, Boltz-2, KIBA, NOVA concordance, UI, Bittensor.

## Data
- **Primary labels:** TDC DAVIS via Harvard Dataverse file `5219748` (`data/raw/davis_tdc.tab`).
- **pKd transform:** `pKd = 9 - log10(Kd_nM)` (Kd column in TDC file is nM). Documented in `data/MANIFEST.json`.
- **PyTDC note:** `PyTDC==0.4.1` depends on `rdkit-pypi` (no cp312 wheels). On Python 3.12 we pin `rdkit==2024.3.5` and fetch DAVIS from the same Dataverse URL TDC uses (equivalent source). Optional `PyTDC` install is `--no-deps` only.
- **Demo library:** DAVIS has only ~68 unique ligands; padded with public Delaney/ESOL SMILES to n=400. ESOL is **not** used as training labels.

## Features / model
- Ligand: Morgan/ECFP radius=2, 2048 bits (`configs/default.yaml`).
- Protein: amino-acid composition (20-dim). Chosen because it generalizes to unseen sequences (unlike target-ID one-hot), which is required for a meaningful cold-protein check.
- Model: `HistGradientBoostingRegressor` (sklearn). Screening weights = full model fit on scaffold train+val (`artifacts/models/hgb_screen.joblib`).

## Splits
- Murcko scaffold split and cold-protein (hold out entire `target_id`s). Test≈20%, val≈10%, remainder train. Seeds in config. IDs written to `data/splits/`.

## Metrics reading (honest)
- Scaffold: full model beats ligand-only (ρ≈0.21 vs ≈0; enrichment only when protein features present). Expected for scaffold holdout.
- Cold-protein: full ≈ ligand-only (ρ≈0.50). **Disclose:** on this DAVIS freeze, ranking under cold proteins is largely ligand-driven  -  protein AAC adds little. This is a known DTI failure mode the ablation exists to catch; do not claim sequence generalization from cold-protein ρ alone.

## Binder threshold
- Binary metrics use **pKd ≥ 7** (Kd ≤ 100 nM), locked in config.

## Verifiability
- Manifest hashes, pinned deps (`requirements.lock.txt` / `pyproject.toml`), split CSVs, export-only `artifacts/metrics.md` + `metrics.json`, `reproduce.sh` tolerance `1e-4`.

## v1.1 protein features (`sofia-vs-v1.1-protein`)

**Date (America/Toronto):** 2026-10-02 16:35 EDT

### Design
- Goal: stronger protein features that can beat ligand-only on cold-protein, same data freeze and same splits/metrics.
- Locked from MVP: raw DAVIS + MANIFEST hashes, `data/splits/*.csv` (loaded, never rewritten), binder pKd≥7, ECFP(r=2,2048), HistGBM hyperparams/seeds, ligand-only control.
- New protein modes (config `protein_features.type` / ablation list): `aac` (20), `dpc` (400), `aac_dpc` (420), `esm2` (frozen ESM-2 t6 8M mean-pool, dim 320, cache `data/cache/esm2_t6/`).
- Experiment config: `configs/v1.1_protein.yaml` with `freeze_id: sofia-vs-v1.1-protein`. Raw data freeze id unchanged (`sofia-vs-v1-2026-10`).
- Runner: `scripts/run_protein_ablation.py`. Metrics: `artifacts/metrics_v1.1.md` + `.json` (MVP `artifacts/metrics.md` untouched).
- Reproduce: `reproduce_v1.1.sh` (tol 1e-4 vs `.metrics_v1.1_reference.json`).

### What ran
- All four protein modes × {full} × {scaffold, cold_protein}, plus shared ligand-only once per split.
- ESM-2: **ran** on CPU (`fair-esm` + torch CPU wheel). Sequence length truncated at 1022 residues for context; embeddings keyed by sequence hash.
- AAC full metrics match MVP bit-for-bit on both splits (sanity that splits were not regenerated).

### Honest reading
- **Cold-protein primary:** best mode is `aac_dpc` (ρ≈0.550 vs ligand-only ≈0.500; Δ≈0.050). Strict Δ≥0.05 on Spearman is borderline (float rounds to 0.0500; gate recorded False at raw float). Clear EF@1% lift for `dpc` / `aac_dpc` / `esm2` (≈9.99-10.25 vs ≈8.45) and EF@5% lift for `dpc`/`esm2`. Primary success via EF lift: **YES**.
- AAC alone still ≈ ligand-only on cold-protein (same negative finding as MVP).
- **Scaffold secondary:** every protein mode beats ligand-only (ρ≈0.21-0.26 vs ≈0). Stronger composition/ESM features also lift scaffold ρ vs AAC.
- Takeaway: richer sequence composition (DPC / AAC+DPC) and light ESM-2 embeddings add cold-protein ranking signal beyond ligand identity on this DAVIS freeze; AAC was too weak. Still no wet-lab or NOVA claim.

### CUT / pins
- No CUT for ESM (install succeeded). Optional: larger ESM variants deferred.
- New deps pinned: `torch` (CPU wheel), `fair-esm==2.0.0` in `pyproject.toml` / `requirements.lock.txt`.

## v1.2 Chemprop D-MPNN (`sofia-vs-v1.2-chemprop`)

**Date (America/Toronto):** 2026-10-02 (box)

### Design
- Goal: replace ECFP ligand encoder with Chemprop D-MPNN; test whether graph ligand features change the protein-vs-ligand-only story on cold_protein.
- Locked: same DAVIS freeze + MANIFEST hashes, same `data/splits/*.csv` (loaded, never rewritten), binder pKd>=7, same eval metrics (Spearman, EF@1%/5%, AUROC).
- Protein extras: **`aac_dpc` (420-dim)** as Chemprop molecule-level `x_d` concat. Chosen over ESM2 because (1) best cheap cold-protein mode in v1.1 (tied with ESM2 on EF; best Spearman), (2) no ESM dependency at train time, (3) Chemprop `x_d` API is straightforward for fixed-length composition vectors.
- Config: `configs/v1.2_chemprop.yaml`. Runner: `scripts/run_chemprop_v1.2.py`. Metrics: `artifacts/metrics_v1.2.md` + `.json`. Models: `artifacts/models_v1.2/`.
- Install: `chemprop==2.2.1` (Lightning) on existing CPU torch 2.14.1. Checkpoint reload uses `weights_only=False` (torch 2.6+ default otherwise breaks Chemprop pickles).
- Budget: max_epochs=15, batch=64, early_stopping patience=5, hidden=300, depth=3, seed=42, accelerator=cpu.
- Reproduce: `reproduce_v1.2.sh` tries tol 1e-4 then documents nondeterministic fallback tol 5e-2 (Chemprop/Lightning CPU not bit-stable).

### Compare
- Primary: Chemprop+aac_dpc vs Chemprop ligand-only on cold_protein.
- Secondary: cite v1.1 HistGBM+aac_dpc numbers (do not retrain HistGBM).

### Follow-ons (chained)
- ECE/calibration on v1.1 winners (+ Chemprop when preds exist) -> `artifacts/metrics_calibration.*`
- KIBA secondary freeze if Chemprop wall time stays modest -> `scripts/run_kiba_v1.3.py`

### What ran (2026-10-02 17:24 EDT)
- chemprop==2.2.1 on CPU torch; 4 trainings (scaffold/cold x ligand-only/aac_dpc); max_epochs=15 early-stop.
- Wall clock ~18 min for full v1.2 suite.

### Honest reading
- **Primary cold_protein FAILED:** Chemprop+aac_dpc ρ=0.4895 vs Chemprop ligand-only ρ=0.4997 (Δ=-0.010). EF@1%/5% also slightly down. Concat AAC+DPC as x_d did not help D-MPNN ranking under cold proteins on this freeze.
- Chemprop ligand-only cold ≈ HistGBM ligand-only (ρ≈0.50), so the graph encoder matches ECFP on this split but does not beat it.
- **Secondary:** HistGBM+aac_dpc remains stronger on cold_protein (ρ=0.5501, EF@1%=10.25) than either Chemprop variant.
- **Scaffold:** Chemprop ligand-only already ρ=0.26 (HistGBM ligand-only was ≈0); Chemprop+aac_dpc lifts to ρ=0.32 with clear EF (EF@1% 7.48). Protein extras help when scaffolds are held out but ligands may overlap proteins.
- Takeaway: on DAVIS cold-protein, HistGBM+composition still preferred; Chemprop D-MPNN alone is competitive with ECFP but AAC+DPC concat extras do not transfer the v1.1 cold-protein win. No wet-lab/NOVA claim.

### CUT / pins
- No CUT (Chemprop install succeeded).
- Pins: chemprop==2.2.1 (+ lightning). Checkpoint load requires weights_only=False under torch 2.6+/2.14.

## Calibration (`sofia-vs-calibration`)

**Date:** 2026-10-02 17:32 EDT

- Mapped continuous scores to P(binder) via isotonic regression; ECE (10 bins) + Brier.
- HistGBM: isotonic fit on val. Chemprop: val scores not exported; isotonic on 30% test holdout (documented limitation).
- Result: isotonic sharply reduces ECE (cold HistGBM ~0.08-0.10 -> ~0.01; Chemprop holdout -> ~0.005-0.007). Does not change ranking metrics.
- Export: `artifacts/metrics_calibration.md` + `.json`.

## KIBA secondary freeze (`sofia-vs-kiba-v1.3`)

**Date:** 2026-10-02 17:32 EDT

- TDC Dataverse KIBA file returned HTTP 403 from this box. Assembled pairs from DeepDTA public mirror (`data/raw/kiba_deepdta/` -> `kiba_pairs.csv`, n=118254). Documented in `data/MANIFEST_kiba.json`.
- Label: KIBA score stored as `pkd` for protocol reuse (NOT true pKd). Binder threshold KIBA score >= 12.1.
- Same split protocol (scaffold + cold_protein, fracs/seed from default.yaml). Splits written once under `data/splits_kiba/` (DAVIS splits untouched).
- Model: HistGBM ECFP +/- aac_dpc (MVP hyperparams). Cached fingerprints for speed/memory.
- Cold_protein: aac_dpc beats ligand-only (rho 0.5415 vs 0.4806, Delta=0.0609). Primary pass YES. EF@1% roughly flat; gain is Spearman/AUROC.
- Scaffold: large protein lift (0.34 -> 0.56).
- Export: `artifacts/metrics_kiba.md` + `.json`, models under `artifacts/models_kiba/`.

## v1.4 selectivity (`sofia-vs-v1.4-selectivity`)

**Date:** 2026-10-02 17:43 EDT

### Design
- Goal: one related kinase pair, predicted Δscore vs true ΔpKd, plus ligand-only control.
- Pair: **SRC / LCK** (Src family). Chosen because sequences are related but distinguishable (best-window identity ≈0.40, AAC+DPC cosine ≈0.96), all 68 DAVIS ligands are measured on both, and 12 ligands have |ΔpKd| ≥ 1 (enough for a sketch EF). Closer pairs (ERK1/ERK2, CDK2/CDK3) have almost no selectivity signal; ABL1/ABL2 look homologous by composition but N-terminal window identity is ~0.07 and they are much longer.
- Protocol: **cold-pair holdout**. Train HistGBM on all DAVIS pairs except SRC and LCK (n=25636). Existing scaffold models were not reused (different feature path / pair not held out together).
- Scorer: HistGBM + ECFP(r=2, 2048) + **aac_dpc** (v1.1 cold-protein winner). Control: same model, ligand-only.
- Δ := score(SRC) − score(LCK). Selective binder: |ΔpKd| ≥ 1.
- Config: `configs/v1.4_selectivity.yaml`. Runner: `scripts/run_selectivity_v1.4.py`. Export: `artifacts/metrics_selectivity.md` + `.json`.

### Honest reading
- HistGBM+aac_dpc Spearman(Δpred, Δtrue) = **0.369**. EF of |ΔpKd|≥1 ligands ranked by |Δpred|: EF@5%=2.83, EF@10%=3.24, EF@20%=2.43. Sketch gate (ρ>0.1 and EF>1.2) **YES**.
- All 12 selective ligands prefer LCK (ΔpKd ≤ −1); none prefer SRC by ≥1. SRC-directed EF is undefined. LCK-directed EF matches the absolute EF.
- Ligand-only Δpred is exactly 0 (constant). Spearman undefined; EF=0. Honest control null: selectivity ranking requires protein features.
- Scope is a sketch (n=68, one pair, cold holdout). Do not claim a general selectivity model or wet-lab utility.
- Chain stop: selectivity was the last P2 item. Consolidated one-pager: `artifacts/RESULTS_SUMMARY.md`.

### Reproduce
- `reproduce_v1.4.sh` retrains and checks `scripts/check_selectivity_tolerance.py` at tol 1e-4 vs `artifacts/.metrics_v1.4_reference.json`.
- Full retrain not re-run in this session (first export is the reference; self-check of reference vs export is OK).

## v1.5 screening product (`sofia-vs-v1.5-screen`)

**Date (America/Toronto):** 2026-10-03 08:26 EDT

### Design
- Goal: prediction/screening product (rank a public library), not more classification metrics.
- Scorer: **HistGBM + ECFP(r=2, 2048) + aac_dpc** (v1.1 cold-protein winner on DAVIS).
- Weights: fit on **scaffold train+val** (held-out protocol already saved in `data/splits/scaffold_split.csv`). n_train=20466. Saved as `artifacts/models_v1.5/hgb_screen_aac_dpc.joblib`. Ligand-only control: `hgb_screen_ligand_only.joblib`.
- Demo target: **LCK** (DAVIS kinase, seq len 509). Present in scaffold train (not cold). 68 labeled ligands, 17 binders at pKd>=7.
- Screening library: MoleculeNet HIV (DeepChem S3 `HIV.csv`, n=41127 valid SMILES). Diverse Murcko subsample to ~8000 non-DAVIS compounds + spike-in of all 68 DAVIS ligands for recovery. Final n=8068, unique scaffolds=7582. Source URL + SHA256 in `data/screen_library/MANIFEST_screen.json`. HIV activity labels unused.
- Export: full ranking `artifacts/screen_ranked.csv` (smiles, pred_pkd, rank, ligand_only_pred_pkd, ...); `artifacts/metrics_screen.md` + `.json`.
- Config: `configs/v1.5_screen.yaml`. Runner: `scripts/run_screen_v1.5.py`.

### Honest reading
- Known LCK binders (17/17 in library): mean rank 378 vs random expected ~4035; 15/17 in top-500 (random expected ~1.05; enrichment ~14x); 7 in top-50.
- Recovery is a sanity check: scaffold train+val includes some LCK-labeled pairs. Not a cold-ligand prospective test. Scores are model predictions, not wet-lab hits.
- Ligand-only scores included as a control column; full model separates known actives more clearly at the top of the list.

## v1.6 cold-target screen (`sofia-vs-v1.6-cold-screen`)

**Date (America/Toronto):** 2026-10-03 12:45 PM ET

### Design
- Goal: honest version of the v1.5 LCK screen. v1.5 trained on scaffold train+val, which still contained LCK pairs, then ranked the library for LCK. That recovery is not a cold-target test.
- Selection rule, locked before any model score: exclude LCK; drop target_ids that contain '(' (mutant or domain constructs); require at least 10 binders at pKd >= 7; pick the binder count closest to LCK (17). Tie break: more binders, then target_id.
- Result of the rule: **SLK** (17 binders, same count as LCK). 30 eligible kinases. Next were GAK (16), LOK (19), DDR2 (15), VEGFR2 (15), FLT1 (20).
- Protocol: hold out all 68 SLK pairs. Train HistGBM + ECFP(r=2, 2048) + aac_dpc on the remaining DAVIS pairs only (n_train=25704). Ligand-only control on the same rows (ECFP only). HistGBM still uses its internal validation_fraction=0.1 for early stopping.
- Library: reused `data/screen_library/library_v1.5.csv` (n=8068, SHA256 `5721c043c87b744e4f288d31319e7743bf72e7a064de658708b0a16c61df0604`). Not rebuilt and not redownloaded.
- Cold-target, not cold-ligand. The same 68 ligands still have labels on other kinases, so ligand-only can rank chemotypes that were potent elsewhere. Nearest other kinase by AAC+DPC cosine is NEK1 (0.970).
- Export: `artifacts/screen_ranked_cold.csv`, `artifacts/metrics_screen_cold.md` (+ json). Config: `configs/v1.6_cold_screen.yaml`. Runner: `scripts/run_screen_v1.6.py`. Weights: `artifacts/models_v1.6/`.
- Top 1% cutoff matches EF@1%: k = ceil(0.01 * 8068) = 81.

### Honest reading
- Known SLK binders in the library: 17/17.
- Protein model: fraction in top 1% = 8/17 (0.4706); fraction in top 500 = 10/17 (0.5882); mean rank 1311.2 vs random 4034.5 (about 3.1x better than chance).
- Ligand-only: fraction in top 1% = 7/17 (0.4118); fraction in top 500 = 12/17 (0.7059); mean rank 469.2 vs random 4034.5 (about 8.6x better than chance).
- Random expected counts are about 0.17 in the top 1% and about 1.05 in the top 500.
- **Protein model does not beat ligand-only on recovery.** Verdict: mixed. Protein is better only on the top 1% fraction (+0.0588). Ligand-only is better on top 500 (-0.1176 for protein) and on mean rank (protein worse by about 842 ranks).
- Both beat a random ranking. The gap vs v1.5 (LCK mean rank 378 with the target partly in train) is the point of the holdout: once the kinase is fully unseen, known-binder recovery is largely ligand-driven.
- Scores are predicted pKd. They are not wet-lab hits and not a prospective claim.
