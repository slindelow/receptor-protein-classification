# Virtual Screening v1.1 protein feature ablation

**Experiment freeze:** `sofia-vs-v1.1-protein`  
**Data freeze:** `sofia-vs-v1-2026-10`  
**Generated (America/Toronto):** 2026-10-02T16:35:33-04:00  
**Binder threshold:** pKd ≥ 7.0  
**Fingerprint:** ECFP radius=2, bits=2048  
**Model:** HistGradientBoostingRegressor  
**Seed:** 42  
**ESM-2:** ran (esm2_t6_8M_UR50D, dim=320, cache=data/cache/esm2_t6)  

All numbers exported by `scripts/run_protein_ablation.py` / `reproduce_v1.1.sh`. Do not hand-edit.
Splits loaded from existing `data/splits/*.csv` (not regenerated).

## Results

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

## Success criteria (auto)

- Primary (cold_protein): best protein mode=`aac_dpc` ρ=0.5501 vs ligand-only ρ=0.5002 (Δ=0.0500). EF@1% 10.2458 vs 8.4528; EF@5% 6.9941 vs 6.5312.
- Primary pass (ΔSpearman≥0.05 OR clear EF lift): YES (rho_gate=False, ef_gate=True).
- Secondary scaffold `aac` ≥ ligand-only: YES (ρ=0.2124 vs -0.0107).
- Secondary scaffold `dpc` ≥ ligand-only: YES (ρ=0.2219 vs -0.0107).
- Secondary scaffold `aac_dpc` ≥ ligand-only: YES (ρ=0.2397 vs -0.0107).
- Secondary scaffold `esm2` ≥ ligand-only: YES (ρ=0.2566 vs -0.0107).

## Notes

- Same raw DAVIS file + MANIFEST hashes as MVP data freeze.
- Same split CSVs under `data/splits/` (loaded, not rewritten).
- Ligand-only control shared across protein modes (trained once per split).
- HistGBM hyperparams and seeds unchanged from MVP.
- Reproduce float tolerance: `0.0001` (absolute).

## Artifacts

- `cold_protein_full_aac`: `artifacts/models_v1.1/hgb_cold_protein_full_aac.joblib`
- `cold_protein_full_aac_dpc`: `artifacts/models_v1.1/hgb_cold_protein_full_aac_dpc.joblib`
- `cold_protein_full_dpc`: `artifacts/models_v1.1/hgb_cold_protein_full_dpc.joblib`
- `cold_protein_full_esm2`: `artifacts/models_v1.1/hgb_cold_protein_full_esm2.joblib`
- `cold_protein_ligand_only`: `artifacts/models_v1.1/hgb_cold_protein_ligand_only.joblib`
- `scaffold_full_aac`: `artifacts/models_v1.1/hgb_scaffold_full_aac.joblib`
- `scaffold_full_aac_dpc`: `artifacts/models_v1.1/hgb_scaffold_full_aac_dpc.joblib`
- `scaffold_full_dpc`: `artifacts/models_v1.1/hgb_scaffold_full_dpc.joblib`
- `scaffold_full_esm2`: `artifacts/models_v1.1/hgb_scaffold_full_esm2.joblib`
- `scaffold_ligand_only`: `artifacts/models_v1.1/hgb_scaffold_ligand_only.joblib`
