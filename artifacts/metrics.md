# Virtual Screening MVP — Metrics

**Freeze:** `sofia-vs-v1-2026-10`  
**Generated (America/Toronto):** 2026-10-02T14:52:41-04:00  
**Binder threshold:** pKd ≥ 7.0  
**Fingerprint:** ECFP radius=2, bits=2048  
**Model:** HistGradientBoostingRegressor  
**Seed:** 42  

All numbers exported by `scripts/run_train_eval.py` / `reproduce.sh`. Do not hand-edit.

## Results

| Split | Model | Spearman ρ | EF@1% | EF@5% | AUROC | n_test | n_binders |
|---|---|---:|---:|---:|---:|---:|---:|
| scaffold | ECFP+AAC | 0.2124 | 6.5072 | 3.3025 | 0.6276 | 5306 | 302 |
| scaffold | ligand-only (ECFP) | -0.0107 | 0.0000 | 0.1982 | 0.4413 | 5306 | 302 |
| cold_protein | ECFP+AAC | 0.4984 | 8.7090 | 6.1712 | 0.8352 | 5168 | 388 |
| cold_protein | ligand-only (ECFP) | 0.5002 | 8.4528 | 6.5312 | 0.8280 | 5168 | 388 |

## Notes

- **Scaffold split:** Murcko scaffolds held out (test ≈ 20%, val ≈ 10%).
- **Cold-protein split:** entire target proteins held out (test ≈ 20%, val ≈ 10%).
- **Ligand-only ablation:** same splits, protein AAC features removed — control for ligand memorization.
- Reproduce float tolerance: `0.0001` (absolute) on Spearman/AUROC/EF.

## Artifacts

- `cold_protein_full`: `artifacts/models/hgb_cold_protein_full.joblib`
- `cold_protein_ligand_only`: `artifacts/models/hgb_cold_protein_ligand_only.joblib`
- `scaffold_full`: `artifacts/models/hgb_scaffold_full.joblib`
- `scaffold_ligand_only`: `artifacts/models/hgb_scaffold_ligand_only.joblib`
- `screen`: `artifacts/models/hgb_screen.joblib`
