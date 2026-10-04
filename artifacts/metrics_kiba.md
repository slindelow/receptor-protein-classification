# Virtual Screening KIBA secondary freeze (v1.3)

**Experiment freeze:** `sofia-vs-kiba-v1.3`  
**Generated (America/Toronto):** 2026-10-02T17:26:49-04:00  
**Binder threshold:** KIBA score >= 12.1  
**Label:** KIBA score stored in `pkd` column for protocol reuse (NOT true pKd).  
**Model:** HistGBM ECFP +/- AAC+DPC (same hyperparams as DAVIS MVP).  

| Split | Model | Spearman rho | EF@1% | EF@5% | AUROC | n_test | n_binders |
|---|---|---:|---:|---:|---:|---:|---:|
| scaffold | HistGBM ligand-only | 0.3394 | 1.6314 | 2.2196 | 0.6671 | 23663 | 4529 |
| scaffold | HistGBM+AAC+DPC | 0.5584 | 4.9602 | 3.8171 | 0.7704 | 23663 | 4529 |
| cold_protein | HistGBM ligand-only | 0.4806 | 3.9041 | 3.3085 | 0.7485 | 24280 | 5912 |
| cold_protein | HistGBM+AAC+DPC | 0.5415 | 3.8872 | 3.4980 | 0.7688 | 24280 | 5912 |

## Success (cold_protein protein vs ligand-only)
- Delta Spearman = 0.0609
- Pass Delta>=0.05 or EF lift: YES

## Notes

- DAVIS splits untouched. KIBA splits under `data/splits_kiba/`.
- No wet-lab / NOVA claims.
