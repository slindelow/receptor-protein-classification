# Calibration / ECE (Sofia VS)

**Generated (America/Toronto):** 2026-10-02T17:25:51-04:00  
**Data freeze:** `sofia-vs-v1-2026-10`  
**Binder threshold:** pKd >= 7.0  
**ECE bins:** 10  

Continuous affinity scores mapped to P(binder) via isotonic regression.
HistGBM: fit on val split. Chemprop: val scores not saved; isotonic on 30% test holdout (documented).

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

## Notes

- Ranking metrics (Spearman/EF) unchanged by calibration; ECE/Brier measure probability quality.
- Honest negative: if ECE stays high after isotonic, scores are poorly ordered for probability.
