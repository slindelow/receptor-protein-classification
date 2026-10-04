# Virtual Screening v1.2 Chemprop (D-MPNN)

**Experiment freeze:** `sofia-vs-v1.2-chemprop`  
**Data freeze:** `sofia-vs-v1-2026-10`  
**Generated (America/Toronto):** 2026-10-02T17:06:25-04:00  
**Binder threshold:** pKd >= 7.0  
**Ligand encoder:** Chemprop BondMessagePassing (D-MPNN)  
**Protein extras:** `aac_dpc` as Chemprop x_d concat  
**Seed:** 42  
**Max epochs:** 15  
**Install:** chemprop==2.2.1  

All numbers exported by `scripts/run_chemprop_v1.2.py` / `reproduce_v1.2.sh`. Do not hand-edit.
Splits loaded from existing `data/splits/*.csv` (not regenerated).

## Results (Chemprop)

| Split | Model | Spearman rho | EF@1% | EF@5% | AUROC | n_test | n_binders |
|---|---|---:|---:|---:|---:|---:|---:|
| scaffold | Chemprop ligand-only | 0.2621 | 0.3254 | 0.5284 | 0.6645 | 5306 | 302 |
| scaffold | Chemprop+aac_dpc | 0.3176 | 7.4833 | 4.5575 | 0.7373 | 5306 | 302 |
| cold_protein | Chemprop ligand-only | 0.4997 | 8.7090 | 6.6341 | 0.8279 | 5168 | 388 |
| cold_protein | Chemprop+aac_dpc | 0.4895 | 8.4528 | 5.9655 | 0.8541 | 5168 | 388 |

## v1.1 HistGBM baselines (cited; not retrained)

| Split | Model | Spearman rho | EF@1% | EF@5% | AUROC |
|---|---|---:|---:|---:|---:|
| scaffold | HistGBM ligand-only (ECFP) | -0.0107 | 0.0000 | 0.1982 | 0.4413 |
| scaffold | HistGBM+AAC+DPC | 0.2397 | 5.2058 | 3.6328 | 0.6437 |
| cold_protein | HistGBM ligand-only (ECFP) | 0.5002 | 8.4528 | 6.5312 | 0.8280 |
| cold_protein | HistGBM+AAC+DPC | 0.5501 | 10.2458 | 6.9941 | 0.8689 |

## Success criteria (auto)

- Primary (cold_protein Chemprop+protein vs Chemprop ligand-only): rho 0.4895 vs 0.4997 (Delta=-0.0102); EF@1% 8.4528 vs 8.7090 (Delta=-0.2561); EF@5% 5.9655 vs 6.6341 (Delta=-0.6686).
- Primary pass (DeltaSpearman>=0.05 OR clear EF lift): NO.
- Secondary vs HistGBM+aac_dpc cold_protein: Chemprop+prot rho=0.4895 vs HistGBM rho=0.5501 (Delta=-0.0606).

## Notes

- Same raw DAVIS + MANIFEST hashes; same split CSVs (loaded, not rewritten).
- Protein mode for x_d: `aac_dpc` (best cheap v1.1 mode).
- Chemprop/Lightning CPU runs are not guaranteed bit-stable; seed fixed; see reproduce_v1.2.sh.
- Reproduce float tolerance (strict): `0.0001`; nondeterministic fallback: `0.05`.

## Artifacts

- `cold_protein_ligand_only`: `artifacts/models_v1.2/chemprop_cold_protein_ligand_only.ckpt`
- `cold_protein_ligand_only_preds`: `artifacts/models_v1.2/preds_cold_protein_ligand_only.npz`
- `cold_protein_protein`: `artifacts/models_v1.2/chemprop_cold_protein_protein.ckpt`
- `cold_protein_protein_preds`: `artifacts/models_v1.2/preds_cold_protein_protein.npz`
- `scaffold_ligand_only`: `artifacts/models_v1.2/chemprop_scaffold_ligand_only.ckpt`
- `scaffold_ligand_only_preds`: `artifacts/models_v1.2/preds_scaffold_ligand_only.npz`
- `scaffold_protein`: `artifacts/models_v1.2/chemprop_scaffold_protein.ckpt`
- `scaffold_protein_preds`: `artifacts/models_v1.2/preds_scaffold_protein.npz`
