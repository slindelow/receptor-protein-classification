# Virtual Screening v1.5 - Screening product metrics

**Experiment freeze:** `sofia-vs-v1.5-screen`  
**Data freeze:** `sofia-vs-v1-2026-10`  
**Generated (America/Toronto):** 2026-10-03T08:26:08-04:00  
**Model:** HistGBM + ECFP(r=2,2048) + aac_dpc  
**Weights:** `artifacts/models_v1.5/hgb_screen_aac_dpc.joblib` (scaffold train+val, n_train=20466)  
**Ligand-only control:** `artifacts/models_v1.5/hgb_screen_ligand_only.joblib`  

Honest scope: these are **model scores** (predicted pKd), not wet-lab hits. No prospective claim.

## Demo target

- **target_id:** `LCK` (DAVIS kinase; present in scaffold train, not a cold protein)
- **sequence length:** 509
- **DAVIS labeled pairs for target:** 68
- **DAVIS binders (pKd >= 7.0):** 17

## Screening library

- **n molecules:** 8068
- **n unique Murcko scaffolds:** 7582
- **source:** MoleculeNet HIV (DeepChem S3) (https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/HIV.csv)
- **raw SHA256:** `9ffa7fe57dc86c342627ee1d5255e937e2ab812393c73c4d16c697022f6e1d22`
- **library CSV:** `data/screen_library/library_v1.5.csv`
- **HIV pool (non-DAVIS) + DAVIS ligand spike-in** for recovery analysis

## Top-20 predicted pKd

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

## Known-binder recovery (DAVIS labels for this target overlapping the library)

- Known binders in library: **17** / 17
- Mean rank of known binders: **378.4118** (random expected mean rank ≈ 4034.5000)
- Median rank: **263.0000**
- In top-50: 7 (random expected ≈ 0.1054)
- In top-100: 7 (random expected ≈ 0.2107)
- In top-500: 15 (random expected ≈ 1.0535); enrichment ≈ 14.2376x

Note: screening weights were fit on scaffold train+val, which includes some LCK-labeled pairs whose scaffolds were not held out. Recovery is a sanity check (known actives should rank high), not a cold-ligand prospective test.

### Known binder detail (by rank)

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

## Artifacts

- Ranked CSV: `artifacts/screen_ranked.csv` (columns: smiles, pred_pkd, ligand_only_pred_pkd, rank, ...)
- Models dir: `artifacts/models_v1.5`
- Library manifest: `data/screen_library/MANIFEST_screen.json`
