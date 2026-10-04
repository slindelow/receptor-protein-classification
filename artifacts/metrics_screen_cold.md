# Virtual Screening v1.6 - cold-target screen

**Experiment freeze:** `sofia-vs-v1.6-cold-screen`  
**Data freeze:** `sofia-vs-v1-2026-10`  
**Generated (America/Toronto):** 2026-10-03T12:45:44-04:00  
**Model:** HistGBM + ECFP(r=2,2048) + aac_dpc  
**Weights:** `artifacts/models_v1.6/hgb_cold_screen_aac_dpc.joblib` (cold-target train, n_train=25704)  
**Ligand-only control:** `artifacts/models_v1.6/hgb_cold_screen_ligand_only.joblib` (n_train=25704)  

Honest scope: these are **model scores** (predicted pKd), not wet-lab hits. No prospective claim.

## Why this is not the v1.5 LCK screen

v1.5 fit HistGBM on scaffold train+val, which still contained LCK pairs, then ranked a library for LCK. Recovery there is not a cold-target test.
v1.6 holds out every pair for one other kinase and trains only on the remaining DAVIS pairs.

## Held-out target

- **target_id:** `SLK`
- **selection rule:** Exclude LCK. Drop target_ids containing '('. Require n_binders >= 10 at pKd >= 7.0. Pick binder count closest to LCK (17). Tie break: more binders, then target_id.
- **eligible kinases:** 30; chosen binders=17 (LCK had 17)
- **runners-up:** GAK (16), LOK (19), DDR2 (15), VEGFR2 (15), FLT1 (20)
- **sequence length:** 1235
- **nearest other kinase by AAC+DPC cosine:** `NEK1` (0.9702)
- **DAVIS labeled pairs held out:** 68
- **DAVIS binders (pKd >= 7.0):** 17
- **train pairs (all other target_ids):** 25704

Cold-target, not cold-ligand: the same 68 ligands still have labels on other kinases, so a ligand-only model can rank chemotypes that were potent elsewhere. The protein model never sees this target_id.

## Screening library (reused, not rebuilt)

- **n molecules:** 8068
- **library CSV:** `data/screen_library/library_v1.5.csv`
- **library SHA256:** `5721c043c87b744e4f288d31319e7743bf72e7a064de658708b0a16c61df0604`
- **source note:** v1.5 MoleculeNet HIV subsample plus DAVIS ligand spike-in. HIV labels unused. Not redownloaded.

## Recovery of known binders

Top 1% cutoff: rank <= 81 (k = max(1, ceil(0.01 * n)); same cutoff as EF@1%). Top 500: rank <= 500.
Fractions use known binders present in the library as the denominator.

| model | n binders in library | fraction top 1% | n top 1% | fraction top 500 | n top 500 | mean rank | random mean rank |
|---|---:|---:|---:|---:|---:|---:|---:|
| protein (ECFP+aac_dpc) | 17 | 0.4706 | 8 | 0.5882 | 10 | 1311.2353 | 4034.5000 |
| ligand-only | 17 | 0.4118 | 7 | 0.7059 | 12 | 469.1765 | 4034.5000 |

Random expected counts: top 1% ~ 0.1707 (fraction 0.0100); top 500 ~ 1.0535 (fraction 0.0620).
Mean rank vs random (random/model, higher is better than chance): protein 3.0769, ligand-only 8.5991.

**Does the protein model beat ligand-only on recovery?** mixed

yes if protein is strictly better on at least one of {top 1% fraction, top 500 fraction, mean rank} and worse on none; no if ligand-only is strictly better on at least one and protein on none; mixed if each wins at least one; tie if all equal. Lower mean rank is better.

Deltas (positive means protein better): top 1% fraction 0.0588, top 500 fraction -0.1176, mean-rank improvement -842.0588.

## Known binder detail (sorted by protein rank)

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

## Top-10 by protein score

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

## Artifacts

- Ranked CSV: `artifacts/screen_ranked_cold.csv`
- Models dir: `artifacts/models_v1.6`

Scores are predicted pKd from a model. They are not confirmed binders and not wet-lab hits.
