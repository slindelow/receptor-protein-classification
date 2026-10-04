# Sofia VS - results summary (v1 to v1.5)

**Data freeze:** `sofia-vs-v1-2026-10` (TDC DAVIS, pKd = 9 − log10(Kd nM), binder = pKd ≥ 7).  
**KIBA exception:** DeepDTA mirror (TDC Dataverse 403); label is KIBA score, not pKd; binder ≥ 12.1.  
**Box only.** Not published. Numbers copied from exported metrics files (2026-10-02, America/Toronto).

## What each step asked

| Step | Question | Answer on this freeze |
|---|---|---|
| v1 | ECFP+HistGBM+AAC vs ligand-only | Scaffold: protein helps (ρ 0.21 vs ~0). Cold-protein: AAC ≈ ligand-only (ρ ~0.50). |
| v1.1 | Stronger protein features on cold-protein? | **Yes.** aac_dpc ρ=0.550 (Δρ≈+0.05) and EF@1% 10.25 vs ligand-only 8.45. ESM2 similar. |
| v1.2 | Chemprop D-MPNN + aac_dpc on cold-protein? | **No.** Chemprop+aac_dpc ρ=0.490 vs Chemprop ligand-only 0.500. HistGBM+aac_dpc still better. |
| cal | Are scores calibrated as P(binder)? | Isotonic cuts ECE sharply (cold HistGBM ~0.08–0.10 → ~0.01). Ranking unchanged. |
| v1.3 | Does aac_dpc cold win transfer to KIBA? | **Yes** on Spearman: Δρ=+0.061 (0.542 vs 0.481). EF@1% flat. |
| v1.4 | Can the scorer rank selectivity on one kinase pair? | **Sketch yes** on SRC/LCK: Spearman(Δ)=0.369, EF@10%=3.24. Ligand-only Δ is constant (null). |

## DAVIS ranking (Spearman ρ / EF@1%)

| Split | Model | ρ | EF@1% |
|---|---|---:|---:|
| scaffold | HistGBM ligand-only | −0.011 | 0.00 |
| scaffold | HistGBM+AAC (v1) | 0.212 | 6.51 |
| scaffold | HistGBM+aac_dpc | 0.240 | 5.21 |
| scaffold | HistGBM+ESM2 | 0.257 | 5.86 |
| scaffold | Chemprop ligand-only | 0.262 | 0.33 |
| scaffold | Chemprop+aac_dpc | 0.318 | 7.48 |
| cold_protein | HistGBM ligand-only | 0.500 | 8.45 |
| cold_protein | HistGBM+AAC | 0.498 | 8.71 |
| cold_protein | HistGBM+aac_dpc | **0.550** | **10.25** |
| cold_protein | HistGBM+ESM2 | 0.548 | 10.25 |
| cold_protein | Chemprop ligand-only | 0.500 | 8.71 |
| cold_protein | Chemprop+aac_dpc | 0.490 | 8.45 |

Preferred deployed scorer for later steps: **HistGBM + aac_dpc**.

## KIBA (v1.3, score not pKd)

| Split | Model | ρ | EF@1% |
|---|---|---:|---:|
| scaffold | HistGBM ligand-only | 0.339 | 1.63 |
| scaffold | HistGBM+aac_dpc | 0.558 | 4.96 |
| cold_protein | HistGBM ligand-only | 0.481 | 3.90 |
| cold_protein | HistGBM+aac_dpc | 0.542 | 3.89 |

## Selectivity (v1.4, SRC / LCK)

Cold-pair holdout (train excludes both). Δ = score(SRC) − score(LCK). Selective = |ΔpKd| ≥ 1 (12/68 ligands, all LCK-preferring).

| Model | Spearman(Δ) | EF@5% | EF@10% | EF@20% |
|---|---:|---:|---:|---:|
| HistGBM+aac_dpc | 0.369 | 2.83 | 3.24 | 2.43 |
| HistGBM ligand-only | n/a (Δ constant) | 0.00 | 0.00 | 0.00 |

## Calibration (ECE, 10 bins)

Cold-protein HistGBM ECE uncalibrated → isotonic-on-val: ligand-only 0.104 → 0.010; aac_dpc 0.085 → 0.014; ESM2 0.094 → 0.015. Chemprop used a 30% test holdout (not val); ECE after isotonic ~0.005–0.007. Do not treat holdout ECE as an unbiased test.

## Limits (do not overclaim)

- DAVIS has only 68 ligands; cold-protein ρ ~0.50 is largely ligand-driven even when protein features add a small lift.
- v1.4 is one pair, n=68, sketch thresholds. Not a selectivity product.
- No wet-lab, NOVA, or prospective screen claim.
- v1.4 was last classification/selectivity P2; v1.5 adds the screening product.

## v1.5 screening product (`sofia-vs-v1.5-screen`)

**Question:** Can HistGBM+aac_dpc rank a ~8k public library for one DAVIS kinase?  
**Answer:** Yes as a model-score product. Target **LCK**. Library n=8068 (MoleculeNet HIV diverse subsample + DAVIS ligand spike-in). Ranked CSV: `artifacts/screen_ranked.csv`. Weights: scaffold train+val `artifacts/models_v1.5/hgb_screen_aac_dpc.joblib`.

Known LCK binders recovered: 15/17 in top-500 (random expected ~1.05; ~14x enrichment). Mean binder rank 378 vs random ~4035. These are predicted pKd scores, not wet-lab hits. Recovery uses labels that partly overlap train (scaffold train+val), so treat as sanity check only.

## v1.6 cold-target screen (`sofia-vs-v1.6-cold-screen`)

**Question:** If the scored kinase is fully held out, does HistGBM+aac_dpc still recover its known binders better than ligand-only?  
**Answer:** No. Recovery stays above random, but it is mixed and ligand-only wins the broader ranks. Target **SLK** (17 binders, pKd >= 7). All 68 SLK pairs held out. n_train=25704 on the remaining DAVIS pairs. Library reused from v1.5 (n=8068), not redownloaded.

| Model | top 1% (k=81) | top 500 | mean rank | vs random mean (4034.5) |
|---|---:|---:|---:|---:|
| HistGBM+aac_dpc | 8/17 (0.4706) | 10/17 (0.5882) | 1311.2 | ~3.1x |
| ligand-only | 7/17 (0.4118) | 12/17 (0.7059) | 469.2 | ~8.6x |

Protein model beats ligand-only on recovery: **no**. Mixed: protein is higher only in the top 1%; ligand-only is better in the top 500 and has the better mean rank. Same ligands remain labeled on other kinases, so this is a cold target, not a cold ligand. Ranked CSV: `artifacts/screen_ranked_cold.csv`. These are predicted pKd scores, not wet-lab hits.
