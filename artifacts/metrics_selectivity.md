# Selectivity metrics — sofia-vs-v1.4-selectivity

**Timestamp (America/Toronto):** 2026-10-02 17:43 EDT  
**Data freeze:** `sofia-vs-v1-2026-10` (MANIFEST verified)  
**Experiment freeze:** `sofia-vs-v1.4-selectivity`

## Kinase pair

| Field | Value |
|---|---|
| Pair | **SRC / LCK** |
| Why | Src-family pair (SRC/LCK); related+distinguishable sequences; dual-measured ligands; non-trivial selective set |
| Best-window sequence identity | 0.4008 |
| AAC+DPC cosine | 0.9621 |
| Dual-measured ligands | 68 |
| Selective (\|ΔpKd\| ≥ 1.0) | 12 |
| Protocol | cold-pair holdout (train excludes SRC and LCK) |
| Scorer | HistGBM + ECFP + **aac_dpc** (v1.1 cold winner) |

Δ := score(SRC) − score(LCK)  (same for true pKd).

## Results

| Model | Spearman(Δpred, Δtrue) | mean\|Δpred\| | EF@5% | EF@10% | EF@20% |
|---|---:|---:|---:|---:|---:|
| HistGBM+aac_dpc | 0.3694 | 0.1938 | 2.83 | 3.24 | 2.43 |
| HistGBM ligand-only | n/a (constant Δ) | 0.0000 | 0.00 | 0.00 | 0.00 |

EF is enrichment of ligands with |ΔpKd| ≥ 1 ranked by |Δpred| (n_selective=12 of 68). All 12 selective ligands prefer LCK (ΔpKd = pKd_SRC − pKd_LCK ≤ −1; none prefer SRC by ≥1), so LCK-directed EF equals the absolute EF above. SRC-directed EF is undefined (zero positives).

Ligand-only Δpred is exactly constant (protein ignored) → cannot rank selectivity. Honest control null.

## Success (sketch)

- Primary gate (ρ>0.1 and EF@10/20%>1.2): **YES**
- Soft signal: **YES**
- Honest null is an acceptable outcome for this P2 sketch.

## Takeaway

See `artifacts/DECISIONS.md` / `LAB_NOTEBOOK.md` (v1.4 section).
