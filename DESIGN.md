# Design: Target → Ranked Candidate Virtual Screening
**Owner:** Sofia · **Status:** locked · **Date:** 2026-10-02 (ET)

## One-liner
Train a protein–ligand **activity predictor**, then wrap it as a **screening tool**: given a receptor, rank a molecule library into probable candidates.

## Architecture (predictor first, screening as product)
```
known (protein, molecule, affinity) data
        ↓ train / validate
   activity predictor  ← scientific core (honest splits + metrics)
        ↓ apply to library
   screen(target, library) → ranked candidates  ← product wrapper
```

| Layer | What it is | How we know it works |
|---|---|---|
| **Predictor** | Score for one (target, molecule) pair | Spearman / EF@k / ECE / AUROC on held-out **experimental** affinities (DAVIS ± KIBA) |
| **Screener** | Sort a library by that score for a new target | Same model; demo CLI `screen → ranked.csv` |

Classification/prediction metrics validate the brain; the ranked list is how you *find* candidates.

## Product shape
| | |
|---|---|
| **Input** | Target protein (UniProt ID or amino-acid sequence); optional antitarget list |
| **Library** | Fixed public compound set (v1: ~10k ±5k SMILES) — *not* billion-scale SAVI |
| **Output** | Ranked table: molecule ID, SMILES, score, optional selectivity vs antitarget, uncertainty/calibration note |
| **Non-goals (v1)** | Wet-lab synthesis, de novo molecule generation, Bittensor mining, UI polish |

## Scientific estimand (what we claim)
Under **scaffold-aware and cold-protein** holdout on public protein–ligand affinity data:
1. **Primary:** Ranking quality = Spearman ρ / enrichment@k of predicted scores vs experimental affinity (or binary binder labels).
2. **Secondary:** Calibration — when score says “top 10%,” what fraction are true binders on held-out pairs? (AUROC secondary.)
3. **Control:** Ligand-only ablation (no protein features) — detect memorization.
4. **Optional:** Selectivity — Δscore(target − antitarget).

Success ≠ “high AUROC on a random split.”

## External checks (not primary labels)
- **NOVA / SN68 public submissions:** optional stress test — for a published weekly target, do our ranks put NOVA-submitted molecules high? Useful narrative; **not gold truth** (those molecules were selected by an AI oracle, not wet-lab assays).
- **Primary ground truth:** experimental affinity matrices (DAVIS ± KIBA; BindingDB later).

## Stack (v1, laptop-tractable)
- **Train/eval freeze:** DAVIS pKd (`sofia-vs-v1-2026-10`); KIBA optional secondary.
- **Screen library:** ~10k diversity SMILES, versioned + hashed.
- **Baseline:** ECFP fingerprints + gradient-boosted / logistic model.
- **Stretch:** Chemprop first; optional PSICHIC inference; Boltz-2 ≤100 shortlist only if GPU free.
- **Splits:** Murcko/Bemiscal scaffold **and** cold-protein/family.
- **Code:** Python + `uv`; `reproduce.sh`; frozen manifests (NeuroAI habits).

## 3-week plan
| Week | Deliverable |
|---|---|
| 1 | Data freeze + splits + baseline predictor; ρ / enrichment on held-out pairs |
| 2 | Chemprop head-to-head; calibration; ligand-only ablation; optional antitarget sketch |
| 3 | CLI screener wrapper; optional NOVA-submission concordance check; methods write-up + failure cases |

## Dual-fit framing
- **Health-AI startup:** screening tool with a measured, calibrated predictor underneath.
- **Mila / research:** leakage-safe DTI evaluation; Goodhart risk of trusting proxy oracles (incl. NOVA-style scores).

## Risks
- Affinity labels are noisy and biased toward studied targets.
- Small libraries look weak next to NOVA marketing — own method + calibration first.
- Treating NOVA winners as “true binders” would bake in oracle bias — keep them as secondary concordance only.

## Decision log
- 2026-10-02: Locked product = target → ranked candidates. Venue = off-chain science, not SN68 mining.
- 2026-10-02: Research pack → DAVIS freeze, dual splits, ligand-only ablation, Chemprop-first stretch.
- 2026-10-02: Clarified architecture = **predictor first, screening tool as product wrapper**.
- 2026-10-02: NOVA public hits = optional external check, not primary ground truth.

## Same-day MVP (2026-10-02)
**Goal:** done today, not 3 weeks.
- MUST: DAVIS download + freeze · ECFP+GBM baseline · scaffold + cold-protein metrics (Spearman, EF@k, AUROC) · ligand-only ablation · `screen` CLI → ranked.csv · short README
- NICE: Chemprop · calibration plot · NOVA concordance
- CUT: Boltz-2 · full KIBA · selectivity pillar · UI

## Verifiability (hard requirement — 2026-10-02)
All steps, data, calculations, and analysis must be independently checkable:
- Frozen inputs with source URL, download date, file hashes (SHA256) in `data/MANIFEST.json`
- Locked random seeds; package versions pinned (`uv.lock` / requirements)
- Split definitions written to disk (indices or IDs), not only described in prose
- Metrics recomputed by `reproduce.sh` → bit-stable (or documented float tolerance)
- Model artifacts + training config logged; inference CLI uses those exact weights
- No silent data edits; any filter/threshold (e.g. pKd≥7) documented in config YAML
- Public companion must let a stranger re-run and match the reported table
