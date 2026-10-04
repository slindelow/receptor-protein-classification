# Design Tweaks Memo — Virtual Screening v1

**From:** research pass on `DESIGN.md` (2026-10-02 ET)  
**To:** Sofia · first-month decisions  
**Companion:** `RESEARCH_PACK.md`

This memo does **not** unlock the product (still: target → ranked candidates, off-chain, not Bittensor). It proposes concrete upgrades so the locked estimands (scaffold holdout, calibration, enrichment) are actually deliverable in ~3 evening-weeks.

---

## A. Gap check vs locked DESIGN.md

| DESIGN says | Research finding | Gap |
|---|---|---|
| BindingDB or PDBbind + MoleculeNet-sized library | Full BindingDB (~3.25M) / PDBbind structures are heavy; DAVIS/KIBA are the proven laptop DTI packs | Freeze **DAVIS (±KIBA)** first; BindingDB as optional expand |
| ECFP + GBM baseline; PSICHIC *or* Chemprop stretch | Both stretches are real; Chemprop is lower ops risk; PSICHIC better story match | Sequence: baseline → Chemprop → optional PSICHIC inference |
| Scaffold splits | Scaffold alone insufficient (CleanSplit / LP-PDBBind) | Add **cold-protein** + ligand-only ablation |
| Boltz-2 “optional week-3+” | Correct—affinity oracle on shortlist only, not library | Keep optional; time-box hard |
| Selectivity optional | PSICHIC has selectivity precedent; cheap if sequence scorer exists | Keep week-2 sketch, not blocker |
| Success ≠ random AUROC | Literature strongly agrees (CASF leakage, DUD-E bias) | Make **EF@k + Spearman + ECE** the headline table |

---

## B. Concrete tweaks (8)

### Tweak 1 — Freeze `sofia-vs-v1` = DAVIS primary (+ KIBA optional), not BindingDB/PDBbind first

| | |
|---|---|
| **Change** | Training labels = DAVIS pKd matrix (68×442 / 30,056 pairs). Optional second freeze: KIBA (~118k). Screening library = separate 5k–15k SMILES CSV. Defer full BindingDB/PDBbind. |
| **Why** | DeepDTA/TDC ecosystem; kinase biology is coherent; laptop-sized; avoids BindingDB cleaning hell and PDBbind structure I/O in week 1. |
| **Cost** | Low (hours to download + sanitize). |
| **Dual-fit** | Startup: clear “kinase screening demo.” Mila: standard DTI benchmark with colder splits than random. |
| **Recommend** | **YES** for first month. |

### Tweak 2 — Dual primary splits: scaffold **and** cold-protein (report both)

| | |
|---|---|
| **Change** | Every model card reports (i) ligand Murcko/Bemiscal scaffold holdout metrics and (ii) leave-kinase / leave-protein-family metrics. Ban leading with random-split AUROC. |
| **Why** | CleanSplit/LP-PDBBind show random & contaminated benches inflate DL; scaffold alone still leaks proteins. |
| **Cost** | Low–medium (split code + 2× eval loops). |
| **Dual-fit** | **Highest Mila signal** in the project; startup gets “we stress-tested generalization.” |
| **Recommend** | **YES**. |

### Tweak 3 — Mandatory ligand-only ablation

| | |
|---|---|
| **Change** | Train ECFP model **without** protein features on same splits. If cold-protein gap ≈ 0, disclose “ligand memorization.” |
| **Why** | CleanSplit shows ligand-only can fake strong CASF-like numbers under leakage; PSICHIC/sequence claims need this control. |
| **Cost** | Low (one extra baseline). |
| **Dual-fit** | Mila honesty; startup trust (“we know failure modes”). |
| **Recommend** | **YES**. |

### Tweak 4 — Metric dashboard: Spearman + EF@1%/5% + AUROC + ECE (headline ≠ AUROC)

| | |
|---|---|
| **Change** | Fixed eval script outputs Spearman ρ (continuous), EF@k & Precision@20 (binary binder threshold, document cut e.g. pKd≥7), AUROC, ECE + reliability plot. Cite Truchon/Bayly + Guo. |
| **Why** | Matches locked estimands; TASC-VS shows decision value is top-k + calibration, not global ROC alone. |
| **Cost** | Low (1–2 evenings). |
| **Dual-fit** | Startup pitch language; Mila metric hygiene. |
| **Recommend** | **YES**. |

### Tweak 5 — Stretch model order: Chemprop before PSICHIC fine-tune

| | |
|---|---|
| **Change** | Week-2 stretch = Chemprop (ligand MPNN ± target ID / simple protein fingerprint). Use **pretrained PSICHIC weights for inference/demo** if time; fine-tune PSICHIC only if Chemprop ships early. |
| **Why** | Chemprop is battle-tested, documented, uncertainty-friendly, evening-trainable. PSICHIC is the glamorous sequence+SMILES story but heavier and easier to burn the week on. |
| **Cost** | Chemprop: medium. PSICHIC fine-tune: high. Inference-only PSICHIC: medium. |
| **Dual-fit** | Chemprop = credible ML engineering; PSICHIC = *Nat. Mach. Intell.* citation gravity. |
| **Recommend** | **YES** (Chemprop first); **MAYBE** (PSICHIC fine-tune in month 1). |

### Tweak 6 — Library size lock: 10k ±5k diversity set; explicitly not SAVI

| | |
|---|---|
| **Change** | Fix screening library at ~10k (range 5k–15k) diversity SMILES; publish scaffold count / Murcko uniqueness. One-liner in README: “Not SAVI-scale (~1.75B); method & calibration first.” |
| **Why** | DESIGN already says 2k–50k; pinning avoids scope creep and NOVA-comparison anxiety. |
| **Cost** | Low. |
| **Dual-fit** | Startup: honest scope. Mila: avoids Goodhart on library size. |
| **Recommend** | **YES**. |

### Tweak 7 — Boltz-2 as shortlist oracle only (cap ≤100 compounds)

| | |
|---|---|
| **Change** | If GPU available in week 3: rescoring **top-50–100** baseline/Chemprop hits with Boltz-2 affinity; report rank agreement vs wet-label proxies if any. Never claim full-library Boltz-2 screen. |
| **Why** | Aligns with Boltz-2 positioning (affinity + structure; still costly per complex) and independent caveats on ranking precision. DiffDock stays context-only. |
| **Cost** | High GPU time; high engineering. |
| **Dual-fit** | Startup “optional deep check”; Mila “foundation-model vs fingerprint concordance.” |
| **Recommend** | **MAYBE** (only if weeks 1–2 green); else backlog. |

### Tweak 8 — Selectivity = one antitarget pair demo, not a product pillar

| | |
|---|---|
| **Change** | Pick one kinase pair (e.g. related off-target) for Δscore table on 50–200 ligands with known selectivity if available; else synthetic demo from DAVIS multi-kinase labels. |
| **Why** | DESIGN optional; PSICHIC Fig. 5 shows narrative value; full antitarget product is out of scope. |
| **Cost** | Low–medium. |
| **Dual-fit** | Nice startup slide; light Mila add-on. |
| **Recommend** | **MAYBE** (week 2 if ahead); not a month-1 blocker. |

---

## C. Priority stack for Sofia’s first month

| Priority | Tweak | Verdict |
|---|---|---|
| P0 | 1 Data freeze DAVIS±KIBA | YES |
| P0 | 2 Scaffold + cold-protein | YES |
| P0 | 3 Ligand-only ablation | YES |
| P0 | 4 Metric dashboard | YES |
| P1 | 5 Chemprop stretch (PSICHIC infer optional) | YES / MAYBE fine-tune |
| P1 | 6 Library ~10k lock | YES |
| P2 | 8 Selectivity sketch | MAYBE |
| P2 | 7 Boltz-2 shortlist | MAYBE |

---

## D. Recommended revised one-pager outline (bullets)

Incorporate best tweaks into an updated DESIGN one-pager (product lock unchanged):

1. **One-liner** — Target sequence → ranked SMILES with defendable scores (enrichment + calibration under hard splits).  
2. **Input / library / output** — UniProt or FASTA; optional antitarget; frozen ~10k SMILES library; CSV with score, optional Δscore, calibration note.  
3. **Non-goals** — Wet lab, de novo gen, SAVI-scale, Bittensor mining, UI polish, full Boltz-2 library screen.  
4. **Estimands** — Primary: Spearman ρ + EF@k on scaffold holdout; Stress: cold-protein; Secondary: ECE; Optional: selectivity Δscore.  
5. **Data freeze `sofia-vs-v1`** — DAVIS (±KIBA); manifests + hashes; binary binder cut documented.  
6. **Models** — Baseline ECFP+GBM; ablation ligand-only; stretch Chemprop; optional PSICHIC pretrained inference; optional Boltz-2 top-100 oracle.  
7. **Splits** — Scaffold + cold-protein committed in repo; no random-split victory laps.  
8. **3-week plan** — W1 freeze+baseline+metrics; W2 Chemprop head-to-head+calibration(+selectivity sketch); W3 CLI `screen` + methods note + failure cases (+ optional Boltz-2).  
9. **Dual-fit** — Startup: calibrated screening API; Mila: leakage-aware evaluation (cite CleanSplit/LP-PDBBind).  
10. **Risks** — Label noise; small-library optics; stretch-model time sink; foundation-model compute.

---

## E. Open decisions (Sofia must choose)

See chat summary; blockers are: **primary matrix (DAVIS vs DAVIS+KIBA)**, **stretch (Chemprop vs PSICHIC priority)**, **binder threshold**, **antitarget yes/no**, **Boltz-2 attempt yes/no**.

---

*End of memo.*
