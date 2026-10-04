# Research Pack: Target → Ranked Small-Molecule Virtual Screening

**Owner:** Sofia Kuttner Lindelow · **Status:** research pack v1 · **Date:** 2026-10-02 (ET)  
**Companion files:** `DESIGN.md` (locked product), `DESIGN_TWEAKS.md` (memo)

---

## 0. Executive takeaway

Lock a **leakage-safe, laptop-tractable** pipeline: freeze a **kinase-focused affinity matrix** (DAVIS ± KIBA) plus a **small public screening library** (≈2k–20k SMILES), train an **ECFP + gradient-boosted** baseline and one stretch scorer (**Chemprop ligand-only** or **PSICHIC sequence+SMILES**), and report **Spearman ρ, enrichment@k, AUROC (binder), ECE**, plus optional **target−antitarget Δscore**. Do **not** chase SAVI-scale (≈1.75B) or Boltz-2/DiffDock as v1 deliverables—use them as optional week-3+ oracles on a shortlist.

Success claim (honest estimand): under **scaffold-aware and/or protein-holdout splits**, predicted ranks enrich true binders and scores are **calibrated enough to say when not to trust top-k**.

---

## 1. Problem definition & prior art

### 1.1 What virtual screening is (and is not)

Virtual screening (VS) prioritizes compounds from a library against a target so that experimental follow-up on the **top-k** recovers more actives than random. Classical framing: docking + scoring as enrichment tools, not perfect affinity oracles ([Kitchen et al., Nat. Rev. Drug Discov. 2004](https://www.nature.com/articles/nrd1549)).

| Paradigm | Input | Strength | Weakness for Sofia v1 |
|---|---|---|---|
| **Structure-based docking** | 3D protein + ligand | Pose + physics priors | Pose error, slow at library scale, scoring ≠ affinity |
| **Ligand-based** | Known actives → similarity / QSAR / fingerprints | Fast, laptop-friendly | Needs actives; weak on novel chemotypes |
| **Structure-based ML scorers** | 3D complex graphs / grids | Can beat classical SFs *on leaked benchmarks* | **Leakage** into CASF/core sets; GPU-heavy |
| **Sequence + SMILES ML** | UniProt sequence + SMILES | No crystal required; matches product input | Easy to overfit scaffolds / popular targets |

**Product alignment:** Sofia’s input is UniProt ID / sequence (± antitargets) and a **fixed SMILES library**—so **sequence+ligand ML** and **ligand-based baselines** fit v1 better than docking-first pipelines.

### 1.2 Landmark / relevant methods (with URLs)

| Work | What it is | Why it matters here | URL |
|---|---|---|---|
| **PSICHIC** (Koh et al., *Nat. Mach. Intell.* 2024) | Physicochemical GNN; affinity + functional effects from **sequence + ligand** only; interpretable interaction fingerprints; experimental A₁R agonist screen | Stretch scorer that matches “target sequence → rank library”; selectivity profiling story | https://www.nature.com/articles/s42256-024-00847-1 · code https://github.com/huankoh/PSICHIC · Zenodo https://zenodo.org/records/10901686 |
| **Chemprop** (Yang et al. JCIM 2019; Heid et al. JCIM 2024; Chemprop v2 Graff et al. JCIM 2025/26) | Directed MPNN for molecular property prediction; mature CLI/API; uncertainty tooling in package lineage | Best **ligand-only** stretch: multi-task affinity / ADMET side heads; laptop GPU optional | https://github.com/chemprop/chemprop · theory https://doi.org/10.1021/acs.jcim.9b00237 · package https://doi.org/10.1021/acs.jcim.3c01250 · docs https://chemprop.readthedocs.io/ |
| **Boltz-2** (Passaro / Corso / Wohlwend et al., bioRxiv 2025) | Co-folding + **affinity** head; claims FEP-like ranking on some benches at ≫1000× FEP speed | Optional **oracle / week-3+** rescoring of top-50–200; **not** full-library screen on laptop | https://www.biorxiv.org/content/10.1101/2025.06.14.659707v1 · code https://github.com/jwohlwend/boltz · PMC https://pmc.ncbi.nlm.nih.gov/articles/PMC12262699/ |
| **DiffDock** (Corso et al., ICLR 2023) | Diffusion docking over translation/rotation/torsion; pose + confidence | **Context only**—pose generation, not affinity ranking; useful later if structure-based shortlist | https://arxiv.org/abs/2210.01776 |
| **MoleculeNet / BBBP** (Wu et al., *Chem. Sci.* 2018) | Standard ML chem benchmarks; BBBP ≈2k compounds, **scaffold split**, AUROC | Template for library size, scaffold split habit, and “physiology” classification metrics | https://pubs.rsc.org/en/content/articlelanding/2018/sc/c7sc02664a · arXiv https://arxiv.org/abs/1703.00564 |
| **DeepDTA** (Öztürk et al., *Bioinformatics* 2018) | CNN on SMILES + protein sequence; DAVIS/KIBA | Classic sequence+SMILES baseline; shows kinase-matrix feasibility | https://pubmed.ncbi.nlm.nih.gov/30423097/ · arXiv https://arxiv.org/abs/1801.10193 |
| **Kitchen et al. 2004** | Docking/scoring review | Language for enrichment vs absolute affinity | https://www.nature.com/articles/nrd1549 |
| **Early-recognition metrics** (Truchon & Bayly, JCIM 2007) | BEDROC / why ROC can mislead early recovery | Use EF@k / BEDROC alongside AUROC | https://doi.org/10.1021/ci600423u |
| **Calibration** (Guo et al., ICML 2017) | ECE for neural nets | Secondary estimand: top-decile “binder” rates | https://proceedings.mlr.press/v70/guo17a.html |
| **TASC-VS** (Zou et al., *PLoS One* 2026) | Calibration-aware top-k VS decision framing | Modern framing aligned with Sofia’s “when not to trust” pitch | https://doi.org/10.1371/journal.pone.0356482 |
| **PDBbind CleanSplit / GEMS** (Graber et al., *Nat. Mach. Intell.* 2025) | Shows CASF↔PDBbind leakage; CleanSplit collapses inflated SOTA | Absolute weapon for Mila honesty narrative | https://www.nature.com/articles/s42256-025-01124-5 |
| **LP-PDBBind** (Li et al., arXiv 2023) | Similarity-controlled PDBbind split; BDB2020+ external set | Complementary leakage control | https://arxiv.org/abs/2308.09639 · https://github.com/THGLab/LP-PDBBind |

### 1.3 Positioning vs industry “billion-scale” marketing

SAVI-2020 enumerates **~1.75 billion** proposed products (~1.53B unique), Plus class alone ~1.09B ([Scientific Data 2020](https://www.nature.com/articles/s41597-020-00727-4); download https://cactus.nci.nih.gov/download/savi_download/). Newer **SAVI-Space** combinatorial encodings go further ([Scientific Data 2025](https://www.nature.com/articles/s41597-025-05384-z)). Sofia’s v1 library is **2k–50k**—own that as **method + calibration science**, not ultra-large docking theater.

---

## 2. Datasets: sizes, access, feasibility

### 2.1 Affinity / interaction training data

| Resource | Rough scale (public figures) | Access | Laptop fit | Notes |
|---|---|---|---|---|
| **BindingDB** | ~**3.25M** measurements, ~**1.45M** compounds, ~**11.5k** targets (site info as of fetch) | Free TSV/SDF downloads; CC licenses (ChEMBL-imported vs curator-curated differ) | Full dump heavy; **slice by target family** | https://www.bindingdb.org/rwd/bind/info.jsp · download https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp · NAR 2025 BindingDB-in-2024 https://doi.org/10.1093/nar/gkae1075 |
| **PDBbind v2020** | General ~**19k** complexes; refined **5,316**; CASF-2016 **285** | Registration historically at pdbbind.org.cn / PDBbind+; mirrors on HF | Structure files ~GBs; refined alone OK | Liu et al. Bioinformatics 2015 https://doi.org/10.1093/bioinformatics/btu626 · CleanSplit paper uses v2020 |
| **DAVIS** | **68** ligands × **442** kinases → **30,056** Kd pairs (DeepDTA packaging) | Bundled in DeepDTA / TDC / Zenodo PGDTA packs | **Excellent v1 freeze** | DeepDTA data notes; TDC DTI https://tdcommons.ai/multi_pred_tasks/dti/ |
| **KIBA** | Common ML pack: **2,111** drugs × **229** kinases → **~118k** scores; Tang et al. integrated larger matrix historically | Same portals | **Excellent stretch freeze** | Tang et al. JCIM 2014 https://doi.org/10.1021/ci400709d |
| **ChEMBL 35** | **~21.1M** activity records (release notes) | FTP dumps (Postgres/MySQL/SQLite) https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_35/ | Full DB heavy; use **filtered extracts** | Blog https://chembl.blogspot.com/2024/12/heres-nice-christmas-gift-chembl-35-is.html |
| **SAVI** | **1.75B** products | NCI CADD download (multi-TB uncompressed) | **Not feasible** for evening laptop full scan | Diversity subsets exist (~3M, ~16M) if needed later |

**TDC / HuggingFace convenience packs:** Therapeutics Data Commons DTI tasks (BindingDB, DAVIS, KIBA) with random / cold-drug / cold-target splits ([TDC](https://tdcommons.ai/multi_pred_tasks/dti/)). Bindwell/PLBA and Zenodo PGDTA bundles mirror KIBA/DAVIS/PDBbind tables for quick start.

### 2.2 Screening library (inference-time compound set)

| Option | Size | Role |
|---|---|---|
| MoleculeNet **BBBP** (~2,039) / Tox21 slice | Tiny | Smoke-test CLI + scaffold-split demo |
| ChEMBL drug-like random / diversity slice | 5k–20k | Realistic demo library |
| Enamine / ZINC diversity subsets | 10k–50k | Still laptop-OK for fingerprint/Chemprop |
| SAVI diversity (~3M SDF on NCI page) | Millions | Only after score is frozen; batch jobs |

### 2.3 Recommended **v1 data freeze** (actionable)

**Freeze name:** `sofia-vs-v1-2026-10`

1. **Train/eval affinity matrix:** DAVIS (primary) + optional KIBA (secondary). Convert Kd → pKd; document assay type.
2. **Protein identifiers:** UniProt accessions + amino-acid sequences (FASTA freeze file).
3. **Ligands:** Canonical SMILES + InChIKey; RDKit sanitize filter; drop failed parses.
4. **Splits (committed in git):**  
   - Scaffold split on ligands (Bemiscal / Murcko via RDKit) for within-target ranking tests.  
   - **Cold-protein** holdout (leave-family or leave-kinase-group) for “new target” claim.  
   - Optional **temporal** slice if BindingDB publication dates used later (BindingDB now exposes pub/curation dates for AI splits).
5. **Screening library freeze:** 5k–15k SMILES diversity set (ChEMBL or ZINC-like), versioned CSV + hash.
6. **Antitarget sketch (optional week 2):** 1–2 related kinases as antitargets for Δscore demos (PSICHIC Fig. 5 selectivity framing).
7. **Do not freeze for v1:** full BindingDB, full ChEMBL, PDBbind structures (unless week-3 Boltz-2 shortlist), SAVI.

Manifest discipline: `data/manifest.json` with URLs, download dates (ET), SHA256, row counts, license notes—NeuroAI-style reproducibility.

---

## 3. Splits & leakage: what papers get wrong

### 3.1 Split types Sofia should use

| Split | Protects against | Misses |
|---|---|---|
| **Random pair** | Nothing meaningful | Scaffold twins, same pocket, ligand memorizers |
| **Scaffold (ligand)** | Near-duplicate chemotypes | Same target family leakage; interaction clones |
| **Cold-protein / family holdout** | “Seen this kinase” | Similar ligands across targets |
| **Temporal** | Future assays | Popular scaffolds re-assayed later |
| **Structure-similarity (CleanSplit / LP)** | Pocket+ligand+pose clones | Needs 3D; overkill for sequence-only v1 but cite for honesty |

### 3.2 Documented failure modes (cite these)

1. **PDBbind ↔ CASF leakage:** Graber et al. (*Nat. Mach. Intell.* 2025) find ~**49%** of CASF complexes have near-identical training mates; simple “copy affinity of top-5 similar train complexes” rivals published DL SFs. Retraining Pafnucy/GenScore on **CleanSplit** **collapses** scores. Ligand-only models can look strong when leakage exists ([article](https://www.nature.com/articles/s42256-025-01124-5)).
2. **LP-PDBBind (Li et al. 2023):** Original general/core overlap; time splits insufficient because old targets get new ligands. Similarity-controlled split + **BDB2020+** external set ([arXiv](https://arxiv.org/abs/2308.09639)).
3. **Ligand memorization in GNNs:** Models can ignore protein and still “score” well on contaminated benches (discussed in CleanSplit; Mastropietro et al. *Nat. Mach. Intell.* 2023 cited therein).
4. **IC50 mixed with Kd/Ki:** Landrum & Riniker (JCIM 2024) warn combining IC50 with Ki/Kd injects noise—flag assay type in freeze ([CleanSplit cites](https://www.nature.com/articles/s42256-025-01124-5)).
5. **DUD-E / LIT-PCBA hidden bias:** Chen et al. *PLoS One* 2019 (DUD-E bias); Huang et al. 2025 leakage audit of LIT-PCBA (arXiv:2507.21404)—don’t claim SOTA on decoy sets without cleaning.
6. **MoleculeNet lesson:** BBBP recommends **scaffold** split and AUROC—copy the habit, not only the dataset ([Wu et al. 2018](https://pubs.rsc.org/en/content/articlelanding/2018/sc/c7sc02664a)).

### 3.3 Practical split policy for v1

- **Primary report:** ligand **scaffold split** within DAVIS (Spearman / EF on held-out scaffolds).  
- **Stress test:** **cold-protein** (hold out entire kinases / groups).  
- **Ablation honesty:** train a **ligand-only** model; if cold-protein barely drops, you’re memorizing ligands—say so.  
- **Do not** lead with random-split AUROC.

---

## 4. Metrics (honest estimands)

### 4.1 Ranking / enrichment

| Metric | Definition / use | Primary? |
|---|---|---|
| **Spearman ρ** | Rank correlation of score vs continuous affinity (pKd / KIBA score) | **Yes** (continuous labels) |
| **EF@k** | (actives in top k%) / (overall active rate); classic early enrichment | **Yes** for binary binder cut (e.g. pKd ≥ 7) |
| **BEDROC** | Early-weighted ROC (Truchon & Bayly 2007) | Nice-to-have |
| **AUROC** | Binder vs non-binder discrimination | Secondary; **not** success alone |
| **Precision@k / hits in top-20** | Fixed experimental budget framing (TASC-VS) | Product-facing |

### 4.2 Calibration

| Metric | Use |
|---|---|
| **ECE** (Expected Calibration Error) | Bin predicted binder probabilities vs empirical frequencies ([Guo et al. 2017](https://proceedings.mlr.press/v70/guo17a.html)) |
| **Reliability diagram** | Visual for write-up / startup pitch |
| **Brier / NLL** | Proper scoring (optional) |

**Claim language:** “When the model puts a pair in the top score decile, empirical binder rate on scaffold-held-out pairs is X% (vs Y% base rate).” Not “our AUROC is 0.92.”

### 4.3 Selectivity (optional secondary)

- Score both **target** and **antitarget** sequences for the same ligand.  
- Report **Δscore = score(target) − score(antitarget)** vs known selective ligands (PSICHIC selectivity profiling precedent: https://www.nature.com/articles/s42256-024-00847-1 Fig. 5).  
- Metric: Spearman of Δscore vs experimental selectivity index, or enrichment of selective ligands in top Δscore.

### 4.4 Uncertainty

Chemprop lineage includes uncertainty / calibration workflows; PSICHIC importance scores are interpretability, not calibrated probs—don’t conflate.

---

## 5. Baseline vs stretch models + compute estimates

Assumptions: consumer laptop (8–16 CPU threads, ≤16 GB RAM, optional single consumer GPU 8–12 GB) · library ≤20k · evenings over ~3 weeks.

| Tier | Model | Train compute (order-of-magnitude) | Screen 10k compounds | Fits week 1–2? |
|---|---|---|---|---|
| **Baseline** | ECFP4/6 (2048–4096) + LightGBM/XGBoost/logistic; protein as target-ID embedding or ESM mean pooled (optional) | Minutes–1 h CPU on DAVIS | Seconds–minutes | **Yes — must ship** |
| **Stretch A** | **Chemprop** D-MPNN on ligand; multi-task per-target heads or global + target features | 0.5–4 h GPU or overnight CPU on DAVIS/KIBA | Minutes on GPU; tens of min CPU | **Yes** |
| **Stretch B** | **PSICHIC** (pretrained weights + optional fine-tune) | Inference-heavy; fine-tune overnight if attempted | Minutes–hours depending on batching | **Maybe** (prefer pretrained inference first) |
| **Oracle (not v1 core)** | **Boltz-2** affinity on top-50–200 | Hours–days GPU per complex batch; not for 10k | Rescore shortlist only | Week 3+ optional |
| **Context only** | DiffDock / Vina | Pose, not affinity | Expensive at scale | Out of scope for ranking claim |

**Head-to-head protocol:** identical splits, identical library, identical metrics; baseline must be reported even if stretch wins.

**Dependency stack suggestion:** Python + `uv`; RDKit; scikit-learn / LightGBM; optional `chemprop`; optional clone of PSICHIC; `reproduce.sh` + frozen manifests.

---

## 6. Dual-fit narrative (startup + Mila)

| Audience | Line |
|---|---|
| **Health-AI startup** | “Screening API: ranked candidates with **measured enrichment** and **calibration**—we know when not to trust top-k.” |
| **Mila / research** | “Oracle scores under **leakage-safe** splits; Goodhart risk of optimizing random-split AUROC; CleanSplit/LP-PDBBind as cautionary priors.” |

Non-goals remain: wet lab, de novo generation, Bittensor/SN68 mining, UI polish, billion-scale SAVI marketing.

---

## 7. Risks & mitigations

| Risk | Mitigation |
|---|---|
| Noisy / biased labels | Prefer Kd-heavy DAVIS; flag IC50; report CIs across targets |
| Leakage inflation | Scaffold + cold-protein; ligand-only ablation |
| Small library looks weak vs SAVI ads | Explicit scope statement + diversity metrics of library |
| Stretch model eats the month | Cap: baseline week 1 mandatory; stretch time-boxed |
| Boltz-2 hype | Cite independent caution (e.g. ranking limits in follow-on evals); use only as shortlist oracle |

---

## 8. References (selected, clickable)

1. Kitchen DB et al. *Nat. Rev. Drug Discov.* **3**, 935–949 (2004). https://www.nature.com/articles/nrd1549  
2. Wu Z et al. MoleculeNet. *Chem. Sci.* **9**, 513–530 (2018). https://doi.org/10.1039/C7SC02664A  
3. Öztürk H et al. DeepDTA. *Bioinformatics* **34**, i821–i829 (2018). https://doi.org/10.1093/bioinformatics/bty593  
4. Tang J et al. KIBA. *J. Chem. Inf. Model.* **54**, 735–743 (2014). https://doi.org/10.1021/ci400709d  
5. Yang K et al. Learned molecular representations. *J. Chem. Inf. Model.* **59**, 3370–3388 (2019). https://doi.org/10.1021/acs.jcim.9b00237  
6. Heid E et al. Chemprop package. *J. Chem. Inf. Model.* **64**, 9–17 (2024). https://doi.org/10.1021/acs.jcim.3c01250  
7. Koh HY et al. PSICHIC. *Nat. Mach. Intell.* **6**, 673–687 (2024). https://doi.org/10.1038/s42256-024-00847-1  
8. Corso G et al. DiffDock. ICLR 2023. https://arxiv.org/abs/2210.01776  
9. Passaro S et al. Boltz-2. bioRxiv (2025). https://doi.org/10.1101/2025.06.14.659707  
10. Graber D et al. PDBbind CleanSplit. *Nat. Mach. Intell.* (2025). https://doi.org/10.1038/s42256-025-01124-5  
11. Li J et al. LP-PDBBind. arXiv:2308.09639. https://arxiv.org/abs/2308.09639  
12. Liu T et al. BindingDB. *Nucleic Acids Res.* (ongoing releases). https://www.bindingdb.org/  
13. Patel H et al. SAVI. *Sci. Data* **7**, 364 (2020). https://doi.org/10.1038/s41597-020-00727-4  
14. Truchon J-F & Bayly CI. BEDROC / early recognition. *J. Chem. Inf. Model.* **47**, 488–508 (2007). https://doi.org/10.1021/ci600423u  
15. Guo C et al. On calibration of modern neural networks. ICML 2017. https://proceedings.mlr.press/v70/guo17a.html  
16. Zou L et al. TASC-VS calibration-aware top-k. *PLoS One* (2026). https://doi.org/10.1371/journal.pone.0356482  
17. Su M et al. CASF-2016. *J. Chem. Inf. Model.* **59**, 895–913 (2019). https://doi.org/10.1021/acs.jcim.8b00545  

---

## 9. Pointers to next artifacts

- Concrete design changes → `DESIGN_TWEAKS.md`  
- Locked product shape → `DESIGN.md` (unchanged product lock; tweaks are advisory)
