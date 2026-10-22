# Publish prep

Date: 2026-10-03 (America/Toronto). This file records local prep only. Nothing was committed, pushed, or given a remote.

Intended public repo, later: https://github.com/slindelow/receptor-protein-classification

## Ready

- `RESEARCH_REPORT.md` no longer says the primary-split training and validation counts were not exported. The sentence now quotes the exported counts.
- Scaffold primary split: n_train 17813, n_val 2653, n_test 5306 (302 binders). Cold-protein primary split: n_train 18020, n_val 2584, n_test 5168 (388 binders). Source: `artifacts/metrics_v1.2.json` `train_info` (all three counts). `artifacts/metrics.json` and `artifacts/metrics_v1.1.json` have the same n_train and n_test and do not have n_val. The `split` column of `data/splits/scaffold_split.csv` and `data/splits/cold_protein_split.csv` matches those integers. Total pairs 25772.
- `README.md` points at `RESEARCH_REPORT.md` and at commands that exist: `reproduce.sh`, `reproduce_v1.1.sh`, `reproduce_v1.2.sh`, `reproduce_v1.4.sh`, `scripts/score_histgbm_cold_protein.py`, and `scripts/uncertainty_bootstrap.py`. It also names the Python entry points that exist without a shell wrapper: `scripts/run_kiba_v1.3.py`, `scripts/run_calibration.py`, `scripts/run_screen_v1.5.py`, and `scripts/run_screen_v1.6.py`.
- KIBA remains a point estimate. The report and README still say those per-pair scores were not saved. No new experiment was run.
- The word "significant" is not used as a statistical claim in `RESEARCH_REPORT.md` or `README.md`. No sentence was changed for that reason.
- No git repository exists under this directory (`git status` fails with "not a git repository"). `git init` was not run. There are no tracked files, so no tracked file is over 100 MB and no tracked file contains the forbidden names.

## Excluded, and why

Listed in `.gitignore` so a later commit can leave them out. Files were not deleted.

| Path | Why |
|---|---|
| `.venv/` | About 1.4 GB. Local environment. |
| `data/raw/*.tab`, `data/raw/davis_pairs.csv` | Raw DAVIS tables. Already ignored. |
| `data/raw/delaney_esol.csv` | Already ignored. Small. Not a DAVIS table. |
| `data/cache/` | ESM cache. Already ignored. |
| `CODEX_MINI_HANDOFF.md` | Names a private host user. Must stay unpublished. Already ignored. |
| `data/splits_kiba/*.csv` | Over GitHub's 100 MB file limit. Already ignored. |
| `data/raw/kiba_pairs.csv` | 101,906,920 bytes. Added to `.gitignore` in this prep. |
| `artifacts/models_v1.1/` | Already ignored. About 11 MB, under the file limit. The cold-protein score script needs these joblibs, so a clone will not be able to run that script unless the weights are supplied another way. |
| `artifacts/*.log`, `artifacts/models_v1.2/*.log` | Absolute machine paths. Added in this prep. |
| `artifacts/WRITING_BRIEF.md` | Personal email addresses, Drive links, and an unrelated writing sample. Added in this prep. |
| `LAB_NOTEBOOK.md` | Absolute machine path, and the notebook says it is not published. Added in this prep. |

Model directories other than `artifacts/models_v1.1/` are about 2 MB to 18 MB. They are not ignored. `artifacts/models_v1.2/` holds the saved graph-model npz scores the README tells a reviewer to check. No model directory is over 100 MB.

`data/raw/moleculenet_hiv.csv` is about 2.1 MB. It is not ignored and is not over 100 MB.

## Corrected training counts

Scaffold: n_train 17,813, n_val 2,653, n_test 5,306.

Cold-protein: n_train 18,020, n_val 2,584, n_test 5,168.

## Remaining blockers

- GitHub write access was not checked. The public repo was not created and was not pushed.
- Three files are over 100,000,000 bytes and must stay out of git: `data/splits_kiba/scaffold_split.csv` (107,037,499), `data/splits_kiba/cold_protein_split.csv` (102,568,146), `data/raw/kiba_pairs.csv` (101,906,920). All three are gitignored. They are still on disk.
- Personal Drive folder `virtual-screening` holds `01-code-report-artifacts.tar.gz`, modified 2026-10-03 16:11 ET. `RESEARCH_REPORT.md` in this working copy was last modified 2026-10-03 19:08 ET. That Drive archive is older than this draft. This prep did not upload a new archive.
- Not gitignored, and would become public if the tree is committed as it stands: absolute machine paths inside `artifacts/metrics_v1.2.json` (`best_ckpt`), `artifacts/metrics_v1.1.json` (`esm2_status`), `artifacts/metrics_v1.1.md`, and one editable-install line in `requirements.lock.txt`. These are machine paths, not an employer name. Research tables were not rewritten.
- `scripts/score_histgbm_cold_protein.py` depends on gitignored DAVIS pairs and `artifacts/models_v1.1/`.

## Not done

No commit. No push. No remote. No `git init`. The Mac was not synced.
