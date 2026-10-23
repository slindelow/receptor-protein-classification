#!/usr/bin/env bash
# Sequence-identity holdouts, censored DAVIS labels, and the additive baseline.
# Does not rewrite artifacts/metrics_v1.1.json or artifacts/metrics_uncertainty.json.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4

if [[ ! -x .venv/bin/python ]]; then
  python3 -m venv .venv
fi
# shellcheck disable=SC1091
source .venv/bin/activate

python -m pip install \
  pandas==2.2.3 \
  numpy==1.26.4 \
  scikit-learn==1.5.2 \
  rdkit==2024.3.5 \
  scipy==1.13.1 \
  joblib==1.4.2 \
  pyyaml==6.0.2

python scripts/run_reviewer_v1.7.py
