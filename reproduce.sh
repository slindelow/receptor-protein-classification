#!/usr/bin/env bash
# Rebuild env (if needed), retrain/eval, verify metrics within documented tolerance.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

if [[ ! -d .venv ]]; then
  uv venv --python 3.12 .venv
fi
# shellcheck disable=SC1091
source .venv/bin/activate

uv pip install -e . --quiet
uv pip install "PyTDC==0.4.1" --no-deps --quiet || true
uv pip install "setuptools==75.1.0" --quiet || true

REF_COPY="$ROOT/artifacts/.metrics_reference.json"
if [[ -f "$ROOT/artifacts/metrics.json" && ! -f "$REF_COPY" ]]; then
  cp "$ROOT/artifacts/metrics.json" "$REF_COPY"
fi

python scripts/run_train_eval.py
python scripts/make_demo_library.py

if [[ -f "$REF_COPY" ]]; then
  python scripts/check_metrics_tolerance.py "$REF_COPY" "$ROOT/artifacts/metrics.json"
else
  # First run: snapshot as reference for subsequent reproduces
  cp "$ROOT/artifacts/metrics.json" "$REF_COPY"
  python scripts/check_metrics_tolerance.py "$REF_COPY" "$ROOT/artifacts/metrics.json"
fi

echo "OK: artifacts/metrics.md and artifacts/metrics.json updated"
echo "Sample ranked: artifacts/sample_ranked.csv"
