#!/usr/bin/env bash
# Re-run sofia-vs-v1.4-selectivity (SRC/LCK cold-pair HistGBM).
# Does NOT regenerate DAVIS splits; does NOT touch prior metrics exports.
# Two HistGBM fits on ~25k pairs; expect ~10 min CPU.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
# shellcheck disable=SC1091
source .venv/bin/activate

REF_COPY="$ROOT/artifacts/.metrics_v1.4_reference.json"
if [[ -f "$ROOT/artifacts/metrics_selectivity.json" && ! -f "$REF_COPY" ]]; then
  cp "$ROOT/artifacts/metrics_selectivity.json" "$REF_COPY"
fi

python scripts/run_selectivity_v1.4.py configs/v1.4_selectivity.yaml

if [[ -f "$REF_COPY" ]]; then
  python scripts/check_selectivity_tolerance.py "$REF_COPY" "$ROOT/artifacts/metrics_selectivity.json" 1e-4
else
  cp "$ROOT/artifacts/metrics_selectivity.json" "$REF_COPY"
  echo "Wrote reference $REF_COPY (first run; no prior comparison)"
fi
echo "OK: artifacts/metrics_selectivity.md + .json"
