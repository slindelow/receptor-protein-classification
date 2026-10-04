#!/usr/bin/env bash
# Re-run sofia-vs-v1.1-protein ablation; verify metrics_v1.1 within tol 1e-4.
# Does NOT regenerate data/splits; does NOT overwrite artifacts/metrics.md (MVP).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

if [[ ! -d .venv ]]; then
  uv venv --python 3.12 .venv
fi
# shellcheck disable=SC1091
source .venv/bin/activate

uv pip install -e . --quiet
# CPU torch preferred for ESM-2
uv pip install torch --index-url https://download.pytorch.org/whl/cpu --quiet || true
uv pip install "fair-esm==2.0.0" --quiet || true
uv pip install "PyTDC==0.4.1" --no-deps --quiet || true
uv pip install "setuptools==78.1.0" --quiet || true

REF_COPY="$ROOT/artifacts/.metrics_v1.1_reference.json"
if [[ -f "$ROOT/artifacts/metrics_v1.1.json" && ! -f "$REF_COPY" ]]; then
  cp "$ROOT/artifacts/metrics_v1.1.json" "$REF_COPY"
fi

python scripts/run_protein_ablation.py configs/v1.1_protein.yaml

if [[ -f "$REF_COPY" ]]; then
  python scripts/check_metrics_tolerance.py "$REF_COPY" "$ROOT/artifacts/metrics_v1.1.json"
else
  cp "$ROOT/artifacts/metrics_v1.1.json" "$REF_COPY"
  python scripts/check_metrics_tolerance.py "$REF_COPY" "$ROOT/artifacts/metrics_v1.1.json"
fi

echo "OK: artifacts/metrics_v1.1.md and artifacts/metrics_v1.1.json updated"
echo "MVP baseline untouched: artifacts/metrics.md"
