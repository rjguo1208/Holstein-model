#!/usr/bin/env bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MPLBACKEND=Agg
exec "${HOLSTEIN_PYTHON:-/anvil/scratch/x-rg47749/nnqs_learning/venv/bin/python}" "$@"
