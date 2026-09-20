#!/bin/bash -l
#SBATCH --job-name=holstein-full-refinement
#SBATCH --account=sipv
#SBATCH --qos=standby
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=104
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
set -euo pipefail
cd /scratch/rjguo/holstein-spectral-fullgrid-20260919/holstein-diagmc
export HOLSTEIN_PYTHON=/scratch/rjguo/holstein-diagmc-kestrel-20260919/venv/bin/python
source scripts/kestrel_environment.sh
date --utc --iso-8601=seconds
hostname
./scripts/python.sh - <<'PYENV'
import json, os, platform, sys
from pathlib import Path
import numpy, scipy, matplotlib
Path('results/environment.json').write_text(json.dumps(dict(python=sys.version,
    platform=platform.platform(),libc=platform.libc_ver(),host=platform.node(),
    numpy=numpy.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__,
    slurm={k:os.environ.get(k) for k in ['SLURM_JOB_ID','SLURM_JOB_PARTITION','SLURM_CPUS_PER_TASK','SLURM_JOB_ACCOUNT']}),indent=2)+'\n')
PYENV
./build/test_kernel
./scripts/python.sh -m pytest -q tests
./scripts/python.sh scripts/full_grid_pipeline.py
date --utc --iso-8601=seconds
