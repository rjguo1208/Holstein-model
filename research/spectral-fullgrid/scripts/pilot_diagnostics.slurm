#!/bin/bash -l
#SBATCH --job-name=holstein-pilot-diagnostics
#SBATCH --account=sipv
#SBATCH --qos=standby
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x-%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export HOLSTEIN_PYTHON=/scratch/rjguo/holstein-diagmc-kestrel-20260919/venv/bin/python
source scripts/kestrel_environment.sh
date --utc --iso-8601=seconds
./scripts/python.sh scripts/pilot_diagnostics.py
date --utc --iso-8601=seconds
