#!/bin/bash -l
#SBATCH --job-name=holstein-pilot-recovery
#SBATCH --account=sipv
#SBATCH --qos=standby
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=logs/%x-%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export HOLSTEIN_PYTHON=/scratch/rjguo/holstein-diagmc-kestrel-20260919/venv/bin/python
source scripts/kestrel_environment.sh
./scripts/python.sh -m pytest -q tests
./scripts/python.sh scripts/pilot_repair_bootstrap.py
