#!/bin/bash -l
#SBATCH --job-name=holstein-solver-check
#SBATCH --account=sipv
#SBATCH --qos=standby
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --output=logs/%x-%j.out
set -euo pipefail
cd /scratch/rjguo/holstein-spectral-fullgrid-20260919/holstein-diagmc
export HOLSTEIN_PYTHON=/scratch/rjguo/holstein-diagmc-kestrel-20260919/venv/bin/python
source scripts/kestrel_environment.sh
date --utc --iso-8601=seconds
./scripts/python.sh recovery/collect_failures.py
date --utc --iso-8601=seconds
