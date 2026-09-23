#!/usr/bin/env bash
#SBATCH --job-name=defect-k-ved
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=04:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
./scripts/python.sh scripts/defect_momentum_ved.py --out results/defect_momentum_ved_01 --jobs "$SLURM_CPUS_PER_TASK"
