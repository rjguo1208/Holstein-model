#!/usr/bin/env bash
#SBATCH --job-name=defect-k-mc
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
./scripts/python.sh scripts/defect_momentum_sample.py --out results/defect_momentum_mc_data_01 --jobs "$SLURM_CPUS_PER_TASK" --validation results/defect_momentum_mc_validation_01.json
