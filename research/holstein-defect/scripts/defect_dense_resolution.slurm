#!/usr/bin/env bash
#SBATCH --job-name=defect-resolution
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=01:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
./scripts/python.sh scripts/defect_dense_resolution.py --fit results/defect_dense_fit_01 --out results/defect_dense_resolution_01 --jobs "$SLURM_CPUS_PER_TASK" --draws 16
