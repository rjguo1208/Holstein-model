#!/usr/bin/env bash
#SBATCH --job-name=defect-dense-fit
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
./scripts/python.sh scripts/defect_dense_continue.py --data results/defect_dense_data_01 --out results/defect_dense_fit_01 --jobs "$SLURM_CPUS_PER_TASK" --bootstrap 64
