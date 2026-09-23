#!/usr/bin/env bash
#SBATCH --job-name=defect-dense-ved
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
./scripts/python.sh scripts/defect_dense_reference.py --fit results/defect_dense_fit_01 --out results/defect_dense_reference_01 --jobs "$SLURM_CPUS_PER_TASK"
