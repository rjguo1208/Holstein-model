#!/usr/bin/env bash
#SBATCH --job-name=defect-dense-sample
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=02:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
./scripts/python.sh scripts/defect_dense_sample.py --out results/defect_dense_data_01 --jobs "$SLURM_CPUS_PER_TASK"
