#!/usr/bin/env bash
#SBATCH --job-name=defect-benchmark
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
./scripts/python.sh scripts/benchmark_defect.py
