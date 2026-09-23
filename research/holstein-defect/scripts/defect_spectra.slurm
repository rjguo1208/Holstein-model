#!/usr/bin/env bash
#SBATCH --job-name=holstein-defect-ldos
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
./scripts/python.sh scripts/defect_spectra.py --ground results/defect_scan_01 \
    --out results/defect_spectra_01 --jobs 4 --blocks 160 --bootstrap 8
