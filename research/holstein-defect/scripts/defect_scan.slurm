#!/usr/bin/env bash
#SBATCH --job-name=holstein-defect-scan
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
./scripts/python.sh scripts/run_defect.py --stage scan --pilot results/defect_pilot_01 \
    --out results/defect_scan_01 --jobs 4 --chains 8 --blocks 320 --steps-per-block 100000 \
    --seed-offset 10000000
