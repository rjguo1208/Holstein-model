#!/usr/bin/env bash
#SBATCH --job-name=holstein-defect
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
./scripts/python.sh scripts/run_defect.py --stage pilot --out results/defect_pilot_01 \
    --jobs 4 --chains 4 --blocks 80 --steps-per-block 100000
