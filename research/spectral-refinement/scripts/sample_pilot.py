"""Acquire an independent four-times-statistics pilot, preserving the baseline."""
import concurrent.futures as cf
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
plan = json.loads((ROOT/'plans/pilot.json').read_text())
assert hashlib.sha256((ROOT/'build/diagmc').read_bytes()).hexdigest() == plan['baseline_executable_sha256']
cpus = int(os.environ.get('SLURM_CPUS_PER_TASK', '1'))
chain_jobs = min(8, cpus)

def run(ik):
    cmd = [sys.executable, str(ROOT/'scripts/run_suite.py'), 'spectral_map',
           '--map-points', '41', '--map-index', str(ik), '--out',
           str(ROOT/f'results/pilot_data/ik{ik:03d}'), '--couplings', '.25', '.5',
           '--chains', str(plan['chains']), '--blocks', str(plan['blocks']),
           '--steps-per-block', str(plan['steps_per_block']),
           '--seed-offset', str(plan['seed_offset']), '--rb-thin', str(plan['rb_thin']),
           '--jobs', str(chain_jobs)]
    with (ROOT/f'logs/pilot-ik{ik:03d}.out').open('x') as log:
        subprocess.run(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
    print(f'Completed independent pilot k-index {ik}', flush=True)

with cf.ThreadPoolExecutor(max_workers=max(1, cpus//chain_jobs)) as pool:
    list(pool.map(run, plan['indices']))
runs = [json.loads(p.read_text()) for p in (ROOT/'results/pilot_data').glob('ik*/lambda*/run.json')]
old = [json.loads(p.read_text()) for p in Path(plan['baseline']).glob('ik*/lambda*/run.json')]
assert len(runs) == 160 and all(r['complete'] for r in runs)
assert len({r['seed'] for r in runs}) == len(runs)
assert not ({r['seed'] for r in runs} & {r['seed'] for r in old})
assert sum(r['blocks']*r['steps_per_block'] for r in runs) == plan['new_production_steps']
record = dict(job_id=os.environ.get('SLURM_JOB_ID'), complete=True, chains=len(runs),
              production_steps=plan['new_production_steps'], unique_disjoint_seeds=True,
              cap_attempts=sum(r['cap_attempts'] for r in runs))
(ROOT/'results/pilot_sampling.json').write_text(json.dumps(record, indent=2)+'\n')
print(json.dumps(record), flush=True)
