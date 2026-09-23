"""Immutable highmem acquisition of open-endpoint momentum Green functions."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--jobs', type=int, default=16)
    parser.add_argument('--pilot', action='store_true')
    parser.add_argument('--validation', type=Path)
    args = parser.parse_args()
    # The user requested VED first, then fresh DiagMC. Never start acquisition
    # on a failed or unfinished VED computation.
    ved = ROOT/'results/defect_momentum_ved_report_01/summary.json'
    assert json.loads(ved.read_text())['complete']
    if not args.pilot:
        if args.validation is None or not json.loads(args.validation.read_text())['passed']:
            raise ValueError('Analytic open-endpoint validation required before production')
    args.out.mkdir(parents=True, exist_ok=False)
    source_dir = args.out/'source'
    source_dir.mkdir()
    names = ['src/defect_momentum.hpp', 'src/defect_momentum_main.cpp', 'src/defect.hpp',
             'src/conditional.hpp', 'src/diagmc.hpp', 'tests/defect_momentum_kernel.cpp',
             'scripts/defect_momentum_sample.py', 'scripts/python.sh', 'Makefile',
             'defect_momentum_mc.py']
    hashes = {}
    for name in names:
        target = source_dir/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT/name, target)
        hashes[name] = sha256(target.read_bytes()).hexdigest()
    executable = source_dir/'defect_momentum_diagmc'
    shutil.copy2(ROOT/'build/defect_momentum_diagmc', executable)
    cases = [('free', 1., 0., 0., -2.5, 4.),
             ('bare_defect', 1., 0., 1., -2.6, 4.),
             ('atomic', 0., .7, 1., -2.2, 4.)] if args.pilot else [
             (f'lambda{lam:.2f}', 1., float(np.sqrt(2*lam)), 1.,
              float(-np.sqrt(5)-2*lam-.1), 12.) for lam in [.25, .5]]
    tasks = []
    for index, (name, t, g, U, mu, T) in enumerate(cases):
        for chain in range(8):
            p = dict(t=t, omega=1., g=g, U=U, mu=mu, window=16,
                     **{'tau-max':T, 'bins':64 if args.pilot else 192,
                        'blocks':160 if args.pilot else 640,
                        'steps-per-block':50000 if args.pilot else 200000,
                        'warmup':500000 if args.pilot else 1000000, 'thin':50,
                        'max-arcs':256, 'max-hops':1024,
                        'seed':1193000001+index*10000000+chain*2003+(0 if args.pilot else 100000000)})
            tasks.append(dict(name=name+f'_s{chain:02d}', case=name, parameters=p))
    manifest = dict(complete=False, started_utc=datetime.now(timezone.utc).isoformat(),
                    job_id=os.environ.get('SLURM_JOB_ID'), jobs=args.jobs, pilot=args.pilot,
                    source_sha256=hashes, executable_sha256=sha256(executable.read_bytes()).hexdigest(),
                    ved_preceded_sampling_sha256=sha256(ved.read_bytes()).hexdigest(),
                    ved_used_as_spectral_prior=False, mu_source='rigorous analytic spectral lower bound minus 0.1t',
                    observable='A_L(k,w), L=33; directly sampled open endpoints on the infinite chain',
                    tasks=tasks, completed_chains=0,
                    production_steps=sum(t['parameters']['blocks']*t['parameters']['steps-per-block'] for t in tasks))
    if args.validation:
        shutil.copyfile(args.validation, args.out/'validation_input.json')
        manifest['validation_sha256'] = sha256(args.validation.read_bytes()).hexdigest()
    def save():
        tmp = args.out/'manifest.tmp'
        tmp.write_text(json.dumps(manifest, indent=2)+'\n')
        tmp.replace(args.out/'manifest.json')
    save()
    (args.out/'chains').mkdir()
    def run(task):
        path = args.out/'chains'/task['name']
        command = [str(executable)]
        for k, v in task['parameters'].items():
            command += ['--'+k, str(v)]
        command += ['--out', str(path)]
        with path.with_suffix('.log').open('w') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        meta = json.loads((path/'run.json').read_text())
        assert meta['complete'] and meta['arc_cap_attempts'] == meta['hop_cap_attempts'] == 0
        for k, v in task['parameters'].items():
            assert meta[k.replace('-', '_')] == v
        assert (path/'blocks.f64').stat().st_size == 8*meta['rows']*meta['columns']
        return task['name']
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for future in as_completed([pool.submit(run, t) for t in tasks]):
            name = future.result()
            manifest['completed_chains'] += 1
            save()
            print('COMPLETE', manifest['completed_chains'], '/', len(tasks), name, flush=True)
    manifest.update(complete=True, finished_utc=datetime.now(timezone.utc).isoformat())
    save()


if __name__ == '__main__':
    main()
