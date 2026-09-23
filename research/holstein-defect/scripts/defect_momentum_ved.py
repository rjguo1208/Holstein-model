"""Parallel projected-resolvent VED with dense, directly contracted momenta."""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import lru_cache
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import sys
import time
import numpy as np
from scipy import sparse
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from defect_reference import Cloud, DefectBasis, Parameters
from defect_momentum import projected_lanczos, source


@lru_cache(maxsize=1)
def operator(nh, radius, coupling, U):
    basis = DefectBasis(nh, radius)
    p = Parameters(g=np.sqrt(2*coupling), U=U)
    cloud = basis.cloud
    sites = len(basis.x)
    diagonal = np.tile(p.omega*cloud.number, sites)+np.repeat(-U*(basis.x == 0), cloud.dim)
    hopping = sparse.diags(np.ones(sites-1), -1, shape=(sites, sites), format='csr')
    forward = sparse.kron(hopping, cloud.shift, format='csr')
    h = (sparse.diags(diagonal, format='csr')
         +p.g*sparse.kron(sparse.eye(sites, format='csr'), cloud.coupling, format='csr')
         -p.t*(forward+forward.T)).tocsr()
    h.sort_indices()
    # The assembled matrix must implement the independently written operator.
    probe = np.random.default_rng(94173).normal(size=basis.dim)
    np.testing.assert_allclose(h@probe, basis.hamiltonian(p)@probe, atol=2e-14, rtol=2e-14)
    return basis, h


def run(task):
    out, coupling, U, nh, radius, window, site, steps, direct = task
    started = time.monotonic()
    basis, h = operator(nh, radius, coupling, U)
    label = f'lambda{coupling:.2f}_U{U:g}_Nh{nh}_R{radius}_W{window}'
    if direct is None:
        indices = (np.arange(-window, window+1)+radius)*basis.cloud.dim
        r = projected_lanczos(h, basis.source(site), indices, steps)
        name = label+f'_i{site:02d}.npz'
        np.savez_compressed(out/name, **r)
        kind = 'column'
    else:
        fields = {}
        for parity in ['cos', 'sin']:
            initial = source(basis, direct*np.pi, window, parity)
            if np.linalg.norm(initial) < 1e-14:
                continue
            r = projected_lanczos(h, initial, np.array([], int), steps)
            fields.update({parity+'_'+k: v for k, v in r.items() if k != 'projection'})
        name = label+f'_k{direct:.3f}.npz'
        np.savez_compressed(out/name, **fields)
        kind = 'direct'
    return dict(kind=kind, coupling=coupling, U=U, generations=nh, radius=radius,
                window=window, site=site, k_over_pi=direct, steps=steps,
                dimension=basis.dim, nnz=h.nnz, seconds=time.monotonic()-started,
                file=name, sha256=sha256((out/name).read_bytes()).hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--jobs', type=int, default=16)
    parser.add_argument('--pilot', action='store_true')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    source_dir = args.out/'source'
    source_dir.mkdir()
    files = ['defect_reference.py', 'defect_momentum.py', 'scripts/defect_momentum_ved.py',
             'scripts/defect_momentum_ved.slurm', 'scripts/python.sh', 'requirements.txt',
             'tests/test_defect_momentum.py']
    hashes = {}
    for name in files:
        target = source_dir/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT/name, target)
        hashes[name] = sha256(target.read_bytes()).hexdigest()
    settings = [(10, 64, 16), (12, 64, 16), (12, 96, 32)]
    tasks = [(args.out, lam, U, nh, R, W, i, 800, None)
             for nh, R, W in settings for lam in [.25, .5] for U in [0., 1.] for i in range(W+1)]
    tasks += [(args.out, lam, U, 12, 96, 16, None, 800, k)
              for lam in [.25, .5] for U in [0., 1.] for k in [0., .125, .5, .875, 1.]]
    if args.pilot:
        tasks = [(args.out, .5, 1., 12, 64, 16, 0, 800, None)]
    record = dict(complete=False, method='VED projected real-space resolvent',
                  job_id=os.environ.get('SLURM_JOB_ID'), jobs=args.jobs,
                  source_sha256=hashes, records=[], tasks=len(tasks), pilot=args.pilot,
                  observable='coherent normalized plane wave on -W,...,W, embedded in the chain',
                  settings=settings, k_points=201, energy_points=1601,
                  eta=[.1, .15, .25, .5, 1.], steps=800, compare_steps=400)
    def save():
        tmp = args.out/'manifest.tmp'
        tmp.write_text(json.dumps(record, indent=2)+'\n')
        tmp.replace(args.out/'manifest.json')
    save()
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        for future in as_completed([pool.submit(run, task) for task in tasks]):
            item = future.result()
            record['records'].append(item)
            save()
            print('COMPLETE', len(record['records']), '/', len(tasks), item['file'],
                  f"{item['seconds']:.2f}s", flush=True)
    record['complete'] = True
    save()


if __name__ == '__main__':
    main()
