"""Freeze CV choices without VED input, then reconstruct all 201 measured k."""
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
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from defect_analysis import public, reblock
from defect_momentum_mc import read_chains, project, prepare, estimate
from defect_momentum_fit import configurations, select, bootstrap
from defect_dense import curves

ETAS = (.1, .15, .25, .5, 1.)


@lru_cache(maxsize=2)
def load(path):
    p = json.loads(path.with_suffix('.json').read_text())
    with np.load(path) as f:
        return f['blocks'], f['owners'], p, sha256(path.read_bytes()).hexdigest()


def point(task):
    input_path, out, lam, index, draws = task
    started = time.monotonic()
    blocks, owners, p, input_hash = load(input_path)
    p = dict(p, k=float(np.pi*index/200))
    rows = project(blocks, p, p['k'])
    path = out/f'lambda{lam:.2f}_k{index:03d}'
    path.mkdir(exist_ok=False)
    results = [select(rows, owners, p, config) for config in configurations()]
    chosen = int(np.argmin([r[1]['score'] for r in results]))
    fit, cv = results[chosen]
    selection = dict(k_index=index, k_over_pi=index/200, coupling=lam,
                     selected=chosen, config=cv['config'], alpha=fit['alpha'],
                     rule='within each representation use largest alpha within one fold standard error of minimum, then minimum selected score across representations',
                     reference_input=False, cv=[r[1] for r in results])
    (path/'selection.json').write_text(json.dumps(public(selection), indent=2)+'\n')
    bs, bsfits = [], []
    for j in range(draws):
        f, score = bootstrap(rows, owners, p, cv['config'],
                             seed=173700000+int(lam*100)*100000+index*101+j)
        bs.append(curves(f, ETAS))
        bsfits.append(dict(alpha=f['alpha'], energies=f['energies'], weights=f['weights'], score=score['score']))
    candidates = np.array([curves(r[0], ETAS) for r in results])
    selected = candidates[chosen]
    bs = np.array(bs)
    tau, grouped, _, pp = prepare(rows, owners, p, 2, 2)
    G, C = estimate(grouped, pp)
    errors = []
    for factor in [1, 2, 4, 8]:
        _, b, _, bp = prepare(rows, owners, p, 2, factor)
        _, cov = estimate(b, bp)
        errors.append(np.sqrt(np.diag(cov)))
    whole = np.array([rows[owners == i].sum(axis=0) for i in range(8)])
    _, b, _, bp = prepare(whole, np.arange(8), p, 2, 1)
    _, cov = estimate(b, bp)
    errors.append(np.sqrt(np.diag(cov)))
    np.savez_compressed(path/'spectra.npz', k=p['k'], energy=np.linspace(-4.5, 5.5, 1601), eta=ETAS,
                        selected=selected, candidates=candidates,
                        bootstrap_quantiles=np.quantile(bs, [.16, .5, .84], axis=0),
                        energies=fit['energies'], weights=fit['weights'],
                        tau=tau, green=G, covariance=C, blocking_errors=errors)
    for f in bsfits:
        np.testing.assert_array_equal(f['energies'], fit['energies'])
    np.savez_compressed(path/'bootstrap.npz', energies=fit['energies'],
                        weights=np.array([f['weights'] for f in bsfits]),
                        alpha=np.array([f['alpha'] for f in bsfits]),
                        score=np.array([f['score'] for f in bsfits]))
    report = dict(complete=True, coupling=lam, k_index=index, k=p['k'], parameters=p,
                  selected=chosen, config=cv['config'], fit=fit, cv=[r[1] for r in results],
                  bootstrap_draws=draws, negative_green_bins=int((G < 0).sum()),
                  green_bins_below_three_sigma=int((G < 3*np.sqrt(np.diag(C))).sum()),
                  source_groups_sha256=input_hash,
                  seconds=time.monotonic()-started, reference_used=False,
                  spectral_resolution_validated=False)
    (path/'summary.json').write_text(json.dumps(public(report), indent=2)+'\n')
    return dict(coupling=lam, k_index=index, config=cv['config']['name'], alpha=fit['alpha'],
                seconds=report['seconds'], selection_sha256=sha256((path/'selection.json').read_bytes()).hexdigest(),
                spectra_sha256=sha256((path/'spectra.npz').read_bytes()).hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--jobs', type=int, default=16)
    parser.add_argument('--bootstrap', type=int, default=32)
    parser.add_argument('--pilot', action='store_true')
    args = parser.parse_args()
    mf = json.loads((args.data/'manifest.json').read_text())
    assert mf['complete'] and not mf['pilot']
    args.out.mkdir(parents=True, exist_ok=False)
    source_dir = args.out/'source'
    source_dir.mkdir()
    files = ['defect_momentum_fit.py', 'defect_momentum_mc.py', 'defect_momentum.py',
             'defect_entropy.py', 'defect_entropy_stable.py', 'defect_dense.py', 'defect_analysis.py',
             'continuation.py', 'spectral_cv.py', 'map_continuation.py', 'spectral_nonlinear.py',
             'scripts/defect_momentum_continue.py', 'requirements.txt']
    source_hashes = {}
    for name in files:
        target = source_dir/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT/name, target)
        source_hashes[name] = sha256(target.read_bytes()).hexdigest()
    tasks, inputs = [], {}
    for lam in [.25, .5]:
        case = f'lambda{lam:.2f}'
        paths = [args.data/'chains'/t['name'] for t in mf['tasks'] if t['case'] == case]
        rows, owners, metas, hashes = read_chains(paths)
        grouped = reblock(rows, owners, 4)
        who = np.concatenate([np.full(np.sum(owners == i)//4, i) for i in range(8)])
        path = args.out/(case+'_groups.npz')
        np.savez_compressed(path, blocks=grouped, owners=who)
        p = dict(metas[0], acquisition_blocks=metas[0]['blocks'], base_block_factor=4,
                 steps_per_block=metas[0]['steps_per_block']*4)
        path.with_suffix('.json').write_text(json.dumps(p, indent=2)+'\n')
        inputs[case] = hashes
        indices = [0, 25, 100, 175, 200] if args.pilot else list(range(201))
        tasks.extend((path, args.out, lam, index, args.bootstrap) for index in indices)
    manifest = dict(complete=False, job_id=os.environ.get('SLURM_JOB_ID'), pilot=args.pilot,
                    source_hashes=source_hashes, inputs=inputs, points=[], bootstrap=args.bootstrap,
                    configurations=configurations(), reference_used=False,
                    source_manifest_sha256=sha256((args.data/'manifest.json').read_bytes()).hexdigest())
    def save():
        tmp = args.out/'manifest.tmp'
        tmp.write_text(json.dumps(public(manifest), indent=2)+'\n')
        tmp.replace(args.out/'manifest.json')
    save()
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        for f in as_completed([pool.submit(point, t) for t in tasks]):
            item = f.result()
            manifest['points'].append(item)
            save()
            print('COMPLETE', len(manifest['points']), '/', len(tasks), item, flush=True)
    manifest['complete'] = True
    save()


if __name__ == '__main__':
    main()
