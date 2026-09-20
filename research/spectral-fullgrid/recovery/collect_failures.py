"""Diagnose failed draws with the original frozen solver, saving each outcome.

This does not change any bootstrap or central-fit output. Heavy work requires
a Slurm allocation. Neither reference spectra nor model-selection rules change.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import sys
import traceback

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'scripts')]
from analysis import read_runs
import refinement
from pilot_compare import curves
from pilot_diagnostics import public

OUT = ROOT / 'results/numerical_recovery/attempt1'


def write(path, data):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(public(data), indent=2) + '\n')
    temporary.replace(path)


def run(task):
    lam, ik, failed = task
    selection = json.loads((ROOT / f'results/comparison/lambda{lam:.2f}_ik{ik:03d}.json').read_text())
    config = selection['primary_cv']['config']
    paths = sorted((ROOT / 'results/pilot_data').glob(f'ik{ik:03d}/lambda{lam:.2f}_*/run.json'))
    rows, owners, metas = read_runs([p.parent for p in paths], 'rb_blocks.csv')
    rng = np.random.default_rng(912731 + int(lam * 100) * 1000 + ik)
    original = refinement.entropy_scan
    current = {}

    def capture(tau, green, cov, p, alphas, **kwargs):
        fits = original(tau, green, cov, p, alphas, **kwargs)
        if kwargs.get('extended_precision'):
            current.update(tau=tau, green=green, covariance=cov, parameters=p,
                           alphas=alphas, options=kwargs, fits=fits)
        return fits

    refinement.entropy_scan = capture
    outcomes = []
    for rep in range(max(failed) + 1):
        blocks = []
        for owner in np.unique(owners):
            group = rows[owners == owner]
            chunks = group.reshape(-1, 16, group.shape[1])
            blocks.append(chunks[rng.integers(len(chunks), size=len(chunks))].reshape(group.shape))
        if rep not in failed:
            continue
        sample = np.concatenate(blocks)
        stem = f'lambda{lam:.2f}_ik{ik:03d}_rep{rep:03d}'
        record = dict(coupling=lam, index=ik, replicate=rep,
                      resampled_rows_sha256=sha256(sample.tobytes()).hexdigest(),
                      config=config, original_parameters=metas[0],
                      job_id=os.environ['SLURM_JOB_ID'])
        current.clear()
        try:
            fit, cv = refinement.cross_validate(sample, owners, metas[0], config, recover=True)
            assert fit['converged'] and max(abs(fit['moment_residuals'])) < 1e-6
            np.savez_compressed(OUT / (stem + '.npz'), curves=curves(fit))
            record.update(success=True, alpha=fit['alpha'], fit=fit, selected_cv=cv['selected'])
        except Exception as error:
            record.update(success=False, error=str(error), traceback=traceback.format_exc())
            if current:
                np.savez_compressed(OUT / (stem + '-input.npz'),
                                    tau=current['tau'], green=current['green'],
                                    covariance=current['covariance'])
                record['attempted_recovery'] = {
                    k: v for k, v in current.items() if k not in ['tau', 'green', 'covariance']}
        write(OUT / (stem + '.json'), record)
        print(json.dumps({k: record[k] for k in ['coupling', 'index', 'replicate', 'success']}), flush=True)
        outcomes.append(record)
    return outcomes


if __name__ == '__main__':
    assert os.environ.get('SLURM_JOB_ID'), 'Compute allocation required'
    OUT.mkdir(parents=True, exist_ok=False)
    frozen = json.loads((ROOT / 'results/inference_source.json').read_text())
    for name, expected in frozen['sha256'].items():
        assert sha256((ROOT / name).read_bytes()).hexdigest() == expected, name
    points = json.loads((ROOT / 'results/bootstrap/summary.json').read_text())['points']
    tasks = [(p['coupling'], p['index'], [f['replicate'] for f in p['failures']])
             for p in points if p['failures']]
    records = []
    with ProcessPoolExecutor(max_workers=min(len(tasks), int(os.environ['SLURM_CPUS_PER_TASK']))) as pool:
        for future in as_completed([pool.submit(run, t) for t in tasks]):
            records.extend(future.result())
    result = dict(complete=True, job_id=os.environ['SLURM_JOB_ID'],
                  utc=datetime.now(timezone.utc).isoformat(),
                  source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                  reference_used=False, frozen_source_unchanged=True,
                  successful=sum(r['success'] for r in records), failed=sum(not r['success'] for r in records),
                  points=[{k: r[k] for k in ['coupling', 'index', 'replicate', 'success']} for r in records])
    write(OUT / 'summary.json', result)
    print(json.dumps(result), flush=True)
