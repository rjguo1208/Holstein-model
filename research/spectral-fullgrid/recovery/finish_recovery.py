"""Complete original failed bootstrap draws after independently checked recovery.

The original 10,478 successful draws and all central fits are immutable. The
original solver is retained wherever it converges; the stable numerical
fallback solves only remaining draws at their original selected alpha.
"""
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'scripts')]
from analysis import read_runs, estimate_green
from refinement import prepare, MAXENT_ALPHAS
from pilot_compare import curves
from pilot_diagnostics import public
from stable_maxent import solve, analytic_checks

DIAG = ROOT / 'results/numerical_recovery/attempt1'
OUT = ROOT / 'results/numerical_recovery/final'


def write(path, data):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(public(data), indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def verify_fit(fit):
    w = np.asarray(fit['weights'])
    assert fit['converged'] and fit['dual_gradient'] < 1e-8
    assert np.isfinite(w).all() and (w >= 0).all()
    assert abs(w.sum() - 1.) < 1e-12
    assert max(abs(np.asarray(fit['moment_residuals']))) < 1e-6


def recover(path):
    record = json.loads(path.read_text())
    attempt = record['attempted_recovery']
    assert len(attempt['alphas']) == 1
    alpha = attempt['alphas'][0]
    options = {k: record['config'][k] for k in ['size', 'upper', 'tau_limit', 'tolerance']}
    with np.load(path.with_name(path.stem + '-input.npz')) as data:
        fit = solve(data['tau'], data['green'], data['covariance'],
                    attempt['parameters'], alpha, **options)
    verify_fit(fit)
    assert fit['alpha'] == attempt['fits'][0]['alpha']
    record.update(success=True, alpha=alpha, fit=fit,
                  selected_cv=int(np.flatnonzero(MAXENT_ALPHAS == alpha)[0]),
                  original_failure=record['error'], fallback_used=True)
    np.savez_compressed(OUT / (path.stem + '.npz'), curves=curves(fit))
    write(OUT / path.name, record)
    print(f'Recovered {path.stem}: gradient={fit["dual_gradient"]:.3g}, '
          f'moment_error={max(abs(fit["moment_residuals"])):.3g}', flush=True)
    return dict(stem=path.stem, alpha=alpha, gradient=fit['dual_gradient'],
                moment_error=float(max(abs(fit['moment_residuals']))))


def control(path):
    record = json.loads(path.read_text())
    lam, ik, rep = record['coupling'], record['index'], record['replicate']
    paths = sorted((ROOT / 'results/pilot_data').glob(f'ik{ik:03d}/lambda{lam:.2f}_*/run.json'))
    rows, owners, metas = read_runs([p.parent for p in paths], 'rb_blocks.csv')
    rng = np.random.default_rng(912731 + int(lam * 100) * 1000 + ik)
    for _ in range(rep + 1):
        blocks = []
        for owner in np.unique(owners):
            group = rows[owners == owner]
            chunks = group.reshape(-1, 16, group.shape[1])
            blocks.append(chunks[rng.integers(len(chunks), size=len(chunks))].reshape(group.shape))
        sample = np.concatenate(blocks)
    assert sha256(sample.tobytes()).hexdigest() == record['resampled_rows_sha256']
    tau, blocks, _, pp = prepare(sample, owners, metas[0], record['config'])
    G, C, _ = estimate_green(blocks, pp)
    options = {k: record['config'][k] for k in ['size', 'upper', 'tau_limit', 'tolerance']}
    fit = solve(tau, G, C, pp, record['alpha'], **options)
    verify_fit(fit)
    difference = float(np.max(abs(curves(fit) - curves(record['fit']))))
    assert difference < 1e-5, (path.name, difference)
    result = dict(stem=path.stem, alpha=record['alpha'],
                  maximum_broadened_curve_difference=difference,
                  gradient=fit['dual_gradient'], moment_error=float(max(abs(fit['moment_residuals']))),
                  original_successful_draw_unchanged=True)
    write(OUT / (path.stem + '-control.json'), result)
    print('Control passed: ' + json.dumps(result), flush=True)
    return result


def main():
    assert os.environ.get('SLURM_JOB_ID'), 'Compute allocation required'
    recovery_source = json.loads((ROOT / 'results/numerical_recovery-source.json').read_text())
    for name, expected in recovery_source['sha256'].items():
        assert sha256((ROOT / name).read_bytes()).hexdigest() == expected, name
    diagnostic = json.loads((DIAG / 'summary.json').read_text())
    assert diagnostic['complete']
    assert diagnostic['source_sha256'] == sha256((ROOT / 'recovery/collect_failures.py').read_bytes()).hexdigest()
    OUT.mkdir(parents=True, exist_ok=False)
    frozen = json.loads((ROOT / 'results/inference_source.json').read_text())
    for name, expected in frozen['sha256'].items():
        assert sha256((ROOT / name).read_bytes()).hexdigest() == expected, name
    analytic = analytic_checks()
    write(OUT / 'analytic-validation.json', analytic)
    paths = sorted(DIAG.glob('lambda*.json'))
    failed = [p for p in paths if not json.loads(p.read_text())['success']]
    successful = [p for p in paths if json.loads(p.read_text())['success']]
    # Select controls spanning available alpha values before solving new cases.
    controls = []
    seen = set()
    for path in successful:
        alpha = json.loads(path.read_text())['alpha']
        if alpha not in seen:
            controls.append(path)
            seen.add(alpha)
    controls = controls[:3]
    assert controls
    with ProcessPoolExecutor(max_workers=min(8, int(os.environ['SLURM_CPUS_PER_TASK']))) as pool:
        recovered = list(pool.map(recover, failed))
        validated_controls = list(pool.map(control, controls))

    initial_path = ROOT / 'results/bootstrap_initial/summary.json'
    initial = json.loads(initial_path.read_text())
    to_write = []
    recovery_points = []
    for point in initial['points']:
        if not point['failures']:
            continue
        lam, ik = point['coupling'], point['index']
        stem = f'lambda{lam:.2f}_ik{ik:03d}'
        backup = ROOT / 'results/bootstrap_initial'
        for ext in ['json', 'npz']:
            assert (ROOT / 'results/bootstrap' / f'{stem}.{ext}').read_bytes() == (backup / f'{stem}.{ext}').read_bytes()
        with np.load(backup / f'{stem}.npz') as saved:
            data = {k: saved[k] for k in saved.files}
        failures = sorted(f['replicate'] for f in point['failures'])
        successes = [r for r in range(128) if r not in failures]
        samples = np.empty((128,) + data['samples'].shape[1:])
        alphas = np.empty(128)
        samples[successes] = data['samples']
        alphas[successes] = data['alphas']
        records = {}
        for rep in failures:
            name = f'{stem}_rep{rep:03d}'
            folder = OUT if (OUT / (name + '.json')).exists() else DIAG
            record = json.loads((folder / (name + '.json')).read_text())
            assert record['success']
            verify_fit(record['fit'])
            with np.load(folder / (name + '.npz')) as saved:
                samples[rep] = saved['curves']
            alphas[rep] = record['alpha']
            records[str(rep)] = {k: record[k] for k in ['alpha', 'fit', 'selected_cv', 'resampled_rows_sha256']}
            records[str(rep)]['numerical_method'] = 'stable fallback' if record.get('fallback_used') else 'original extended-precision recovery'
        np.testing.assert_array_equal(samples[successes], data['samples'])
        np.testing.assert_array_equal(alphas[successes], data['alphas'])
        data.update(samples=samples, alphas=alphas,
                    quantiles=np.quantile(samples, [.025, .16, .84, .975], axis=0))
        meta = json.loads((backup / (stem + '.json')).read_text())
        meta.update(initial_failures=meta['failures'], failures=[], successful=128,
                    numerical_recoveries=records)
        to_write.append((stem, data, meta))
        point.update(meta)
        recovery_points.append(dict(coupling=lam, index=ik, replicates=failures))
    assert len(paths) == sum(len(p['replicates']) for p in recovery_points) == 18
    # Mutate final outputs only after every fit and every control has passed.
    for stem, data, meta in to_write:
        np.savez_compressed(ROOT / 'results/bootstrap' / (stem + '.npz'), **data)
        write(ROOT / 'results/bootstrap' / (stem + '.json'), meta)
    initial.update(recovery_job_id=os.environ['SLURM_JOB_ID'], numerical_recoveries=recovery_points)
    write(ROOT / 'results/bootstrap/summary.json', initial)
    source = {str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest()
              for p in (ROOT / 'recovery').glob('*') if p.is_file()}
    recovery = dict(complete=True, job_id=os.environ['SLURM_JOB_ID'],
                    utc=datetime.now(timezone.utc).isoformat(), recoveries=recovery_points,
                    source_hashes=source, unchanged_central_selections=True,
                    unchanged_successful_bootstrap_draws=True, reference_used=False,
                    original_extended_precision_successes=len(successful),
                    stable_solver_recoveries=recovered, analytic_validation=analytic,
                    independent_controls=validated_controls)
    write(ROOT / 'results/bootstrap_recovery.json', recovery)
    write(OUT / 'summary.json', recovery)
    # The unmodified pipeline rechecks source/reuse/raw data and all 82 draws.
    from full_grid_pipeline import main as original_pipeline
    original_pipeline()


if __name__ == '__main__':
    main()
