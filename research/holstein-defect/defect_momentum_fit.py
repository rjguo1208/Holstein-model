"""Reference-free continuation of signed, nonlocal DiagMC measurements."""
from __future__ import annotations
import numpy as np
from continuation import bin_kernel
from defect_entropy import entropy_scan
from defect_momentum import window_moments
from defect_momentum_mc import prepare, estimate, variance_whitening

ALPHAS = np.logspace(-2, 4, 7)


def configurations():
    # Input blocks already contain four original acquisition blocks (0.8M
    # steps). The default factor two below makes 1.6M-step statistical blocks.
    base = dict(factor=2, block=2, size=641, upper=10., tau_limit=6., tolerance=1e-8)
    return [dict(base, name=name)|changes for name, changes in [
        ('primary', {}), ('time_coarse', dict(factor=4)),
        ('time_long', dict(tau_limit=8.)), ('energy_coarse', dict(size=321)),
        ('block_long', dict(block=4)), ('covariance_cut', dict(tolerance=1e-6))]]


def fit_scan(tau, G, C, p, config):
    moments = window_moments([p['k']], p['window'], p['t'], p['g'], p['U'])[0]
    options = {key:config[key] for key in ['size', 'upper', 'tau_limit', 'tolerance']}
    fs = entropy_scan(tau, G, C, p, ALPHAS, target_moments=moments,
                      whitening_mode='variance', max_iterations=600, **options)
    assert all(f['converged'] for f in fs)
    for f in fs:
        f['method'] = 'maximum entropy with exact coherent-window moments and variance-scaled full covariance'
    return fs


def select_prepared(training, p, config, validation):
    tau, blocks, who, pp = training
    tv, vb, vw, vp = validation
    use = tv <= 6.
    folds = []
    for fold in range(4):
        G, C = estimate(blocks[who//2 != fold], pp)
        test, VC = estimate(vb[vw//2 == fold], vp)
        W = variance_whitening(VC[np.ix_(use, use)], 1e-8)
        fs = fit_scan(tau, G, C, pp, config)
        K = bin_kernel(tv[use], vp['tau_max']/vp['bins'], fs[0]['energies'], pp['mu'])
        predicted = np.array([f['weights'] for f in fs])@K.T
        residual = (predicted-test[use])@W.T
        folds.append(np.sum(residual**2, axis=1)/len(W))
    folds = np.array(folds)
    scores = folds.mean(axis=0)
    best = int(np.argmin(scores))
    se = float(folds[:, best].std(ddof=1)/2)
    selected = int(np.flatnonzero(scores <= scores[best]+se)[-1])
    G, C = estimate(blocks, pp)
    fs = fit_scan(tau, G, C, pp, config)
    return fs[selected], dict(config=config, alphas=ALPHAS, fold_scores=folds,
                             score=float(scores[selected]), mean_scores=scores,
                             selected_alpha_index=selected, minimum_alpha_index=best,
                             fold_standard_error=se)


def select(rows, owners, p, config):
    training = prepare(rows, owners, p, config['factor'], config['block'])
    validation = prepare(rows, owners, p, 4, 2)
    return select_prepared(training, p, config, validation)


def bootstrap(rows, owners, p, config, seed):
    _, base, who, pp = prepare(rows, owners, p, 1, config['block'])
    rng = np.random.default_rng(seed)
    sampled = np.concatenate([b[rng.integers(len(b), size=len(b))]
                              for b in [base[who == i] for i in range(8)]])
    changed = dict(config, block=1)
    # Blocking has already occurred before resampling. Preserve chain labels
    # so each validation fold still holds out two complete independent chains.
    return select_prepared(prepare(sampled, who, pp, changed['factor'], 1), pp, changed,
                           prepare(sampled, who, pp, 4, 1))
