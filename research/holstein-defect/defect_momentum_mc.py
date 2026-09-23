"""Read open-endpoint bare DiagMC and retain signed momentum measurements."""
from __future__ import annotations
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
from defect_analysis import reblock, jackknife


def read_chains(paths):
    rows, owners, metadata, hashes = [], [], [], {}
    keys = ['t', 'omega', 'g', 'U', 'mu', 'window', 'tau_max', 'bins', 'thin']
    for owner, path in enumerate(map(Path, paths)):
        m = json.loads((path/'run.json').read_text())
        if not m['complete'] or m['arc_cap_attempts'] or m['hop_cap_attempts']:
            raise ValueError('Incomplete or capped sampling')
        if metadata and any(m[k] != metadata[0][k] for k in keys):
            raise ValueError('Incompatible chains')
        a = np.fromfile(path/'blocks.f64', dtype='<f8').reshape(m['rows'], m['columns'])
        if not np.isfinite(a).all() or np.any(a < 0):
            raise ValueError('Invalid positive real-space histograms')
        np.testing.assert_allclose(a[:, 2:].sum(axis=1), a[:, 1], atol=1e-8, rtol=1e-11)
        assert np.all(a[:, 1] == m['steps_per_block']//m['thin'])
        rows.append(a)
        owners.extend([owner]*len(a))
        metadata.append(m)
        for name in ['run.json', 'blocks.f64']:
            hashes[path.name+'/'+name] = sha256((path/name).read_bytes()).hexdigest()
    assert len({m['seed'] for m in metadata}) == len(metadata)
    return np.concatenate(rows), np.array(owners), metadata, hashes


def reference_integral(p):
    rates = np.array([-p['mu'], -p['mu']-p['U']])
    integrals = -np.expm1(-rates*p['tau_max'])/rates
    return (2*p['window']*integrals[0]+integrals[1])/(2*p['window']+1)


def project(rows, p, k):
    """Cosine transform of measured endpoint displacement; keeps every bin."""
    count = rows[:, 2:].reshape(len(rows), 2*p['window']+1, p['bins'])
    signed = np.einsum('bdt,d->bt', count, np.cos(k*np.arange(2*p['window']+1)))
    return np.column_stack([rows[:, :2], signed])


def prepare(rows, owners, p, factor=2, block=4):
    n = p['bins']
    if n % factor or rows.shape[1] != n+2:
        raise ValueError('Unexpected time bins')
    compact = np.column_stack([rows[:, :2], rows[:, 2:].reshape(len(rows), n//factor, factor).sum(axis=2)])
    grouped = reblock(compact, owners, block)
    who = np.concatenate([np.full(np.sum(owners == i)//block, i) for i in np.unique(owners)])
    pp = dict(p, bins=n//factor)
    tau = (np.arange(pp['bins'])+.5)*pp['tau_max']/pp['bins']
    return tau, grouped, who, pp


def normalize(total, p):
    if np.any(total[..., 0] <= 0):
        raise ValueError('Insufficient zero-vertex normalization')
    return total[..., 2:]/total[..., 0, None]*reference_integral(p)/(p['tau_max']/p['bins'])


def estimate(rows, p):
    return jackknife(rows, lambda total: normalize(total, p))[:2]


def variance_whitening(cov, tolerance=1e-8):
    """Correlation-matrix conditioning works even for noisy negative G bins."""
    scale = np.sqrt(np.maximum(np.diag(cov), np.finfo(float).tiny))
    correlation = cov/np.outer(scale, scale)
    values, vectors = np.linalg.eigh(correlation)
    use = values > max(values[-1]*tolerance, 0.)
    return (vectors[:, use]/np.sqrt(values[use])).T/scale[None, :]
