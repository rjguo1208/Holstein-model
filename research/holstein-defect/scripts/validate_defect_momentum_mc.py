"""Validate normalization and nonlocal momentum phases using analytic limits."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import numpy as np
from scipy.linalg import eigh_tridiagonal
from numpy.polynomial.legendre import leggauss
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from defect_momentum_mc import read_chains, project, prepare, estimate
from defect_analysis import public


def expected(p, k, tau, dt):
    nodes, weights = leggauss(24)
    time = tau[:, None]+dt*nodes/2
    if p['t'] == 0:
        exponent = (p['mu']+p['g']**2/p['omega'])*time-(p['g']/p['omega'])**2*(-np.expm1(-p['omega']*time))
        green = np.exp(exponent)*(2*p['window']+np.exp(p['U']*time))/(2*p['window']+1)
    else:
        radius = 80
        x = np.arange(-radius, radius+1)
        energy, states = eigh_tridiagonal(-p['U']*(x == 0).astype(float), -p['t']*np.ones(2*radius))
        initial = np.zeros(len(x), complex)
        use = abs(x) <= p['window']
        initial[use] = np.exp(1j*k*x[use])/np.sqrt(2*p['window']+1)
        residue = abs(states.T@initial)**2
        green = np.einsum('n,nbq->bq', residue, np.exp(-(energy[:, None, None]-p['mu'])*time))
    return green@weights/2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    mf = json.loads((args.data/'manifest.json').read_text())
    assert mf['complete'] and mf['pilot']
    records = []
    for case in ['free', 'bare_defect', 'atomic']:
        paths = [args.data/'chains'/t['name'] for t in mf['tasks'] if t['case'] == case]
        rows, owners, metas, hashes = read_chains(paths)
        p = metas[0]
        for kk in [0., .125, .5, .875, 1.]:
            tau, blocks, _, pp = prepare(project(rows, p, kk*np.pi), owners, p, factor=2, block=4)
            G, C = estimate(blocks, pp)
            exact = expected(pp, kk*np.pi, tau, pp['tau_max']/pp['bins'])
            sigma = np.sqrt(np.diag(C))
            z = (G-exact)/sigma
            # A predeclared 8-sigma maximum allows many strongly correlated
            # bins without accepting a persistent percent-level discrepancy.
            assert np.isfinite(z).all() and abs(z).max() < 8., (case, kk, abs(z).max())
            records.append(dict(case=case, k_over_pi=kk, max_standardized_error=float(abs(z).max()),
                                median_standardized_error=float(np.median(abs(z))),
                                tau=tau, green=G, exact=exact, standard_error=sigma))
        if p['t']:
            assert all(m['nonlocal_samples'] > 0 and m['odd_displacement_samples'] > 0 for m in metas)
        else:
            assert all(m['nonlocal_samples'] == m['odd_displacement_samples'] == 0 for m in metas)
    result = dict(passed=True, analytic_cases=3, independent_chains=24,
                  production_steps=mf['production_steps'], maximum_allowed_standardized_error=8., records=records)
    if args.out.exists():
        raise ValueError('Refusing to overwrite validation')
    args.out.write_text(json.dumps(public(result), indent=2)+'\n')
    print('PASSED analytic open-endpoint normalization and phases;',
          'max |z| =', max(r['max_standardized_error'] for r in records), flush=True)


if __name__ == '__main__':
    main()
