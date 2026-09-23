"""Positive full-spectrum continuation without an imposed quasiparticle pole.

This is a preliminary inverse-Laplace estimator, not a resolution guarantee.
All k values use the same physical energy origin and exact unbroadened moments.
The reference map is never an input to selection or inference.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import nnls
from analysis import estimate_green
from continuation import bin_kernel
from spectral_cv import whitening


def lorentz_map(axis, energies, weights, eta):
    """Return (..., energy) densities, preserving absolute spectral weights."""
    if eta <= 0:
        raise ValueError('eta must be positive')
    kernel = eta/np.pi/((np.asarray(axis)[:, None]-energies)**2+eta**2)
    return np.asarray(weights) @ kernel.T


def rebin_conditional(blocks, parameters, factor=2):
    n = parameters['bins']
    if n % factor or blocks.shape[1] != 3+n:
        raise ValueError('Expected conditional count blocks with divisible time bins')
    result = np.column_stack([blocks[:, :3],
        blocks[:, 3:].reshape(len(blocks), n//factor, factor).sum(axis=2)])
    p = dict(parameters, bins=n//factor)
    tau = (np.arange(p['bins'])+.5)*p['tau_max']/p['bins']
    return tau, result, p


def scan(tau, green, cov, parameters, alphas, size=281, upper=10.,
         energies=None, tau_limit=8., tolerance=1e-8, lower_bound=None, moments=None):
    p = parameters
    # H >= -2t - g^2/omega in the one-electron sector, by completing the square.
    # This lower bound is independent of VED and of any fitted E(k), Z(k).
    lower = -2*p['t']-p['g']**2/p['omega'] if lower_bound is None else float(lower_bound)
    energies = np.linspace(lower, upper, size) if energies is None else np.asarray(energies)
    if not np.all(np.diff(energies) > 0):
        raise ValueError('Need an increasing spectral grid')
    if not np.allclose(np.diff(energies), np.diff(energies)[0]):
        raise ValueError('Regularizer requires a uniform energy grid')
    use = tau <= tau_limit
    if np.any(green[use] <= 0) or np.sum(use) < 4:
        raise ValueError('Insufficient positive imaginary-time data')
    W = whitening(green[use], cov[np.ix_(use, use)], tolerance)
    if len(W) < 3:
        raise ValueError('Insufficient covariance rank')
    K = bin_kernel(tau, p['tau_max']/p['bins'], energies, p['mu'])
    M = np.vstack([np.ones(len(energies)), energies, energies**2])
    if moments is None:
        eps = -2*p['t']*np.cos(p['k'])
        moments = np.array([1., eps, eps**2+p['g']**2])
    moments = np.asarray(moments, dtype=float)
    if moments.shape != (3,) or not np.isfinite(moments).all():
        raise ValueError('Require three finite unbroadened spectral moments')
    scale = np.array([1e-7, 2e-7, 5e-7])
    spacing = energies[1]-energies[0]
    R = np.diff(np.eye(len(energies)), n=2, axis=0)/spacing**2.5
    base = np.vstack([W@K[use], M/scale[:, None]])
    target = np.r_[W@green[use], moments/scale]
    fits = []
    for alpha in alphas:
        w, _ = nnls(np.vstack([base, np.sqrt(alpha)*R]),
            np.r_[target, np.zeros(len(R))], maxiter=40*len(energies))
        pred = K@w
        residual = W@(pred[use]-green[use])
        fits.append(dict(alpha=float(alpha), energies=energies, weights=w,
            predicted=pred, rank=len(W), chi2_per_mode=float(residual@residual/len(W)),
            moment_residuals=M@w-moments, lower_bound=lower,
            upper_support=float(energies[-1]), tau_limit=tau_limit,
            edge_weight=float(w[-3:].sum())))
    return fits


def select(tau, blocks, owners, parameters, alphas, **kwargs):
    green, covariance, _ = estimate_green(blocks, parameters)
    folds = []
    use = tau <= kwargs.get('tau_limit', 8.)
    if len(np.unique(owners)) < 3:
        raise ValueError('At least three independent chains are required for selection')
    for owner in np.unique(owners):
        train, train_cov, _ = estimate_green(blocks[owners != owner], parameters)
        test, test_cov, _ = estimate_green(blocks[owners == owner], parameters)
        W = whitening(test[use], test_cov[np.ix_(use, use)], kwargs.get('tolerance', 1e-8))
        if len(W) < 3:
            raise ValueError('Insufficient held-out covariance rank')
        fits = scan(tau, train, train_cov, parameters, alphas, **kwargs)
        residuals = [W@(fit['predicted'][use]-test[use]) for fit in fits]
        folds.append([float(r@r/len(W)) for r in residuals])
    folds = np.asarray(folds)
    selected = int(np.argmin(folds.mean(axis=0)))
    fits = scan(tau, green, covariance, parameters, alphas, **kwargs)
    return fits[selected], dict(alphas=np.asarray(alphas), folds=folds, selected=selected,
        boundary=selected in [0, len(alphas)-1],
        rule='minimum mean independent held-out-chain score; no VED inputs')
