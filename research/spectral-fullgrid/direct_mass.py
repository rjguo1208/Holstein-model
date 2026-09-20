"""Curvature from the slope of <K-J^2> at k=0, with a free intercept.

The intercept is -d_k^2 log Z, so dividing a finite-tau measurement by tau
without extrapolation is biased. Uniform-bin averaging only changes the
intercept of the asymptotic straight line. This analysis uses the saved bare-G
diagrams and does not need a reference solution or finite-k energy differences.
"""
from __future__ import annotations
import numpy as np
from analysis import combine_blocks, jackknife_error


def curvature(totals, bins):
    counts=totals[...,3:3+bins]
    # Reflection symmetry makes <J>=0 exactly at k=0. Do not square a noisy
    # sample mean, which would introduce a positive finite-sample bias.
    return -totals[...,3+3*bins:3+4*bins]/counts


def fit_direct_mass(rows, owners, meta, lower=6., upper=20., tolerance=1e-6,
                    time_bin_factor=1,shrinkage=0.):
    if abs(meta['k'])>1e-12:
        raise ValueError('This estimator uses k=0 reflection symmetry')
    n=meta['bins']
    if n%time_bin_factor:raise ValueError('Time rebinning must divide the bin count')
    if time_bin_factor>1:
        fields=rows[:,3:].reshape(len(rows),4,n//time_bin_factor,time_bin_factor).sum(axis=-1)
        rows=np.column_stack([rows[:,:3],fields.reshape(len(rows),-1)])
        n//=time_bin_factor
    tau=(np.arange(n)+.5)*meta['tau_max']/n
    use=(tau>=lower)&(tau<=upper)
    grouped=combine_blocks(rows,owners,4)
    total=grouped.sum(axis=0)
    point=curvature(total,n)
    leave=curvature(total[None,:]-grouped,n)
    centered=leave[:,use]-leave[:,use].mean(axis=0)
    C=(len(grouped)-1)/len(grouped)*centered.T@centered
    C=(1-shrinkage)*C+shrinkage*np.diag(np.diag(C))
    eig,vec=np.linalg.eigh(C)
    precision=(vec/np.maximum(eig,eig[-1]*tolerance))@vec.T
    X=np.column_stack([np.ones(use.sum()),tau[use]])
    operator=np.linalg.solve(X.T@precision@X,X.T@precision)
    beta=operator@point[use]
    def parameters(values):
        b=values[...,use]@operator.T
        return np.stack([b[...,0],b[...,1],2*meta['t']/b[...,1]],axis=-1)
    errors={}
    for factor in [1,2,4]:
        blocks=combine_blocks(rows,owners,factor)
        samples=curvature(blocks.sum(axis=0)[None,:]-blocks,n)
        errors[str(factor)]=jackknife_error(parameters(samples)).tolist()
    error=np.max(list(errors.values()),axis=0)
    chain_error=None
    if len(np.unique(owners))>=4:
        chain_sums=np.array([rows[owners==i].sum(axis=0) for i in np.unique(owners)])
        samples=curvature(chain_sums.sum(axis=0)[None,:]-chain_sums,n)
        chain_error=jackknife_error(parameters(samples))
        error=np.maximum(error,chain_error)
    residual=point[use]-X@beta
    current=total[3+2*n:3+3*n]/total[3:3+n]
    return dict(tau_min=lower,tau_max=upper,intercept=float(beta[0]),
        inverse_mass=float(beta[1]),inverse_mass_error=float(error[1]),
        mass_ratio=float(2*meta['t']/beta[1]),mass_ratio_error=float(error[2]),
        chi2_per_dof=float(residual@precision@residual/(use.sum()-2)),
        reblocking_errors=errors,leave_chain_error=None if chain_error is None else chain_error.tolist(),
        time_bin_factor=time_bin_factor,covariance_shrinkage=shrinkage,covariance_relative_floor=tolerance,
        max_absolute_mean_current=float(abs(current[use]).max()))


def mass_curve(rows,owners,meta):
    blocks=combine_blocks(rows,owners,4); total=blocks.sum(axis=0); n=meta['bins']
    leave=curvature(total[None,:]-blocks,n)
    return curvature(total,n),jackknife_error(leave)
