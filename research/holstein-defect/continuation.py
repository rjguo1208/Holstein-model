"""Positive, moment-constrained spectral reconstruction with explicit sensitivity.

This is a regularized inverse Laplace transform, not an exact real-frequency
solver. Spectral features must be judged from noise and regularization scans.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import nnls


def bin_kernel(tau,dt,energy,mu):
    rate=np.atleast_1d(energy)-mu
    x=rate*dt
    factor=np.ones_like(rate)
    nonzero=np.abs(x)>1e-10
    factor[nonzero]=-np.expm1(-x[nonzero])/x[nonzero]
    factor[~nonzero]=1-x[~nonzero]/2+x[~nonzero]**2/6
    return np.exp(-np.outer(tau-dt/2,rate))*factor


def reconstruct(tau,green,cov,parameters,E,Z,cutoff=12.,size=180,alphas=None):
    dt=parameters['tau_max']/parameters['bins'];mu=parameters['mu']
    use=tau<=10.
    values,vectors=np.linalg.eigh(cov[np.ix_(use,use)])
    whiten=(vectors/np.sqrt(np.maximum(values,values[-1]*1e-6))).T
    omega=E+np.linspace(.075,cutoff,size)
    K=bin_kernel(tau[use],dt,omega,mu)
    ground=bin_kernel(tau[use],dt,[E],mu)[:,0]
    target=green[use]-Z*ground
    eps=-2*parameters['t']*np.cos(parameters['k'])
    # Exact unbroadened moments of the Hamiltonian, not fitted to VED data.
    moments=np.array([1-Z,eps-Z*E,eps*eps+parameters['g']**2-Z*E*E])
    moments_matrix=np.vstack([np.ones(size),omega,omega**2])
    moment_sigma=np.array([1e-5,2e-5,5e-5])
    # A second-difference penalty on the density; same continuum scaling as
    # the grid is refined (weights = density * domega).
    spacing=omega[1]-omega[0]
    regularizer=np.diff(np.eye(size),n=2,axis=0)/spacing**2.5
    A=np.vstack([whiten@K,moments_matrix/moment_sigma[:,None]])
    b=np.r_[whiten@target,moments/moment_sigma]
    alphas=np.logspace(-5,5,17) if alphas is None else np.asarray(alphas)
    fits=[]
    for alpha in alphas:
        aug=np.vstack([A,np.sqrt(alpha)*regularizer])
        weights,_=nnls(aug,np.r_[b,np.zeros(size-2)],maxiter=20*size)
        residual=whiten@(K@weights-target)
        fits.append(dict(alpha=float(alpha),chi2_per_bin=float(residual@residual/use.sum()),
                         moment_residuals=(moments_matrix@weights-moments).tolist(),weights=weights))
    valid=[i for i,f in enumerate(fits) if f['chi2_per_bin']<=1.2]
    selected=max(valid) if valid else int(np.argmin([f['chi2_per_bin'] for f in fits]))
    return dict(omega=omega,selected=selected,fits=fits,E=float(E),Z=float(Z),
                discrepancy_met=bool(valid),size=size,cutoff=cutoff)


def broaden(axis,omega,weights,E,Z,eta):
    continuum=(eta/np.pi/((axis[:,None]-omega[None,:])**2+eta**2))@weights
    pole=Z*eta/np.pi/((axis-E)**2+eta**2)
    return pole+continuum,continuum


def perturb_ez(delta,green,cov,tau,dt,mu,E,Z):
    use=(tau>=6)&(tau<=20)
    X=np.column_stack([np.ones(use.sum()),tau[use]])
    C=cov[np.ix_(use,use)]/np.outer(green[use],green[use])
    precision=np.linalg.pinv(C,rcond=1e-6)
    estimator=np.linalg.solve(X.T@precision@X,X.T@precision)
    beta=estimator@(delta[use]/green[use])
    newE=E-beta[1]
    def logfactor(e):
        x=(e-mu)*dt/2
        return np.log(np.sinh(x)/x) if abs(x)>1e-8 else x*x/6
    newZ=Z*np.exp(beta[0]-logfactor(newE)+logfactor(E))
    return newE,newZ
