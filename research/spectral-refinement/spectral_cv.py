"""Spectral reconstruction selected by held-out independent Monte Carlo chains.

The long-time E/Z likelihood comes from a separate simulation. A linearized
ground-pole energy shift propagates that correlated prior in the fit. The
reference spectrum is never an input to reconstruction or model selection.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import nnls
from continuation import bin_kernel,broaden


def whitening(green,cov,tolerance=1e-8):
    relative=cov/np.outer(green,green)
    values,vectors=np.linalg.eigh(relative)
    use=values>values[-1]*tolerance
    return (vectors[:,use]/np.sqrt(values[use])).T/green[None,:]


def reconstruct_scan(tau,green,cov,parameters,prior,prior_cov,
                     alphas,cutoff=12.,size=180,tolerance=1e-8):
    E,Z=prior; sigma=np.sqrt(prior_cov[0,0])
    dt=parameters['tau_max']/parameters['bins'];mu=parameters['mu']
    omega=E+np.linspace(.075,cutoff,size)
    energies=np.r_[E,omega]
    kernel=bin_kernel(tau,dt,energies,mu)
    derivative=Z*(bin_kernel(tau,dt,[E+1e-5],mu)[:,0]-bin_kernel(tau,dt,[E-1e-5],mu)[:,0])/2e-5
    # x[-1]>=0 encodes delta_E=sigma*(x[-1]-5). Five sigma is
    # a numerical bound only; flag solutions approaching it.
    K=np.column_stack([kernel,sigma*derivative])
    target=green+5*sigma*derivative
    W=whitening(green,cov,tolerance)
    eps=-2*parameters['t']*np.cos(parameters['k'])
    moments=np.array([1,eps,eps*eps+parameters['g']**2])
    M=np.column_stack([np.vstack([np.ones(size+1),energies,energies**2]),
                       sigma*np.array([0,Z,2*Z*E])])
    moment_sigma=np.array([1e-6,2e-6,5e-6])
    P=np.zeros((2,size+2));P[0,-1]=sigma;P[1,0]=1
    pe,pu=np.linalg.eigh(prior_cov)
    PW=(pu/np.sqrt(pe)).T
    spacing=omega[1]-omega[0]
    regularizer=np.zeros((size-2,size+2))
    regularizer[:,1:-1]=np.diff(np.eye(size),n=2,axis=0)/spacing**2.5
    A=np.vstack([W@K,M/moment_sigma[:,None],PW@P])
    b=np.r_[W@target,(moments+5*M[:,-1])/moment_sigma,PW@np.array([5*sigma,Z])]
    fits=[]
    for alpha in alphas:
        x,_=nnls(np.vstack([A,np.sqrt(alpha)*regularizer]),
                 np.r_[b,np.zeros(size-2)],maxiter=30*(size+2))
        e=E+sigma*(x[-1]-5);z=x[0];w=x[1:-1]
        predicted=bin_kernel(tau,dt,[e],mu)[:,0]*z+kernel[:,1:]@w
        residual=W@(predicted-green)
        exact_moments=np.array([z+w.sum(),z*e+w@omega,z*e**2+w@omega**2])
        fits.append(dict(alpha=float(alpha),E=float(e),Z=float(z),weights=w,
                         predicted=predicted,chi2_per_mode=float(residual@residual/len(W)),
                         rank=len(W),prior_shift_sigma=float((e-E)/sigma),
                         moment_residuals=(exact_moments-moments).tolist(),
                         linearization_error=float(np.max(np.abs(predicted-(K@x-5*sigma*derivative))))))
    return dict(omega=omega,fits=fits,tolerance=tolerance,size=size,cutoff=cutoff)


def evaluate(tau,green,cov,parameters,omega,fit,tolerance=1e-8):
    dt=parameters['tau_max']/parameters['bins'];mu=parameters['mu']
    prediction=bin_kernel(tau,dt,[fit['E']],mu)[:,0]*fit['Z']+bin_kernel(tau,dt,omega,mu)@fit['weights']
    W=whitening(green,cov,tolerance)
    residual=W@(prediction-green)
    return float(residual@residual/len(W))


def prior_covariance(meta,saved):
    # Correlation from the long-time jackknife, marginal errors from the
    # largest of the reported block sizes.
    from continuation import perturb_ez
    primary=meta['primary'];p=meta['parameters'];E,Z=primary['E'],primary['Z']
    tau,G,C,leave=(saved[k] for k in ['tau','G','covariance','jackknife'])
    samples=np.array([perturb_ez(j-G,G,C,tau,p['tau_max']/p['bins'],p['mu'],E,Z) for j in leave])
    correlation=np.corrcoef(samples.T)[0,1]
    e,z=primary['E_error'],primary['Z_error']
    return np.array([[e*e,correlation*e*z],[correlation*e*z,z*z]])
