"""Blind synthetic resolution checks with the measured absolute covariance.

The synthetic spectra have exactly the Hamiltonian's first three moments.
They contain a low-energy rectangle and either one narrow peak, a narrow
doublet, or a broad main continuum. They are test cases, not spectral priors.
"""
import numpy as np
from continuation import bin_kernel
from refinement import exact_moments,fit_config,choose_alpha
from spectral_cv import whitening

def synthetic_spectrum(p,separation=0.,main_width=.06):
    mean=exact_moments(p)[1];variance=p['g']**2
    lam=p['g']**2/(2*p['t']*p['omega'])
    low=-2*p['t']-lam*p['t'];low_width=.2*p['t']
    low_variance=low_width**2/12
    main_variance=main_width**2/12+separation**2/4
    D=mean-low
    a=main_variance-low_variance;b=variance-low_variance+D*D
    discriminant=b*b-4*a*D*D
    if discriminant<=0:raise ValueError('Infeasible synthetic spectrum')
    main_weight=2*D*D/(b+np.sqrt(discriminant))
    center=low+D/main_weight
    centers=np.array([low,center-separation/2,center+separation/2])
    widths=np.array([low_width,main_width,main_width])
    weights=np.array([1-main_weight,main_weight/2,main_weight/2])
    if np.any(weights<=0):raise ValueError('Nonpositive synthetic spectral weight')
    lower=-2*p['t']-p['g']**2/p['omega']
    if np.any(centers-widths/2<lower):raise ValueError('Synthetic support below physical lower bound')
    nodes,quad_weights=np.polynomial.legendre.leggauss(24)
    energies=(centers[:,None]+widths[:,None]*nodes/2).reshape(-1)
    quadrature=(weights[:,None]*quad_weights/2).reshape(-1)
    moments=np.array([quadrature.sum(),energies@quadrature,energies**2@quadrature])
    np.testing.assert_allclose(moments,exact_moments(p),atol=2e-12,rtol=0)
    return dict(centers=centers,widths=widths,weights=weights,energies=energies,
                quadrature=quadrature,main_center=float(center),main_weight=float(main_weight),
                separation=separation,main_width=main_width)

def broaden_truth(axis,truth,eta):
    lo=truth['centers']-truth['widths']/2;hi=truth['centers']+truth['widths']/2
    return ((np.arctan((axis[:,None]-lo)/eta)-np.arctan((axis[:,None]-hi)/eta))
            /(np.pi*truth['widths']))@truth['weights']

def rebin_gaussian(chains,cov,p,factor):
    n=p['bins'];new=n//factor
    if n%factor:raise ValueError('Nondivisible synthetic time bins')
    g=chains.reshape(len(chains),new,factor).mean(axis=2)
    c=cov.reshape(new,factor,new,factor).mean(axis=(1,3))
    pp=dict(p,bins=new)
    tau=(np.arange(new)+.5)*p['tau_max']/new
    return tau,g,c,pp

def select_gaussian(chains,cov,p,config):
    # cov is the covariance of the average of all independent synthetic chains.
    n=len(chains)
    tau,values,C,pp=rebin_gaussian(chains,cov,p,config['factor'])
    tt,tv,TC,tp=rebin_gaussian(chains,cov,p,2)
    use=tt<=8.;folds=[]
    for i in range(n):
        train=np.delete(values,i,axis=0).mean(axis=0)
        fits=fit_config(tau,train,C*n/(n-1),pp,config)
        test=tv[i];W=whitening(test[use],(TC*n)[np.ix_(use,use)],1e-8)
        scores=[]
        for fit in fits:
            pred=bin_kernel(tt[use],tp['tau_max']/tp['bins'],fit['energies'],p['mu'])@fit['weights']
            r=W@(pred-test[use])
            scores.append(float(r@r/len(W)) if fit.get('converged',True) else float('inf'))
        folds.append(scores)
    choice=choose_alpha(folds,[f['alpha'] for f in fits])
    fits=fit_config(tau,values.mean(axis=0),C,pp,config)
    selected=fits[choice['selected']]
    if not selected.get('converged',True):raise ValueError('Synthetic fit failed convergence')
    return selected
