"""Covariance-aware continuation comparisons, with an exact-moment MaxEnt dual.

MaxEnt minimizes ||W(Kw-G)||^2/2 + alpha*sum(w*log(w/m)), w>=0,
subject to the three Hamiltonian moments. The default is a uniform density.
The small convex dual uses normalized exponential weights and exact moment
multipliers. Projection removes the moment directions before the kernel SVD;
it does not change the moment-constrained least-squares objective.
No reference spectra enter this module.
"""
from __future__ import annotations
import numpy as np
from scipy.special import logsumexp
from analysis import estimate_green, combine_blocks
from continuation import bin_kernel
from map_continuation import rebin_conditional, scan
from spectral_cv import whitening

TIKHONOV_ALPHAS=np.r_[0.,np.logspace(-9,3,13)]
MAXENT_ALPHAS=np.logspace(-3,4,8)

def exact_moments(p):
    eps=-2*p['t']*np.cos(p['k'])
    return np.array([1.,eps,eps*eps+p['g']**2])

def entropy_scan(tau,green,cov,p,alphas,size=281,upper=10.,tau_limit=8.,tolerance=1e-8,
                 target_moments=None,energies=None,max_iterations=250,extended_precision=False):
    lower=-2*p['t']-p['g']**2/p['omega']
    en=np.linspace(lower,upper,size) if energies is None else np.asarray(energies)
    use=tau<=tau_limit
    W=whitening(green[use],cov[np.ix_(use,use)],tolerance)
    K=bin_kernel(tau,p['tau_max']/p['bins'],en,p['mu'])
    B=W@K[use];y=W@green[use]
    M=np.vstack([np.ones(len(en)),en,en**2])
    moments=exact_moments(p) if target_moments is None else np.asarray(target_moments)
    if not np.isclose(moments[0],1.):raise ValueError('MaxEnt requires normalized moments')
    Q,R=np.linalg.qr(M.T,mode='reduced')
    c=np.linalg.solve(R.T,moments)
    projected=B-(B@Q)@Q.T
    yp=y-(B@Q)@c
    U,s,V=np.linalg.svd(projected,full_matrices=False)
    keep=s>s[0]*1e-12
    s=s[keep];U=U[:,keep];V=V[keep]
    F=np.vstack([V,Q[:,1:].T]);d=np.r_[(U.T@yp)/s,c[1:]]
    logm=np.full(len(en),-np.log(len(en)))
    theta=np.zeros(len(d));fits=[]
    if extended_precision:
        F=F.astype(np.longdouble);d=d.astype(np.longdouble)
        logm=logm.astype(np.longdouble);theta=theta.astype(np.longdouble)
    # Warm starts move from the smoothest spectrum toward the least penalized.
    for alpha in sorted(alphas,reverse=True):
        if alpha<=0:raise ValueError('Entropy alpha must be positive')
        reg=np.r_[alpha/s**2,0.,0.]
        if extended_precision:reg=reg.astype(np.longdouble)
        def evaluate(x,hessian=False):
            z=logm-F.T@x;logz=logsumexp(z);w=np.exp(z-logz)
            f=logz+d@x+.5*np.dot(reg*x,x)
            fw=F@w;grad=d-fw+reg*x
            if not hessian:return f,grad,w
            centered=F-fw[:,None]
            H=(centered*w)@centered.T+np.diag(reg)
            return f,grad,w,H
        converged=False
        for iteration in range(max_iterations):
            f,grad,w,H=evaluate(theta,True)
            if np.max(abs(grad))<2e-11:
                converged=True;break
            scale=np.sqrt(np.maximum(np.diag(H),1e-30))
            eig,vec=np.linalg.eigh(np.asarray(H/np.outer(scale,scale),dtype=float))
            step=-((vec/np.maximum(eig,max(eig[-1]*1e-13,1e-18)))@(vec.T@(grad/scale)))/scale
            slope=float(grad@step)
            fraction=1.
            for _ in range(50):
                trial=theta+fraction*step
                ft,gt,wt=evaluate(trial)
                if ft<=f+1e-4*fraction*slope+1e-14:
                    theta=trial;break
                fraction*=.5
            else:break
            if fraction*np.linalg.norm(step)<1e-12:break
        f,grad,w=evaluate(theta)
        w=np.asarray(w,dtype=float)
        mr=M@w-moments
        converged=bool(converged or (np.max(abs(grad))<1e-8 and np.max(abs(mr))<1e-6))
        predicted=K@w;r=W@(predicted[use]-green[use])
        fits.append(dict(alpha=float(alpha),energies=en,weights=w,predicted=predicted,
            rank=len(W),chi2_per_mode=float(r@r/len(W)),moment_residuals=mr,
            lower_bound=lower,upper_support=float(en[-1]),tau_limit=tau_limit,
            edge_weight=float(w[-3:].sum()),method='maximum entropy; exact moments',
            converged=converged,dual_gradient=float(np.max(abs(grad))),iterations=iteration+1))
    return sorted(fits,key=lambda f:f['alpha'])

def candidates():
    base=dict(method='tikhonov',factor=2,block=4,size=281,upper=10.,tau_limit=8.,tolerance=1e-8)
    variants=[('baseline',{}),('time_fine',dict(factor=1)),('time_short',dict(tau_limit=4.)),
        ('time_long',dict(tau_limit=12.)),('cov_cut_1e-6',dict(tolerance=1e-6)),
        ('cov_cut_1e-10',dict(tolerance=1e-10)),('energy_141',dict(size=141)),
        ('energy_561',dict(size=561)),('energy_shift',dict(shift=.5)),
        ('upper_6',dict(upper=6.)),('upper_14',dict(upper=14.)),
        ('block_8',dict(block=8)),('block_16',dict(block=16)),
        ('maxent',dict(method='maxent'))]
    return [dict(base,name=name,**{})|change for name,change in variants]

def prepare(rows,owners,p,config):
    grouped=combine_blocks(rows,owners,config['block'])
    grouped_owners=np.concatenate([np.full(np.sum(owners==owner)//config['block'],owner)
                                  for owner in np.unique(owners)])
    tau,grouped,pp=rebin_conditional(grouped,p,config['factor'])
    return tau,grouped,grouped_owners,pp

def fit_config(tau,G,C,p,config):
    options={k:config[k] for k in ['size','upper','tau_limit','tolerance']}
    if config.get('shift'):
        lower=-2*p['t']-p['g']**2/p['omega']
        en=np.linspace(lower,config['upper'],config['size'])
        options['energies']=en+config['shift']*(en[1]-en[0])
    if config['method']=='maxent':return entropy_scan(tau,G,C,p,MAXENT_ALPHAS,**options)
    return scan(tau,G,C,p,TIKHONOV_ALPHAS,**options)

def choose_alpha(folds,alphas):
    """One-standard-error heuristic within one parameterization, favoring smoothing.

    Folds share training observations; this is a model-selection heuristic,
    not a calibrated confidence interval or a resolution guarantee.
    """
    a=np.asarray(folds);means=np.mean(a,axis=0)
    if not np.isfinite(means).any():raise ValueError('No candidate converged in every fold')
    best=int(np.argmin(means));se=float(np.std(a[:,best],ddof=1)/np.sqrt(len(a)))
    eligible=np.flatnonzero(means<=means[best]+se)
    chosen=int(eligible[-1])
    return dict(best=best,selected=chosen,mean=means,se=se,threshold=float(means[best]+se),
                alphas=np.asarray(alphas),score=float(means[chosen]),
                rule='largest alpha within one fold-standard-error of minimum; heuristic')

def cross_validate(rows,owners,p,config,recover=False):
    """All configurations use identical test bins, cutoff and test covariance."""
    tau,blocks,who,pp=prepare(rows,owners,p,config)
    test_config=dict(config,factor=2,block=4,tau_limit=8.,tolerance=1e-8)
    tt,tb,tw,tp=prepare(rows,owners,p,test_config)
    use=tt<=8.;folds=[];alphas=None
    for owner in np.unique(owners):
        G,C,_=estimate_green(blocks[who!=owner],pp)
        test,tc,_=estimate_green(tb[tw==owner],tp)
        W=whitening(test[use],tc[np.ix_(use,use)],1e-8)
        fits=fit_config(tau,G,C,pp,config)
        alphas=[f['alpha'] for f in fits]
        scores=[]
        for fit in fits:
            pred=bin_kernel(tt[use],tp['tau_max']/tp['bins'],fit['energies'],p['mu'])@fit['weights']
            residual=W@(pred-test[use])
            scores.append(float(residual@residual/len(W)) if fit.get('converged',True) else float('inf'))
        folds.append(scores)
    choice=choose_alpha(folds,alphas)
    G,C,_=estimate_green(blocks,pp)
    fits=fit_config(tau,G,C,pp,config)
    selected=fits[choice['selected']]
    if recover and not selected.get('converged',True) and config['method']=='maxent':
        original=dict(alpha=selected['alpha'],iterations=selected['iterations'],
                      dual_gradient=selected['dual_gradient'],moment_residuals=selected['moment_residuals'])
        options={k:config[k] for k in ['size','upper','tau_limit','tolerance']}
        selected=entropy_scan(tau,G,C,pp,[selected['alpha']],max_iterations=4000,extended_precision=True,**options)[0]
        selected['numerical_recovery']=dict(original=original,method='same objective and selected alpha; fresh start, extended-precision dual accumulation, up to 4000 Newton iterations')
    if not selected.get('converged',True):raise ValueError('Selected MaxEnt solution did not converge')
    return selected,dict(choice,folds=np.asarray(folds),config=config,minimum_fit=fits[choice['best']])
