"""Local, exact-moment MaxEnt; adapted from the published clean-system dual.

The defect potential changes both the support bound and local moments.
Whitened kernel SVD is done once per alpha scan, with smooth-to-sharp warm starts.
No reference spectrum or ground-state pole is used in continuation.
"""
from __future__ import annotations
from functools import lru_cache
import numpy as np
from scipy.special import logsumexp
from continuation import bin_kernel
from spectral_cv import whitening


def local_moments(p):
    V=-p['U'] if p['origin']==0 else 0.
    return np.array([1.,V,V*V+2*p['t']**2+p['g']**2])


@lru_cache(maxsize=96)
def operators(bins,tau_max,mu,lower,upper,size):
    en=np.linspace(lower,upper,size)
    tau=(np.arange(bins)+.5)*tau_max/bins
    K=bin_kernel(tau,tau_max/bins,en,mu)
    M=np.vstack([np.ones(size),en,en**2])
    Q,R=np.linalg.qr(M.T,mode='reduced')
    for a in [en,K,M,Q,R]:a.setflags(write=False)
    return en,K,M,Q,R


def entropy_scan(tau,green,cov,p,alphas,size=641,upper=10.,tau_limit=8.,tolerance=1e-8,
                 target_moments=None,energies=None,max_iterations=400,extended_precision=True,
                 whitening_mode='relative_green'):
    from defect_entropy_stable import long_solve
    lower=-np.hypot(2*p['t'],p['U'])-p['g']**2/p['omega']
    en,K,M,Q,R=operators(p['bins'],p['tau_max'],p['mu'],lower,upper,size)
    if energies is not None:
        en=np.asarray(energies);K=bin_kernel(tau,p['tau_max']/p['bins'],en,p['mu'])
        M=np.vstack([np.ones(len(en)),en,en**2]);Q,R=np.linalg.qr(M.T,mode='reduced')
    use=tau<=tau_limit
    if whitening_mode == 'variance':
        from defect_momentum_mc import variance_whitening
        W=variance_whitening(cov[np.ix_(use,use)],tolerance)
    elif whitening_mode == 'relative_green':
        W=whitening(green[use],cov[np.ix_(use,use)],tolerance)
    else:
        raise ValueError('Unknown covariance whitening mode')
    if len(W)<3:raise ValueError('Insufficient covariance rank')
    moments=local_moments(p) if target_moments is None else np.asarray(target_moments)
    c=np.linalg.solve(R.T,moments);B=W@K[use];y=W@green[use]
    projected=B-(B@Q)@Q.T;yp=y-(B@Q)@c
    U,s,V=np.linalg.svd(projected,full_matrices=False)
    keep=s>s[0]*1e-12;s=s[keep];U=U[:,keep];V=V[keep]
    F=np.asarray(np.vstack([V,Q[:,1:].T]),dtype=np.longdouble)
    d=np.asarray(np.r_[(U.T@yp)/s,c[1:]],dtype=np.longdouble)
    logm=np.full(len(en),-np.log(len(en)),dtype=np.longdouble)
    theta=np.zeros(len(d),dtype=np.longdouble);result=[]
    # Warm starts and centered objective differences avoid singular trial
    # spectra and cancellation at the solution. Damping affects only the
    # Newton linear solve; every line search uses the unchanged objective.
    requested=set(float(a) for a in alphas)
    # Start where entropy dominates the largest data mode, then decrease it
    # gradually. This homotopy is essential for long-time defect kernels,
    # whose rigorous support bound lies well below the physical ground state.
    highest=max(float(max(alphas)),float(s[0]**2))
    bridge=10.**np.arange(np.ceil(np.log10(highest)),np.log10(max(alphas)),-1.)
    pending=sorted(requested|set(bridge),reverse=True)
    last_alpha=None;last_theta=theta.copy();refinements=0
    while pending:
        alpha=pending.pop(0)
        if alpha<=0:raise ValueError('Require positive entropy alpha')
        reg=np.asarray(np.r_[alpha/s**2,0.,0.],dtype=np.longdouble)
        def evaluate(x):
            z=logm-F.T@x;logw=z-logsumexp(z);w=np.exp(logw)
            normalization=w.sum();logw-=np.log(normalization);w/=normalization
            fw=F@w;grad=d-fw+reg*x;centered=F-fw[:,None]
            H=(centered*w)@centered.T+np.diag(reg)
            return grad,w,H,centered,logw
        damped=0
        for iteration in range(max_iterations):
            grad,w,H,centered,logw=evaluate(theta)
            mr=M@np.asarray(w,float)-moments
            if max(abs(grad))<1e-9 and max(abs(mr))<1e-7:break
            scale=np.sqrt(np.maximum(np.diag(H),np.longdouble('1e-30')))
            normalized=H/np.outer(scale,scale);accepted=False
            for damping in [0.,1e-18,1e-16,1e-14,1e-12,1e-10,1e-8,1e-6,1e-4,.01]:
                try:step=long_solve(normalized+damping*np.eye(len(grad)),-grad/scale)/scale
                except np.linalg.LinAlgError:continue
                slope=grad@step
                if not np.isfinite(step).all() or slope>=0:continue
                fraction=np.longdouble(1)
                for _ in range(65):
                    delta=fraction*step;change=-(centered.T@delta)
                    if max(abs(change))<.1:
                        cumulant=np.log1p(np.sum(w*(np.expm1(change)-change)))
                    else:cumulant=logsumexp(logw+change)
                    decrease=grad@delta+cumulant+.5*np.dot(reg*delta,delta)
                    if decrease<=np.longdouble('1e-4')*fraction*slope:
                        candidate=theta+delta
                        if np.array_equal(candidate,theta):break
                        theta=candidate;accepted=True;damped+=damping>0;break
                    fraction*=.5
                if accepted:break
            if not accepted:break
        grad,w,_,_,_=evaluate(theta);w=np.asarray(w,float)
        mr=M@w-moments;pred=K@w;r=W@(pred[use]-green[use])
        converged=bool(max(abs(grad))<1e-8 and max(abs(mr))<1e-6)
        if not converged:
            if last_alpha is None or refinements>=80 or last_alpha/alpha<1.000001:
                raise ValueError(f'Convex dual failed at alpha={alpha}, gradient={max(abs(grad))}, moment={max(abs(mr))}')
            # Never carry a failed Newton iterate into the next alpha. Insert
            # a smaller continuation step and retry the identical objective.
            pending=[float(np.sqrt(last_alpha*alpha)),alpha]+pending
            theta=last_theta.copy();refinements+=1
            continue
        last_alpha=alpha;last_theta=theta.copy()
        if float(alpha) not in requested:continue
        result.append(dict(alpha=float(alpha),energies=en,weights=w,predicted=pred,rank=len(W),
            chi2_per_mode=float(r@r/len(W)),moment_residuals=mr,lower_bound=lower,upper_support=float(en[-1]),
            tau_limit=tau_limit,edge_weight=float(w[-3:].sum()),method='maximum entropy; exact local moments',
            converged=converged,dual_gradient=float(max(abs(grad))),iterations=iteration+1,damped_iterations=damped,homotopy_refinements=refinements,
            numerical_method='warm-started long-double Newton; centered objective differences; adaptive linear-solve damping'))
    return sorted(result,key=lambda f:f['alpha'])
