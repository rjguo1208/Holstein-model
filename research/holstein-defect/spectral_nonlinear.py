"""Exact ground-pole profiling and a continuous-rectangle SOM variant.

The SOM implementation follows the positive random-solution ensemble idea of
Mishchenko et al. (2000), Appendix B. It uses variable projection for nonnegative
areas, a covariance-weighted objective, and an independent E/Z likelihood. It is
not a verbatim implementation of that paper's seven rectangular updates.
No function in this module reads a reference spectrum.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import nnls,minimize_scalar
from scipy.special import expi
from continuation import bin_kernel,broaden
from spectral_cv import whitening


def rectangle_kernel(tau,dt,centers,widths,mu):
    """Exact bin average of the Laplace transform of unit-area rectangles.

All supports must lie above mu. The first time bin has a logarithmic limiting
term rather than Ei(0); narrow rectangles use their point-mass limit.
"""
    c=np.atleast_1d(centers);w=np.atleast_1d(widths)
    lo=c-w/2-mu;hi=c+w/2-mu
    if np.any(lo<=0) or np.any(w<=0):raise ValueError('Rectangle support must be above mu')
    a=np.maximum(tau-dt/2,0)[:,None];b=(tau+dt/2)[:,None]
    # Avoid evaluating Ei(0)-Ei(0) before selecting the analytic limit.
    aa=np.where(a==0,1.,a)
    first=expi(-aa*hi)-expi(-aa*lo)
    first=np.where(a==0,np.log(hi/lo),first)
    second=expi(-b*hi)-expi(-b*lo)
    return (first-second)/(dt*w)


class SpectralProblem:
    def __init__(self,tau,G,C,parameters,prior,prior_cov,tolerance=1e-8):
        self.tau=np.asarray(tau);self.G=np.asarray(G);self.p=parameters
        self.dt=parameters['tau_max']/parameters['bins'];self.prior=np.asarray(prior)
        self.prior_cov=np.asarray(prior_cov);self.sigma=np.sqrt(prior_cov[0,0])
        self.W=whitening(G,C,tolerance)
        self.PW=np.linalg.inv(np.linalg.cholesky(prior_cov))
        eps=-2*parameters['t']*np.cos(parameters['k'])
        self.moments=np.array([1.,eps,eps*eps+parameters['g']**2])
        self.moment_sigma=np.array([1e-7,2e-7,5e-7])

    def solve(self,E,centers,widths=None,alpha=0.,regularizer=None):
        c=np.asarray(centers)
        kernel=(bin_kernel(self.tau,self.dt,c,self.p['mu']) if widths is None else
                rectangle_kernel(self.tau,self.dt,c,widths,self.p['mu']))
        kernel=np.column_stack([bin_kernel(self.tau,self.dt,[E],self.p['mu'])[:,0],kernel])
        energies=np.r_[E,c]
        second=energies**2
        if widths is not None:second[1:]+=np.asarray(widths)**2/12
        moments=np.vstack([np.ones(len(energies)),energies,second])
        prior_A=np.zeros((2,len(energies)));prior_A[:,0]=self.PW[:,1]
        prior_b=self.PW@np.array([self.prior[0]-E,self.prior[1]])
        A=np.vstack([self.W@kernel,moments/self.moment_sigma[:,None],prior_A])
        b=np.r_[self.W@self.G,self.moments/self.moment_sigma,prior_b]
        if alpha and regularizer is not None:
            R=np.pad(regularizer,((0,0),(1,0)))
            A=np.vstack([A,np.sqrt(alpha)*R]);b=np.r_[b,np.zeros(len(R))]
        x,_=nnls(A,b,maxiter=100*len(energies))
        prediction=kernel@x
        residual=self.W@(prediction-self.G)
        prior_residual=self.PW@(np.array([E,x[0]])-self.prior)
        moment_residual=moments@x-self.moments
        objective=float(residual@residual+prior_residual@prior_residual+
                        np.sum((moment_residual/self.moment_sigma)**2))
        if alpha and regularizer is not None:objective+=float(alpha*np.sum((regularizer@x[1:])**2))
        return dict(E=float(E),Z=float(x[0]),centers=c.copy(),
            widths=None if widths is None else np.asarray(widths).copy(),weights=x[1:],
            predicted=prediction,objective=objective,rank=len(self.W),
            chi2_per_mode=float(residual@residual/len(self.W)),
            prior_shift_sigma=float((E-self.prior[0])/self.sigma),
            prior_mahalanobis=float(prior_residual@prior_residual),
            moment_residuals=moment_residual.tolist())

    def profile(self,centers,widths=None,alpha=0.,regularizer=None):
        cache={}
        def evaluate(shift):
            if shift not in cache:
                cache[shift]=self.solve(self.prior[0]+self.sigma*shift,
                                       centers,widths,alpha,regularizer)
            return cache[shift]['objective']
        result=minimize_scalar(evaluate,bounds=(-6.,6.),method='bounded',
                               options={'xatol':.005,'maxiter':36})
        evaluate(0.)
        fit=min(cache.values(),key=lambda f:f['objective'])
        fit['energy_profile_success']=bool(result.success)
        fit['energy_bound_hit']=bool(abs(fit['prior_shift_sigma'])>5.8)
        return fit


def tikhonov_scan(problem,alphas,cutoff=12.,size=180):
    centers=problem.prior[0]+np.linspace(.075,cutoff,size)
    spacing=centers[1]-centers[0]
    R=np.diff(np.eye(size),n=2,axis=0)/spacing**2.5
    fits=[]
    for alpha in alphas:
        fit=problem.profile(centers,alpha=alpha,regularizer=R)
        fit['alpha']=float(alpha);fits.append(fit)
    return fits


def stochastic_rectangles(problem,seed=1234,solutions=16,updates=250,
                          rectangles=14,cutoff=12.):
    """Independent annealed searches, followed by an objective-window average.

    Areas are reoptimized by NNLS after each center/width, birth or death move.
    Hence these are stochastic optimizations, not posterior Monte Carlo draws.
    The exact E profile is optimized initially, periodically, and finally.
    """
    rng=np.random.default_rng(seed)
    lower=problem.prior[0]+.075;upper=problem.prior[0]+cutoff
    finals=[];accepted=0;attempted=0
    for run in range(solutions):
        widths=np.exp(rng.uniform(np.log(.04),np.log(2.),rectangles))
        centers=lower+widths/2+(upper-lower-widths)*rng.beta(1.,2.,rectangles)
        current=problem.profile(centers,widths);best=current
        for step in range(updates):
            centers=current['centers'].copy();widths=current['widths'].copy()
            move=rng.integers(5);i=rng.integers(len(centers))
            if move==0 and len(centers)<2*rectangles:
                width=np.exp(rng.uniform(np.log(.025),np.log(2.)))
                center=rng.uniform(lower+width/2,upper-width/2)
                centers=np.r_[centers,center];widths=np.r_[widths,width]
            elif move==1 and len(centers)>4:
                centers=np.delete(centers,i);widths=np.delete(widths,i)
            elif move==2:
                widths[i]=np.clip(widths[i]*np.exp(rng.normal(0,.7)),.015,min(3.,upper-lower))
            else:
                centers[i]+=rng.normal(0,.15 if move==3 else 1.)
            centers=np.clip(centers,lower+widths/2,upper-widths/2)
            trial=problem.solve(current['E'],centers,widths)
            attempted+=1
            temperature=max(.03,2.*(1-step/max(updates,1))**2)
            delta=trial['objective']-current['objective']
            if delta<0 or rng.random()<np.exp(-min(700.,max(0.,delta)/temperature)):
                current=trial;accepted+=1
                if current['objective']<best['objective']:best=current
            if (step+1)%80==0:
                current=problem.profile(current['centers'],current['widths'])
                if current['objective']<best['objective']:best=current
        best=problem.profile(best['centers'],best['widths'])
        finals.append(best)
    objectives=np.array([f['objective'] for f in finals])
    # The window depends only on data errors; reference spectra are absent.
    keep=objectives<=objectives.min()+max(2.,np.sqrt(2*len(problem.W)))
    chosen=[f for f,k in zip(finals,keep) if k]
    return dict(method='continuous-rectangle stochastic optimization with variable projection',
        solutions=chosen,all_objectives=objectives,retained=int(keep.sum()),
        attempted=attempted,accepted=accepted,seed=seed,updates=updates,
        predicted=np.mean([f['predicted'] for f in chosen],axis=0),
        E=float(np.mean([f['E'] for f in chosen])),Z=float(np.mean([f['Z'] for f in chosen])),
        moments_residual=np.mean([f['moment_residuals'] for f in chosen],axis=0))


def broaden_fit(axis,fit,eta):
    if fit.get('widths') is None:
        return broaden(axis,fit['centers'],fit['weights'],fit['E'],fit['Z'],eta)
    widths=fit['widths'];centers=fit['centers'];weights=fit['weights']
    upper=centers+widths/2;lower=centers-widths/2
    rest=((np.arctan((axis[:,None]-lower)/eta)-np.arctan((axis[:,None]-upper)/eta))
          /(np.pi*widths))@weights
    pole=fit['Z']*eta/np.pi/((axis-fit['E'])**2+eta**2)
    return pole+rest,rest


def broaden_ensemble(axis,ensemble,eta):
    curves=np.asarray([broaden_fit(axis,f,eta) for f in ensemble['solutions']])
    return curves.mean(axis=0)
