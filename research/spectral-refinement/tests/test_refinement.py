"""Independent checks of entropy inference and held-out selection."""
import numpy as np
import pytest
from scipy.optimize import minimize
from continuation import bin_kernel
from refinement import entropy_scan, choose_alpha

@pytest.mark.parametrize('extended',[False,True])
def test_entropy_matches_independent_constrained_primal(extended):
    p=dict(t=1.,omega=1.,g=1.,k=1.,mu=-2.8,tau_max=4.,bins=16)
    tau=(np.arange(16)+.5)*.25
    en=np.linspace(-2,4,15)
    truth=np.exp(-.5*((en-.3)/.7)**2);truth/=truth.sum()
    M=np.vstack([np.ones(len(en)),en,en**2]);mom=M@truth
    K=bin_kernel(tau,.25,en,p['mu']);G=K@truth
    C=np.diag((.015*G)**2);alpha=3.
    fit=entropy_scan(tau,G,C,p,[alpha],energies=en,target_moments=mom,extended_precision=extended)[0]
    assert fit['converged']
    np.testing.assert_allclose(M@fit['weights'],mom,atol=1e-7)
    W=np.diag(1/np.sqrt(np.diag(C)));B=W@K;y=W@G
    def fun(w):
        r=B@w-y
        return .5*r@r+alpha*np.sum(w*np.log(np.maximum(w,1e-300)*len(w)))
    def jac(w):return B.T@(B@w-y)+alpha*(np.log(np.maximum(w,1e-300)*len(w))+1)
    result=minimize(fun,truth,jac=jac,method='SLSQP',bounds=[(1e-14,1)]*len(en),
        constraints={'type':'eq','fun':lambda w:M@w-mom,'jac':lambda w:M},
        options={'maxiter':500,'ftol':1e-10})
    assert result.success,result.message
    np.testing.assert_allclose(fun(fit['weights']),result.fun,rtol=1e-7,atol=1e-7)

def test_one_standard_error_favors_larger_penalty_but_excludes_bad_fit():
    result=choose_alpha([[1.,1.1,9.],[2.,2.1,9.],[1.,1.1,9.],[2.,2.1,9.]],[0.,1.,10.])
    assert result['best']==0 and result['selected']==1
