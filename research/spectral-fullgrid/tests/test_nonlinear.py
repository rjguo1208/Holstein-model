import numpy as np
from scipy.integrate import quad
from numpy.testing import assert_allclose
from continuation import bin_kernel
from spectral_nonlinear import SpectralProblem,rectangle_kernel,broaden_fit,stochastic_rectangles


def test_rectangle_kernel_against_nested_quadrature():
    tau=np.array([.0625,.1875,1.0625,8.0625]);dt=.125;mu=-2.8
    c=np.array([-2.,.5,6.]);w=np.array([.02,1.7,3.])
    expected=np.array([[quad(lambda t:quad(lambda e:np.exp(-(e-mu)*t),a-b/2,a+b/2,
                          epsabs=1e-13)[0]/b,t0-dt/2,t0+dt/2,epsabs=1e-13)[0]/dt
                        for a,b in zip(c,w)] for t0 in tau])
    assert_allclose(rectangle_kernel(tau,dt,c,w,mu),expected,rtol=3e-11,atol=1e-13)


def test_exact_profile_propagates_pole_and_exact_moments():
    E,Z,excited=-2.2,.8,-1.2
    variance=Z*E**2+(1-Z)*excited**2-4
    p=dict(t=1,k=0,g=np.sqrt(variance),mu=-2.8,tau_max=8,bins=32)
    tau=(np.arange(32)+.5)*.25
    G=bin_kernel(tau,.25,[E,excited],p['mu'])@np.array([Z,1-Z])
    C=np.diag((G*1e-4)**2)
    problem=SpectralProblem(tau,G,C,p,[E+.0005,Z-.001],np.diag([.001**2,.005**2]))
    fit=problem.profile(np.array([excited]))
    assert_allclose(fit['predicted'],G,rtol=2e-5)
    assert abs(fit['E']-E)<1e-5
    assert max(abs(np.array(fit['moment_residuals'])))<1e-6
    assert fit['weights'].min()>=0


def test_rectangle_broadening_preserves_unit_area_and_point_limit():
    axis=np.linspace(-10000,10000,200001)
    fit=dict(E=-2.,Z=.6,centers=np.array([1.]),widths=np.array([.2]),weights=np.array([.4]))
    full,rest=broaden_fit(axis,fit,.25)
    assert_allclose(np.trapezoid(full,axis),1,atol=3e-5)
    assert_allclose(np.trapezoid(rest,axis),.4,atol=1e-5)


def test_stochastic_rectangle_search_respects_free_electron_limit():
    p=dict(t=1,k=0,g=0,mu=-2.8,tau_max=4,bins=16)
    tau=(np.arange(16)+.5)*.25
    G=bin_kernel(tau,.25,[-2.],p['mu'])[:,0]
    problem=SpectralProblem(tau,G,np.diag((G*1e-3)**2),p,
                            [-2.,1.],np.diag([1e-3**2,1e-3**2]))
    ensemble=stochastic_rectangles(problem,seed=11,solutions=2,updates=15,rectangles=6)
    assert_allclose(ensemble['predicted'],G,rtol=1e-6)
    assert_allclose(ensemble['Z'],1.,atol=1e-7)
    assert max(abs(ensemble['moments_residual']))<1e-7
