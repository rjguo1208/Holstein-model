import numpy as np
from scipy.integrate import quad
from continuation import bin_kernel
from resolution import synthetic_spectrum,rebin_gaussian,broaden_truth

def test_synthetic_moments_and_independent_rectangle_integral():
    p=dict(t=1.,omega=1.,g=1.,k=.9*np.pi,mu=-2.8,tau_max=12.,bins=96)
    truth=synthetic_spectrum(p,.6)
    tau=np.array([.0625,2.0625,7.9375,11.9375])
    computed=bin_kernel(tau,.125,truth['energies'],p['mu'])@truth['quadrature']
    expected=[]
    for time in tau:
        val=0.
        for center,width,weight in zip(truth['centers'],truth['widths'],truth['weights']):
            def green(u):
                return quad(lambda energy:np.exp(-(energy-p['mu'])*u),center-width/2,center+width/2)[0]/width
            val+=weight*quad(green,time-.0625,time+.0625)[0]/.125
        expected.append(val)
    np.testing.assert_allclose(computed,expected,rtol=2e-12)
    axis=np.linspace(-1000,1000,100001)
    assert abs(np.trapezoid(broaden_truth(axis,truth,.25),axis)-1)<2e-4

def test_gaussian_time_rebin_covariance_matches_linear_transform():
    rng=np.random.default_rng(88);a=rng.normal(size=(8,8));C=a@a.T
    chains=rng.normal(size=(4,8));p=dict(bins=8,tau_max=2.)
    tau,g,c,pp=rebin_gaussian(chains,C,p,2)
    T=np.kron(np.eye(4),np.ones((1,2))/2)
    np.testing.assert_allclose(g,chains@T.T)
    np.testing.assert_allclose(c,T@C@T.T)
