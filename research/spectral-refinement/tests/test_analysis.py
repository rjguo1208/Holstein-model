import numpy as np
from scipy.integrate import quad
from numpy.testing import assert_allclose
from analysis import bin_log_factor,estimate_green,reference_integral
from continuation import bin_kernel
from spectral_cv import reconstruct_scan


def test_bin_kernel_against_independent_quadrature():
    tau=np.array([.125,1.125,5.125]);dt=.25;mu=-2.5
    energies=np.array([-3.,-2.5,-2.4999999999,-2.,0.,4.])
    exact=np.array([[quad(lambda x:np.exp(-(e-mu)*x),t-dt/2,t+dt/2)[0]/dt
                     for e in energies] for t in tau])
    assert_allclose(bin_kernel(tau,dt,energies,mu),exact,rtol=2e-13)


def test_finite_bin_correction_preserves_pole_residue():
    E,mu,Z,dt=-2.2,-2.7,.73,.5
    tau=np.arange(.25,15,dt)
    G=Z*bin_kernel(tau,dt,[E],mu)[:,0]
    slope,intercept=np.polyfit(tau,np.log(G),1)
    assert_allclose(mu-slope,E,atol=1e-13)
    assert_allclose(np.exp(intercept-bin_log_factor(E-mu,dt)),Z,rtol=1e-13)


def test_shared_reference_produces_correlated_errors():
    p=dict(t=1,k=0,mu=-2.5,tau_max=2,bins=2)
    rows=np.array([[10,30,0,10,20],[12,30,0,10,20],[8,30,0,10,20],[11,30,0,10,20]],float)
    G,C,_=estimate_green(rows,p)
    assert_allclose(G,reference_integral(p)*np.array([40,80])/41)
    assert C[0,1]>0
    assert_allclose(C[0,1]**2,C[0,0]*C[1,1],rtol=1e-12)


def test_two_pole_synthetic_spectrum_preserves_moments_and_data():
    # Synthetic, exactly known two-pole measure: its first moment is -2.
    # No Holstein VED input is used here.
    E,Z=-2.2,.8;excited=-1.2
    variance=Z*E**2+(1-Z)*excited**2-4
    p=dict(t=1,k=0,g=np.sqrt(variance),mu=-2.8,tau_max=8,bins=32)
    tau=(np.arange(32)+.5)*.25
    G=bin_kernel(tau,.25,[E,excited],p['mu'])@np.array([Z,1-Z])
    C=np.diag((G*1e-3)**2);prior_cov=np.diag([1e-4**2,1e-3**2])
    result=reconstruct_scan(tau,G,C,p,[E,Z],prior_cov,[0.],cutoff=4.075,size=81)
    f=result['fits'][0]
    assert_allclose(f['predicted'],G,rtol=5e-4)
    assert np.max(np.abs(f['moment_residuals']))<1e-5
    assert abs(f['E']-E)<2e-4 and abs(f['Z']-Z)<3e-3
    assert np.all(f['weights']>=0)
