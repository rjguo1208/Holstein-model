import numpy as np
from numpy.testing import assert_allclose
from continuation import bin_kernel
from map_continuation import scan
from defect_analysis import normalize, reference_integral


def test_real_space_reference_is_zero_hops_and_zero_phonons():
    for origin in [0,3]:
        p=dict(U=1.,mu=-2.7,origin=origin,tau_max=4.,bins=16)
        tau=(np.arange(16)+.5)*.25
        V=-1. if origin==0 else 0.
        G=bin_kernel(tau,.25,[V],p['mu'])[:,0]
        # Expected counts in the pure zero-vertex sector have this exact law.
        total=np.r_[1000,1000,1000*G*.25/reference_integral(p)]
        assert_allclose(normalize(total,p),G,rtol=1e-14)
        assert_allclose(total[2:].sum(),1000,rtol=1e-14)


def test_local_spectral_moments_do_not_use_fixed_momentum_moments():
    p=dict(t=1.,omega=1.,g=1.,U=1.,origin=0,mu=-2.8,tau_max=8.,bins=32)
    tau=(np.arange(32)+.5)*.25
    poles=np.array([-2.5,1.]);weights=np.array([4/7,3/7])
    G=bin_kernel(tau,.25,poles,p['mu'])@weights
    cov=np.diag((G*.001)**2)
    fit=scan(tau,G,cov,p,[0.],energies=np.linspace(-4,6,201),lower_bound=-4.,moments=[1,-1,4])[0]
    assert_allclose(fit['predicted'],G,rtol=2e-5)
    assert_allclose(fit['moment_residuals'],0,atol=1e-8)
    assert fit['lower_bound']==-4.
