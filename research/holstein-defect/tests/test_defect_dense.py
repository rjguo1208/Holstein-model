import numpy as np
from numpy.testing import assert_allclose
from continuation import bin_kernel
from defect_entropy import entropy_scan,local_moments
from defect_entropy_stable import solve
from defect_dense import prepare,estimate
from defect_analysis import reference_integral


def test_entropy_known_analytic_optimum_and_stable_solver():
    p=dict(t=1.,omega=1.,g=np.sqrt(.5),U=1.,origin=0,tau_max=6.,bins=24,mu=-3.)
    en=np.linspace(-3.5,5.,121);tau=(np.arange(24)+.5)*.25
    weights=np.exp(-.5*(en-.2)**2);weights/=weights.sum()
    moments=np.array([1.,en@weights,en**2@weights])
    G=bin_kernel(tau,.25,en,p['mu'])@weights;C=np.diag((G*.01)**2)
    for f in entropy_scan(tau,G,C,p,[.1,10.,1000.],energies=en,target_moments=moments):
        assert f['converged'];assert_allclose(f['lower_bound'],-np.sqrt(5)-.5)
        assert_allclose(f['weights'],weights,atol=2e-7)
        assert_allclose(f['moment_residuals'],0,atol=1e-6)
    f=solve(tau,G,C,p,10.,energies=en,target_moments=moments)
    assert f['converged'];assert_allclose(f['weights'],weights,atol=2e-7)
    assert_allclose(local_moments(p),[1.,-1.,3.5])
    assert_allclose(local_moments(dict(p,origin=4)),[1.,0.,2.5])


def test_compact_time_rebin_keeps_exact_bin_integrals_and_chain_boundaries():
    p=dict(U=1.,mu=-2.7,origin=0,tau_max=12.,bins=192)
    tau=(np.arange(192)+.5)*12/192;G=bin_kernel(tau,12/192,[-1.],p['mu'])[:,0]
    one=np.r_[100.,100.,100*G*(12/192)/reference_integral(p)]
    owners=np.repeat(np.arange(8),16);rows=np.tile(one,(128,1))
    for factor in [1,2,4]:
        tt,b,who,pp=prepare(rows,owners,p,dict(factor=factor,block=4))
        actual,C=estimate(b,pp)
        assert_allclose(actual,bin_kernel(tt,pp['tau_max']/pp['bins'],[-1.],p['mu'])[:,0],rtol=2e-14)
        assert_allclose(C,0,atol=1e-26)
        assert_allclose(np.bincount(who),4)


def test_long_time_local_covariance_regression():
    import json
    from pathlib import Path
    folder=Path(__file__).parent/'data'
    data=np.load(folder/'defect_long_time_covariance.npz')
    p=json.loads((folder/'defect_long_time_covariance.json').read_text())['parameters']
    for f in entropy_scan(data['tau'],data['G'],data['covariance'],p,np.logspace(-3,4,8),size=641,tau_limit=12.):
        assert f['converged'] and f['dual_gradient']<1e-8
        assert_allclose(f['moment_residuals'],0,atol=1e-6)
        assert np.min(f['weights'])>=0
