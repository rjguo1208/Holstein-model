import numpy as np
from numpy.testing import assert_allclose
from direct_mass import fit_direct_mass


def test_direct_mass_retains_residue_curvature_intercept():
    rng=np.random.default_rng(77)
    n=32; blocks=160
    p=dict(k=0,t=1,bins=n,tau_max=16)
    tau=(np.arange(n)+.5)*.5
    counts=rng.uniform(800,1200,(blocks,n))
    # A nonzero residue curvature: naive <K-J^2>/tau is NOT 1/m.
    slope=1.4;intercept=.9
    deviations=rng.normal(0,.1,(blocks//2,n))
    noise=np.concatenate([deviations,-deviations])
    rows=np.zeros((blocks,3+4*n));rows[:,3:3+n]=counts
    rows[:,3+3*n:]=-counts*(intercept+slope*tau)-noise
    result=fit_direct_mass(rows,np.zeros(blocks,int),p,6,14)
    assert_allclose(result['inverse_mass'],slope,atol=1e-10)
    assert_allclose(result['intercept'],intercept,atol=1e-10)
    assert_allclose(result['mass_ratio'],2/slope,atol=1e-10)
