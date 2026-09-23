import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import eigh_tridiagonal
from defect_momentum_mc import reference_integral, project, normalize, variance_whitening
from defect_momentum import window_moments
from defect_entropy import entropy_scan
from continuation import bin_kernel


def test_zero_vertex_window_normalization_and_signed_projection():
    p = dict(window=2, bins=8, tau_max=4., mu=-3., U=1.)
    edges = np.linspace(0, 4, 9)
    integral = sum(n*np.diff(-np.exp(-r*edges)/r) for n, r in [(4, 3.), (1, 2.)])
    Z = integral.sum()
    rows = np.zeros((1, 2+5*8))
    rows[0, :2] = [1000, 1000]
    rows[0, 2:10] = 1000*integral/Z
    for k in [0., .4, np.pi]:
        value = normalize(project(rows, p, k)[0], p)
        assert_allclose(value, integral/5/.5, rtol=1e-14)
    rows[0, 10:18] = 10
    assert_allclose(project(rows, p, np.pi)[0, 2:], rows[0, 2:10]-10)


def test_variance_whitening_preserves_negative_data_and_covariance():
    rng = np.random.default_rng(977)
    a = rng.normal(size=(11, 11))
    c = a@a.T+np.eye(11)
    scale = np.logspace(-6, 1, 11)
    c *= np.outer(scale, scale)
    W = variance_whitening(c)
    assert_allclose(W@c@W.T, np.eye(11), atol=2e-14)


def test_window_moment_entropy_uses_physical_moments_without_reference_poles():
    p = dict(window=2, t=1., omega=1., g=0., U=1., bins=24, tau_max=4., mu=-2.8)
    x = np.arange(-20, 21)
    en, v = eigh_tridiagonal(-(x == 0).astype(float), -np.ones(40))
    initial = np.zeros(len(x), complex)
    initial[abs(x) <= 2] = np.exp(.7j*x[abs(x) <= 2])/np.sqrt(5)
    weights = abs(v.T@initial)**2
    tau = (np.arange(24)+.5)/6
    G = bin_kernel(tau, 1/6, en, p['mu'])@weights
    C = np.diag((.002*G+.0001)**2)
    moments = window_moments([.7], 2, g=0., U=1.)[0]
    fits = entropy_scan(tau, G, C, p, [1., 100.], size=161, upper=6.,
                        tau_limit=4., target_moments=moments, whitening_mode='variance')
    for f in fits:
        assert f['converged'] and np.min(f['weights']) >= 0
        assert_allclose(np.array([np.ones(161), f['energies'], f['energies']**2])@f['weights'], moments, atol=1e-7)
