import numpy as np
from numpy.testing import assert_allclose
from defect_reference import Cloud, DefectBasis, Parameters, lanczos


def test_free_impurity_energy_density_and_pole_weight():
    p = Parameters(g=0., U=1.)
    basis = DefectBasis(0, 30)
    result = basis.ground(p)
    kappa = np.arcsinh(p.U/(2*p.t))
    density = np.tanh(kappa)*np.exp(-2*kappa*abs(basis.x))
    assert_allclose(result['E'], -np.sqrt(4*p.t**2+p.U**2), atol=1e-11)
    assert_allclose(result['p'], density, atol=1e-10)
    assert_allclose(result['Zi'], result['p'], atol=1e-13)
    assert_allclose(result['Zb'], 1, atol=1e-13)


def test_atomic_displaced_oscillator_and_local_sum_rules():
    p = Parameters(t=0., g=.7, omega=1.2, U=.8)
    basis = DefectBasis(12, 0)
    result = basis.ground(p)
    assert_allclose(result['E'], -p.U-p.g**2/p.omega, atol=2e-10)
    assert_allclose(result['Zb'], np.exp(-(p.g/p.omega)**2), atol=2e-9)
    assert_allclose(result['nph'], (p.g/p.omega)**2, atol=2e-9)
    energy, weight, _, _ = lanczos(basis.hamiltonian(p), basis.source(), 100)
    assert_allclose([weight.sum(), energy@weight, energy**2@weight], [1., -p.U, p.U**2+p.g**2], atol=1e-12)


def test_hermiticity_site_dependent_moments_and_hellmann_feynman():
    basis = DefectBasis(4, 5)
    p = Parameters()
    h = basis.hamiltonian(p)
    rng = np.random.default_rng(193)
    a, b = rng.normal(size=(2,basis.dim))
    assert_allclose(a@(h@b), b@(h@a), atol=1e-11)
    for site in [0, 1, -2]:
        v = basis.source(site); w = h@v
        potential = -p.U if site == 0 else 0
        assert_allclose(v@w, potential, atol=1e-13)
        assert_allclose(w@w, potential**2+2*p.t**2+p.g**2, atol=1e-13)
    step=1e-4
    low=basis.ground(Parameters(U=1-step)); high=basis.ground(Parameters(U=1+step))
    mid=basis.ground(p)
    assert_allclose(-(high['E']-low['E'])/(2*step), mid['p0'], atol=2e-8)
