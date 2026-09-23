import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import eigh
from defect_reference import DefectBasis, Parameters
from defect_momentum import (source, window_moments, projected_lanczos, columns,
                             complete_columns, fourier_window, analytic_free_matrix)


def test_projected_resolvent_and_coherent_fourier_against_full_diagonalization():
    basis = DefectBasis(2, 4)
    p = Parameters(g=.7, U=1.2)
    h = basis.hamiltonian(p)
    matrix = np.column_stack([h@v for v in np.eye(basis.dim)])
    en, v = eigh(matrix)
    axis = np.linspace(-3.3, 5., 43)
    w = 2
    indices = (np.arange(-w, w+1)+basis.radius)*basis.cloud.dim
    rec = [projected_lanczos(h, basis.source(j), indices, basis.dim) for j in range(w+1)]
    responses = [columns(r, axis, .25) for r in rec]
    green = complete_columns(np.stack([r[0] for r in responses]))
    exact = np.einsum('in,jn,ne->ije', v[indices], v[indices], 1/(axis[None]+.25j-en[:, None]))
    assert_allclose(green, exact, atol=2e-11)
    k = np.linspace(0, np.pi, 13)
    projected = fourier_window(green, k, w)
    for ik, kk in enumerate(k):
        initial = source(basis, kk, w, 'cos')+1j*source(basis, kk, w, 'sin')
        weights = abs(v.T@initial)**2
        expected = weights@(1/(axis[None]+.25j-en[:, None]))
        assert_allclose(projected[ik], expected, atol=2e-11)
        assert_allclose(window_moments([kk], w, p.t, p.g, p.U)[0],
                        [1., en@weights, en**2@weights], atol=1e-12)
    assert (-projected.imag/np.pi >= 0).all()


def test_free_infinite_chain_and_impurity_against_large_box():
    axis = np.linspace(-3.2, 3.2, 65)
    for U in [0., 1.]:
        basis = DefectBasis(0, 120)
        h = basis.hamiltonian(Parameters(g=0, U=U))
        indices = np.arange(-3, 4)+basis.radius
        rec = projected_lanczos(h, basis.source(2), indices, 241)
        green, moments, residual = columns(rec, axis, .3)
        exact = analytic_free_matrix(axis, .3, 3, U=U)[:, 5]
        assert_allclose(green, exact, atol=3e-13)
        assert residual.max() < 1e-12


def test_subwindow_contraction_and_window_edge_moments():
    basis = DefectBasis(2, 7)
    p = Parameters(g=.9, U=.6)
    for w in [0, 1, 3, 5]:
        for k in [0., .37, np.pi]:
            v = source(basis, k, w, 'cos')+1j*source(basis, k, w, 'sin')
            hv = basis.hamiltonian(p)@v
            assert_allclose(window_moments([k], w, p.t, p.g, p.U)[0],
                            [np.vdot(v, v).real, np.vdot(v, hv).real, np.vdot(hv, hv).real], atol=1e-13)
    rng = np.random.default_rng(332)
    matrix = rng.normal(size=(9, 9, 3))
    k = np.array([0., .8, 2.])
    assert_allclose(fourier_window(matrix, k, 2), fourier_window(matrix[2:7, 2:7], k, 2))
