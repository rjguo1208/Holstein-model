"""Map conventions and the full-spectrum inverse transform's solvable limit."""
import numpy as np
from continuation import bin_kernel
from map_continuation import lorentz_map, rebin_conditional, scan


def test_lorentz_map_absolute_weight_and_axes():
    # Unequal weights detect accidental per-column normalization.
    axis = np.linspace(-500, 500, 100001)
    energies = np.array([-1., 2.])
    weights = np.array([[.2, .3], [.7, .3]])
    A = lorentz_map(axis, energies, weights, .25)
    assert A.shape == (2, len(axis))
    np.testing.assert_allclose(np.trapezoid(A, axis, axis=1), [.5, 1.], atol=4e-4)
    assert axis[np.argmax(A[1])] == -1.


def test_full_grid_free_particle_with_no_pole_prior():
    p = dict(t=1., omega=1., g=0., k=np.pi/2, mu=-2.8, tau_max=6., bins=24)
    tau = (np.arange(p['bins'])+.5)*.25
    G = bin_kernel(tau, .25, [0.], p['mu'])[:, 0]
    cov = np.diag((.002*G)**2)
    energy = np.linspace(-2, 6, 81)
    fit = scan(tau, G, cov, p, [0.], energies=energy)[0]
    assert abs(fit['energies']@fit['weights']) < 1e-7
    assert abs(fit['weights'][20]-1.) < 1e-5
    assert max(abs(fit['moment_residuals'])) < 1e-7
    np.testing.assert_allclose(fit['predicted'], G, rtol=1e-5)


def test_time_rebin_preserves_integrated_green():
    from analysis import estimate_green
    rng = np.random.default_rng(309)
    rows = rng.uniform(10, 100, (40, 3+16))
    p = dict(bins=16, tau_max=4., t=1., k=.7, mu=-2.8)
    tau, rebinned, p2 = rebin_conditional(rows, p, factor=2)
    old = estimate_green(rows, p)[0]
    new = estimate_green(rebinned, p2)[0]
    np.testing.assert_allclose(new, old.reshape(-1, 2).mean(axis=1))
    assert len(tau) == 8
