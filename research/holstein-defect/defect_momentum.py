"""Coherent finite-window momentum spectra of a single Holstein defect.

The infinite-chain observable is <v_k|(z-H)^-1|v_k>, where v_k has
zero phonons and amplitudes exp(ikx)/sqrt(L) on -W,...,W (L=2W+1).
The electron boundary R and phonon cloud Nh are independent VED cutoffs.
No periodic repetition of the impurity and no Fourier transform of LDOS.
"""
from __future__ import annotations
import numpy as np
from scipy.linalg import eigh_tridiagonal


def window_moments(k, radius, t=1., g=0., U=1.):
    """Exact unbroadened M0,M1,M2, including hopping outside the window."""
    k = np.atleast_1d(k)
    x = np.arange(-radius-1, radius+2)
    v = np.zeros((len(k), len(x)), complex)
    v[:, 1:-1] = np.exp(1j*k[:, None]*x[None, 1:-1])/np.sqrt(2*radius+1)
    hv = -U*(x == 0)*v
    hv[:, 1:] -= t*v[:, :-1]
    hv[:, :-1] -= t*v[:, 1:]
    return np.stack([np.ones(len(k)), np.real(np.sum(v.conj()*hv, axis=1)),
                     np.sum(abs(hv)**2, axis=1)+g*g], axis=1)


def source(basis, k, radius, parity):
    if radius >= basis.radius or parity not in ['cos', 'sin']:
        raise ValueError('Need a window strictly inside the electron boundary')
    v = np.zeros((len(basis.x), basis.cloud.dim))
    x = np.arange(-radius, radius+1)
    wave = np.cos(k*x) if parity == 'cos' else np.sin(k*x)
    v[x+basis.radius, 0] = wave/np.sqrt(len(x))
    return v.ravel()


def projected_lanczos(h, initial, indices, steps):
    """Keep only selected entries of Krylov vectors, not the large vectors.

    The last beta is retained to estimate the resolvent residual. Short
    recurrence loss of global orthogonality is controlled by comparison to
    doubled recursion and independently injected cosine/sine states.
    """
    q = np.array(initial, float)
    norm = np.linalg.norm(q)
    if norm <= 1e-14:
        raise ValueError('Zero initial vector')
    q /= norm
    previous = np.zeros_like(q)
    beta_previous = 0.
    alpha, beta, projection = [], [], []
    for _ in range(min(steps, len(q))):
        projection.append(q[indices].copy())
        v = h@q-beta_previous*previous
        a = float(q@v)
        v -= a*q
        # A second local subtraction reduces roundoff without storing Q.
        correction = float(q@v)
        a += correction
        v -= correction*q
        b = float(np.linalg.norm(v))
        alpha.append(a)
        beta.append(b)
        if b < 1e-13:
            break
        previous, q, beta_previous = q, v/b, b
    return dict(alpha=np.asarray(alpha), beta=np.asarray(beta),
                projection=np.asarray(projection), initial_norm=float(norm))


def pole_columns(recursion, steps=None):
    """Reuse the tridiagonal eigensystem for all energies and broadenings."""
    n = min(steps or len(recursion['alpha']), len(recursion['alpha']))
    alpha = recursion['alpha'][:n]
    beta = recursion['beta'][:n]
    en, vec = eigh_tridiagonal(alpha, beta[:-1])
    residue = (recursion['projection'][:n].T@vec)*vec[0]
    residue *= recursion['initial_norm']
    moments = residue@np.stack([np.ones(n), en, en**2], axis=1)
    tail = beta[-1]*vec[-1]*vec[0]*recursion['initial_norm']
    return en, residue, moments, tail


def columns(recursion, energy, eta, steps=None):
    """Projected resolvent column, its first three moments and residual norm."""
    en, residue, moments, tail = pole_columns(recursion, steps)
    kernel = 1./(np.asarray(energy)[None, :]+1j*eta-en[:, None])
    green = residue@kernel
    residual = abs(tail@kernel)
    return green, moments, residual


def complete_columns(positive_columns):
    """Reflect source sites j>=0 into j<0, retaining target-site ordering."""
    # shape (..., source>=0, target=-W..W, energy/moment)
    a = np.asarray(positive_columns)
    w = a.shape[-2]//2
    if a.shape[-3] != w+1:
        raise ValueError('Incompatible source and target windows')
    full = np.concatenate([a[..., 1:, ::-1, :][..., ::-1, :, :], a], axis=-3)
    return full


def fourier_window(matrix, k, radius):
    """Contract the full G_ij or A_ij; no momentum interpolation is used."""
    a = np.asarray(matrix)
    length = 2*radius+1
    w = a.shape[-2]//2
    if radius > w or a.shape[-3] != 2*w+1:
        raise ValueError('Requested window outside projected matrix')
    a = a[..., w-radius:w+radius+1, w-radius:w+radius+1, :]
    # For a real, time-reversal invariant Hamiltonian G_ij=G_ji. Only the
    # symmetric part contributes; retain its antisymmetric error separately.
    diagonal = [np.trace(a, axis1=-3, axis2=-2)]
    for d in range(1, length):
        diagonal.append(np.diagonal(a, offset=d, axis1=-3, axis2=-2).sum(axis=-1)
                        +np.diagonal(a, offset=-d, axis1=-3, axis2=-2).sum(axis=-1))
    correlations = np.stack(diagonal, axis=-2)/length
    return np.einsum('kd,...de->...ke', np.cos(np.outer(k, np.arange(length))), correlations)


def analytic_free_matrix(energy, eta, radius, t=1., U=1.):
    """Exact infinite-chain g=0 resolvent including one attractive defect."""
    z = np.asarray(energy)+1j*eta
    root = np.sqrt(z*z-4*t*t)
    root = np.where(root.imag < 0, -root, root)
    q = (-z+root)/(2*t)
    x = np.arange(-radius, radius+1)
    free = q[None, None, :]**abs(x[:, None]-x[None, :])[:, :, None]/root
    to_zero = q[None, :]**abs(x[:, None])/root
    return free-U/(1+U/root)*to_zero[:, None, :]*to_zero[None, :, :]
