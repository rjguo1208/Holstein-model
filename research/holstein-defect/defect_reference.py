"""Real-space variational reference for one on-site defect.

Independent cutoffs: the electron lies in [-R,R]; phonon configurations are
generated relative to it with Nh operations. Neither translation reduction
nor a phonon occupation truncation in the Monte Carlo sampler is used here.
The -g convention is unitarily equivalent to +g by b_i -> -b_i.
"""
from __future__ import annotations
from collections import Counter
from dataclasses import dataclass
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import LinearOperator, eigsh


@dataclass(frozen=True)
class Parameters:
    t: float = 1.
    omega: float = 1.
    g: float = np.sqrt(.5)
    U: float = 1.

    def __post_init__(self):
        if not all(np.isfinite([self.t, self.omega, self.g, self.U])) or min(self.t, self.g, self.U) < 0 or self.omega <= 0:
            raise ValueError('Require finite t,g,U >= 0 and omega > 0')


class Cloud:
    def __init__(self, generations):
        if generations < 0: raise ValueError('Negative generation cutoff')
        states, index, frontier = [()], {(): 0}, [()]
        for _ in range(generations):
            new = []
            for state in frontier:
                for target in (tuple(sorted(state+(0,))), tuple(x-1 for x in state), tuple(x+1 for x in state)):
                    if target not in index:
                        index[target] = len(states); states.append(target); new.append(target)
            frontier = new
        self.states, self.dim = states, len(states)
        self.number = np.array([len(s) for s in states], float)
        row, col, val, sr, sc = [], [], [], [], []
        for i, state in enumerate(states):
            target = index.get(tuple(sorted(state+(0,))))
            if target is not None:
                row += [i, target]; col += [target, i]
                val += [-np.sqrt(state.count(0)+1)]*2
            target = index.get(tuple(x-1 for x in state))
            if target is not None: sr.append(target); sc.append(i)
        self.coupling = sparse.csr_matrix((val, (row, col)), shape=(self.dim, self.dim))
        self.shift = sparse.csr_matrix((np.ones(len(sr)), (sr, sc)), shape=(self.dim, self.dim))

    def bulk(self, p, k=0.):
        phase = np.exp(1j*k) if k else 1.
        return (sparse.diags(p.omega*self.number)+p.g*self.coupling
                -p.t*(phase*self.shift+np.conj(phase)*self.shift.T)).tocsr()


class DefectBasis:
    def __init__(self, generations=10, radius=16, cloud=None):
        if radius < 0: raise ValueError('Negative electron radius')
        self.cloud = cloud or Cloud(generations)
        self.radius = radius
        self.x = np.arange(-radius, radius+1)
        self.dim = len(self.x)*self.cloud.dim

    def hamiltonian(self, p):
        cloud = self.cloud
        internal = sparse.diags(p.omega*cloud.number)+p.g*cloud.coupling
        potential = -p.U*(self.x == 0)
        def mv(vector):
            v = np.asarray(vector).reshape(len(self.x), cloud.dim)
            out = (internal@v.T).T+potential[:, None]*v
            out[1:] -= p.t*(cloud.shift@v[:-1].T).T
            out[:-1] -= p.t*(cloud.shift.T@v[1:].T).T
            return out.ravel()
        return LinearOperator((self.dim, self.dim), matvec=mv, rmatvec=mv, dtype=np.float64)

    def source(self, site=0):
        if abs(site) > self.radius: raise ValueError('Source outside electron window')
        v = np.zeros(self.dim); v[(site+self.radius)*self.cloud.dim] = 1
        return v

    def ground(self, p, tolerance=1e-10):
        h = self.hamiltonian(p)
        if self.dim == 1:
            vector = np.ones(1); energy = float((h@vector)[0])
        else:
            guess = np.zeros((len(self.x), self.cloud.dim))
            guess[:, 0] = np.exp(-np.maximum(.15, np.arcsinh(p.U/(2*p.t))) * abs(self.x)) if p.t else (self.x == 0)
            values, vectors = eigsh(h, k=1, which='SA', v0=guess.ravel(), tol=tolerance, ncv=min(28,self.dim), maxiter=5000)
            energy, vector = float(values[0]), vectors[:, 0]
        residual = float(np.linalg.norm(h@vector-energy*vector))
        probability = vector.reshape(len(self.x), self.cloud.dim)**2
        density = probability.sum(axis=1)
        zero_ph = probability[:, 0]
        return dict(E=energy, residual=residual, dimension=self.dim, x=self.x,
            p=density, Zi=zero_ph, Zb=float(zero_ph.sum()), Z0=float(zero_ph[self.radius]),
            p0=float(density[self.radius]), r2=float(self.x**2@density),
            nph=float(probability.sum(axis=0)@self.cloud.number),
            edge_probability=float(density[0]+density[-1]) if self.radius else 0., vector=vector)


def lanczos(h, initial, steps=200):
    """Short recursion for local spectra; steps/convergence are checked externally."""
    q = np.array(initial, float); q /= np.linalg.norm(q)
    previous = np.zeros_like(q); beta_previous = 0.
    alpha, beta = [], []
    for _ in range(min(steps, len(q))):
        v = h@q-beta_previous*previous
        a = float(q@v); v -= a*q
        b = float(np.linalg.norm(v)); alpha.append(a)
        if b < 1e-13: break
        beta.append(b); previous, q, beta_previous = q, v/b, b
    from scipy.linalg import eigh_tridiagonal
    energy, vectors = eigh_tridiagonal(np.asarray(alpha), np.asarray(beta[:len(alpha)-1]))
    return energy, vectors[0]**2, np.asarray(alpha), np.asarray(beta[:len(alpha)-1])


def projected_at_time(h, initial, tau, mu=0.):
    # Exact in the supplied finite variational space; useful for small tests.
    from scipy.linalg import eigh
    matrix = np.column_stack([h@v for v in np.eye(h.shape[0])])
    energy, vectors = eigh(matrix)
    overlap = vectors.T@initial
    return np.array([np.sum(overlap**2*np.exp(-(energy-mu)*t)) for t in np.atleast_1d(tau)])
