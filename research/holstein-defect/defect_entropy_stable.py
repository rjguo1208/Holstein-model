"""Numerical fallback for the SAME frozen moment-constrained MaxEnt objective.

Only the Newton linear solve and objective-decrease evaluation differ from
refinement.entropy_scan. Grids, whitening, projected SVD, moments and alpha
are identical. Extended-precision elimination avoids clipping small Hessian
eigenvalues; a centered cumulant avoids subtracting nearly equal objectives.
"""
import numpy as np
from scipy.special import logsumexp
from defect_entropy import local_moments
from continuation import bin_kernel
from spectral_cv import whitening


def long_solve(matrix, rhs):
    """Partial-pivot Gaussian elimination retaining long-double arithmetic."""
    a = np.array(matrix, dtype=np.longdouble, copy=True)
    b = np.array(rhs, dtype=np.longdouble, copy=True)
    for k in range(len(b)):
        pivot = k + int(np.argmax(abs(a[k:, k])))
        if not np.isfinite(a[pivot, k]) or abs(a[pivot, k]) < 1e-28:
            raise np.linalg.LinAlgError('Singular extended-precision Newton system')
        if pivot != k:
            a[[k, pivot]] = a[[pivot, k]]
            b[[k, pivot]] = b[[pivot, k]]
        factors = a[k + 1:, k] / a[k, k]
        a[k + 1:, k + 1:] -= factors[:, None] * a[k, k + 1:]
        b[k + 1:] -= factors * b[k]
    x = np.empty_like(b)
    for k in range(len(b) - 1, -1, -1):
        x[k] = (b[k] - a[k, k + 1:] @ x[k + 1:]) / a[k, k]
    return x


def solve(tau, green, cov, p, alpha, size=281, upper=10., tau_limit=8.,
          tolerance=1e-8, target_moments=None, energies=None, max_iterations=800):
    lower = -np.hypot(2*p['t'],p['U']) - p['g'] ** 2 / p['omega']
    en = np.linspace(lower, upper, size) if energies is None else np.asarray(energies)
    use = tau <= tau_limit
    W = whitening(green[use], cov[np.ix_(use, use)], tolerance)
    K = bin_kernel(tau, p['tau_max'] / p['bins'], en, p['mu'])
    B = W @ K[use]
    y = W @ green[use]
    M = np.vstack([np.ones(len(en)), en, en ** 2])
    moments = local_moments(p) if target_moments is None else np.asarray(target_moments)
    assert np.isclose(moments[0], 1.) and alpha > 0
    Q, R = np.linalg.qr(M.T, mode='reduced')
    c = np.linalg.solve(R.T, moments)
    projected = B - (B @ Q) @ Q.T
    yp = y - (B @ Q) @ c
    U, s, V = np.linalg.svd(projected, full_matrices=False)
    keep = s > s[0] * 1e-12
    s, U, V = s[keep], U[:, keep], V[keep]
    F = np.asarray(np.vstack([V, Q[:, 1:].T]), dtype=np.longdouble)
    d = np.asarray(np.r_[(U.T @ yp) / s, c[1:]], dtype=np.longdouble)
    reg = np.asarray(np.r_[alpha / s ** 2, 0., 0.], dtype=np.longdouble)
    logm = np.full(len(en), -np.log(len(en)), dtype=np.longdouble)
    theta = np.zeros(len(d), dtype=np.longdouble)

    def evaluate(x):
        z = logm - F.T @ x
        logw = z - logsumexp(z)
        w = np.exp(logw)
        normalization = w.sum()
        logw -= np.log(normalization)
        w /= normalization
        fw = F @ w
        grad = d - fw + reg * x
        centered = F - fw[:, None]
        hessian = (centered * w) @ centered.T + np.diag(reg)
        return grad, w, hessian, centered, logw

    history = []
    for iteration in range(max_iterations):
        grad, w, hessian, centered, logw = evaluate(theta)
        mr = M @ np.asarray(w, dtype=float) - moments
        gmax = float(np.max(abs(grad)))
        history.append(dict(iteration=iteration, gradient=gmax, moment_error=float(max(abs(mr)))))
        if gmax < 2e-11 and max(abs(mr)) < 1e-6:
            break
        scale = np.sqrt(np.maximum(np.diag(hessian), np.longdouble('1e-30')))
        normalized=hessian/np.outer(scale,scale)
        for damping in [0.,1e-16,1e-14,1e-12,1e-10,1e-8,1e-6]:
            try:
                step=long_solve(normalized+damping*np.eye(len(grad)),-grad/scale)/scale
            except np.linalg.LinAlgError:
                continue
            slope=grad@step
            if np.isfinite(step).all() and slope<0:break
        else:raise ValueError('No finite descent direction for identical convex objective')
        logchange=centered.T@step
        slope=grad@step
        fraction = np.longdouble(1)
        for _ in range(80):
            delta = fraction * step
            centered_change = -(centered.T @ delta)
            if max(abs(centered_change)) < .1:
                # E[centered_change] = 0; evaluate only the quadratic remainder.
                remainder = np.sum(w * (np.expm1(centered_change) - centered_change))
                cumulant = np.log1p(remainder)
            else:
                # Retain log weights even when exp(logw) underflows: a large
                # trial step can return that support to appreciable weight.
                cumulant = logsumexp(logw + centered_change)
            change = grad @ delta + cumulant + .5 * np.dot(reg * delta, delta)
            if change <= np.longdouble('1e-4') * fraction * slope:
                candidate = theta + delta
                if np.array_equal(candidate, theta):
                    raise ValueError('Extended-precision Newton step underflowed before convergence')
                theta = candidate
                break
            fraction *= .5
        else:
            raise ValueError('Centered objective line search failed')
    grad, w, _, _, _ = evaluate(theta)
    weights = np.asarray(w, dtype=float)
    mr = M @ weights - moments
    predicted = K @ weights
    residual = W @ (predicted[use] - green[use])
    converged = bool(max(abs(grad)) < 1e-8 and max(abs(mr)) < 1e-6)
    positive = w > 0
    return dict(alpha=float(alpha), energies=en, weights=weights, predicted=predicted,
                rank=len(W), chi2_per_mode=float(residual @ residual / len(W)),
                moment_residuals=mr, lower_bound=lower, upper_support=float(en[-1]),
                tau_limit=tau_limit, edge_weight=float(weights[-3:].sum()),
                method='maximum entropy; exact moments', converged=converged,
                dual_gradient=float(max(abs(grad))), iterations=iteration + 1,
                numerical_method='same frozen dual; long-double pivoted Newton solve and centered objective differences',
                primal_objective=float(.5 * (residual @ residual) + alpha * np.sum(w[positive] * (np.log(w[positive]) - logm[positive]))),
                iteration_tail=history[-10:])

