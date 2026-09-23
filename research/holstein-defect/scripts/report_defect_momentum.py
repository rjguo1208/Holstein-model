"""Audit and plot dense VED momentum spectra; retain absolute spectral weights."""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from hashlib import sha256
import json
import os
from pathlib import Path
import sys
import numpy as np
from scipy.linalg import eigh_tridiagonal
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import Normalize, TwoSlopeNorm
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from defect_momentum import pole_columns, complete_columns, fourier_window, window_moments
from defect_analysis import public
from spectral_plot import vector_map, save_figure

ENERGY = np.linspace(-4.5, 5.5, 1601)
K = np.linspace(0, np.pi, 201)
ETA = np.array([.1, .15, .25, .5, 1.])


def relative(a, b):
    return np.trapezoid(abs(a-b), ENERGY, axis=-1)/np.trapezoid(b, ENERGY, axis=-1)


def assemble(task):
    data, out, lam, U, nh, radius, window = task
    label = f'lambda{lam:.2f}_U{U:g}_Nh{nh}_R{radius}_W{window}'
    shape = (len(ETA), window+1, 2*window+1, len(ENERGY))
    greens = np.empty(shape, complex)
    half = np.empty(shape, complex)
    moments = []
    residuals = []
    for j in range(window+1):
        with np.load(data/(label+f'_i{j:02d}.npz')) as f:
            rec = dict(f)
        en, residue, moment, tail = pole_columns(rec)
        en2, residue2, _, _ = pole_columns(rec, 400)
        moments.append(moment)
        residuals.append([])
        for ie, eta in enumerate(ETA):
            kernel = 1/(ENERGY[None]+1j*eta-en[:, None])
            greens[ie, j] = residue@kernel
            half[ie, j] = residue2@(1/(ENERGY[None]+1j*eta-en2[:, None]))
            residuals[-1].append(float(abs(tail@kernel).max()))
    full = complete_columns(greens)
    del greens
    half_full = complete_columns(half)
    del half
    full_moments = complete_columns(np.array(moments))
    fields = dict(k=K, energy=ENERGY, eta=ETA, coupling=lam, U=U, generations=nh,
                  electron_radius=radius, window_radius=window, t=1., omega0=1.,
                  method='VED coherent finite-window momentum spectrum', steps=800)
    audit = dict(label=label, coupling=lam, U=U, generations=nh, radius=radius, window=window,
                 max_residual_by_eta=np.max(residuals, axis=0),
                 max_spectral_error_bound_from_residual_by_eta=np.sqrt(2*window+1)*np.max(residuals, axis=0)/ETA/np.pi,
                 max_reciprocity_error_by_eta=np.max(abs(full-full.swapaxes(1, 2)), axis=(1, 2, 3)))
    for w in sorted({16, window}):
        A = -fourier_window(full, K, w).imag/np.pi
        A_half = -fourier_window(half_full, K, w).imag/np.pi
        mm = fourier_window(full_moments, K, w)
        expected = window_moments(K, w, g=np.sqrt(2*lam), U=U)
        error = float(abs(mm-expected).max())
        if error > 2e-9 or A.min() < -1e-8 or not np.isfinite(A).all():
            raise ValueError(f'Nonphysical assembled map {label}: moments={error}, min={A.min()}')
        fields['A_W'+str(w)] = A
        fields['moments_W'+str(w)] = mm
        fields['steps_relative_l1_W'+str(w)] = relative(A_half, A)
        audit['W'+str(w)] = dict(moment_max_abs_error=error, minimum_spectrum=float(A.min()),
                                 steps_max_relative_l1_by_eta=relative(A_half, A).max(axis=1),
                                 displayed_weight_range_by_eta=np.stack([np.trapezoid(A, ENERGY, axis=2).min(axis=1),
                                                                        np.trapezoid(A, ENERGY, axis=2).max(axis=1)], axis=1))
    path = out/(label+'.npz')
    np.savez_compressed(path, **fields)
    audit.update(file=path.name, sha256=sha256(path.read_bytes()).hexdigest())
    return audit


def direct_spectrum(path):
    A = np.zeros((len(ETA), len(ENERGY)))
    with np.load(path) as f:
        for parity in ['cos', 'sin']:
            if parity+'_alpha' not in f:
                continue
            en, v = eigh_tridiagonal(f[parity+'_alpha'], f[parity+'_beta'][:-1])
            weights = v[0]**2*float(f[parity+'_initial_norm'])**2
            for ie, eta in enumerate(ETA):
                A[ie] += np.sum(weights[:, None]*eta/np.pi/((ENERGY[None]-en[:, None])**2+eta**2), axis=0)
    return A


def draw(out, bundles, eta, difference=False):
    ie = np.flatnonzero(ETA == eta).item()
    if difference:
        fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.8), layout='constrained', sharex=True, sharey=True)
        values = [b[1]['A_W16'][ie]-b[0]['A_W16'][ie] for b in bundles]
        vmax = max(np.max(abs(a)) for a in values)
        norm = TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)
        panels = [(axes[i], values[i], rf'$\lambda={lam:g}$: defect $-$ clean') for i, lam in enumerate([.25, .5])]
        cmap = 'RdBu_r'
    else:
        fig, axes = plt.subplots(2, 2, figsize=(10.4, 8.8), layout='constrained', sharex=True, sharey=True)
        vmax = max(b[u]['A_W16'][ie].max() for b in bundles for u in [0, 1])
        norm = Normalize(0, vmax)
        panels = [(axes[row, col], b[U]['A_W16'][ie],
                   rf'$\lambda={lam:g}$: '+('single defect, $U=t$' if U else 'clean, $U=0$'))
                  for col, (lam, b) in enumerate(zip([.25, .5], bundles)) for row, U in enumerate([1, 0])]
        cmap = 'magma'
    meshes = []
    for ax, a, title in panels:
        mesh = vector_map(ax, K/np.pi, ENERGY, a, cmap, norm)
        meshes.append(mesh)
        ax.plot(K/np.pi, -2*np.cos(K), '--', color='white', lw=1.15,
                path_effects=[pe.Stroke(linewidth=2., foreground='black'), pe.Normal()])
        ax.set(title=title, xlim=(0, 1), ylim=(-4.5, 5.5), xlabel=r'$k/\pi$',
               xticks=[0, .25, .5, .75, 1], xticklabels=['0', '1/4', '1/2', '3/4', '1'])
    for ax in np.asarray(axes).reshape((-1, 2))[:, 0]:
        ax.set_ylabel(r'Energy $\omega/t$ (vacuum = 0)')
    label = r'$t\,[A_{L,U=t}(k,\omega)-A_{L,U=0}(k,\omega)]$' if difference else r'$t A_{L,\eta}(k,\omega)$'
    bar = fig.colorbar(mesh, ax=np.asarray(axes).ravel().tolist(), label=label, shrink=.88, pad=.02)
    bar.solids.set_rasterized(False)
    fig.suptitle(rf'VED: coherent window $L=33$, 201 momenta, $\eta={eta:g}t$'+'\n'+
                 r'White dashed: bare host $\epsilon_k/t=-2\cos k$; $t=\omega_0=1$', fontsize=11)
    name = f'defect-momentum-ved-eta{eta:g}'+('-difference' if difference else '')
    save_figure(fig, out/name, meshes+[bar.solids])
    plt.close(fig)
    return dict(name=name, eta=eta, difference=difference, vmax=vmax,
                vector_color_runs=[m.vector_audit for m in meshes])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--jobs', type=int, default=8)
    args = parser.parse_args()
    mf = json.loads((args.data/'manifest.json').read_text())
    assert mf['complete'] and not mf['pilot']
    for r in mf['records']:
        assert sha256((args.data/r['file']).read_bytes()).hexdigest() == r['sha256']
    args.out.mkdir(parents=True, exist_ok=False)
    summary = dict(complete=False, job_id=os.environ.get('SLURM_JOB_ID'), method='VED',
                   observable='A_L(k,w) from all off-diagonal G_ij in a normalized coherent window',
                   window=33, comparison_window=65, k_points=201, energy_points=1601,
                   k_path=[0., float(np.pi)], eta=ETA, energy_range=[-4.5, 5.5],
                   source_manifest_sha256=sha256((args.data/'manifest.json').read_bytes()).hexdigest(),
                   report_source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(), records=[], plots=[])
    tasks = [(args.data, args.out, lam, U, nh, R, W) for nh, R, W in mf['settings']
             for lam in [.25, .5] for U in [0., 1.]]
    def save():
        (args.out/'summary.json').write_text(json.dumps(public(summary), indent=2)+'\n')
    save()
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        for f in as_completed([pool.submit(assemble, t) for t in tasks]):
            item = f.result()
            summary['records'].append(item)
            save()
            print('ASSEMBLED', item['label'], flush=True)
    all_data = {}
    for r in summary['records']:
        with np.load(args.out/r['file']) as f:
            all_data[(r['coupling'], r['U'], r['generations'], r['radius'])] = dict(f)
    summary['convergence'] = []
    for lam in [.25, .5]:
        for U in [0., 1.]:
            base = all_data[lam, U, 12, 96]
            coarse = all_data[lam, U, 10, 64]['A_W16']
            middle = all_data[lam, U, 12, 64]['A_W16']
            final = base['A_W16']
            checks = dict(coupling=lam, U=U,
                          cloud_max_relative_l1_by_eta=relative(coarse, middle).max(axis=1),
                          boundary_max_relative_l1_by_eta=relative(middle, final).max(axis=1),
                          steps_max_relative_l1_by_eta=base['steps_relative_l1_W16'].max(axis=1),
                          window_change_max_relative_l1_by_eta=relative(base['A_W32'], final).max(axis=1), direct=[])
            for r in mf['records']:
                if r['kind'] != 'direct' or r['coupling'] != lam or r['U'] != U:
                    continue
                k_index = int(round(r['k_over_pi']*200))
                assert abs(K[k_index]/np.pi-r['k_over_pi']) < 1e-14
                direct = direct_spectrum(args.data/r['file'])
                error = relative(direct, final[:, k_index])
                if error.max() > 1e-5:
                    raise ValueError(f'Independent momentum injection disagrees: {r["file"]} {error}')
                checks['direct'].append(dict(k_over_pi=r['k_over_pi'], relative_l1_by_eta=error))
            assert len(checks['direct']) == 5
            summary['convergence'].append(checks)
    bundles = [{U: all_data[lam, float(U), 12, 96] for U in [0, 1]} for lam in [.25, .5]]
    plt.rcParams.update({'font.size':10, 'mathtext.fontset':'stix', 'svg.fonttype':'path',
                         'svg.hashsalt':'holstein-defect-momentum-v1', 'savefig.facecolor':'white'})
    for eta in [.1, .15, .25]:
        summary['plots'].append(draw(args.out, bundles, eta))
    summary['plots'].append(draw(args.out, bundles, .15, difference=True))
    summary.update(complete=True, no_k_interpolation=True, no_column_normalization=True,
                   no_energy_shift=True, contains_diagmc=False,
                   notes='Window changes define a different observable, not a numerical convergence error. VED cutoff differences are not rigorous bounds on the infinite phonon-space error.')
    save()
    print('COMPLETE VED report', args.out, flush=True)


if __name__ == '__main__':
    main()
