"""Compare frozen nonlocal DiagMC reconstructions with independent VED."""
from __future__ import annotations
import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import Normalize
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from defect_analysis import public
from defect_momentum import pole_columns, complete_columns, fourier_window
from continuation import bin_kernel
from spectral_plot import vector_map, save_figure


def green_reference(lam, tau, dt, mu, ved):
    columns = []
    for j in range(33):
        with np.load(ved/f'lambda{lam:.2f}_U1_Nh12_R96_W32_i{j:02d}.npz') as f:
            en, residue, _, _ = pole_columns(dict(f))
        columns.append(residue@bin_kernel(tau, dt, en, mu).T)
    return fourier_window(complete_columns(np.array(columns)), np.linspace(0, np.pi, 201), 16)


def draw(out, bundles, eta):
    ie = np.flatnonzero(bundles[0]['eta'] == eta).item()
    vmax = max(b[key][ie].max() for b in bundles for key in ['A', 'A_ved'])
    norm = Normalize(0, vmax)
    fig, axes = plt.subplots(2, 2, figsize=(10.4, 8.8), layout='constrained', sharex=True, sharey=True)
    meshes = []
    for row, b in enumerate(bundles):
        for col, (key, title) in enumerate([('A', 'New nonlocal DiagMC'), ('A_ved', 'VED reference')]):
            ax = axes[row, col]
            mesh = vector_map(ax, b['k']/np.pi, b['energy'], b[key][ie], 'magma', norm)
            meshes.append(mesh)
            ax.plot(b['k']/np.pi, -2*np.cos(b['k']), '--', color='white', lw=1.15,
                    path_effects=[pe.Stroke(linewidth=2., foreground='black'), pe.Normal()])
            ax.set(title=rf'$\lambda={float(b["coupling"]):g}$: '+title,
                   xlim=(0, 1), ylim=(-4.5, 5.5), xticks=[0, .25, .5, .75, 1],
                   xticklabels=['0', '1/4', '1/2', '3/4', '1'])
            if row == 1:
                ax.set_xlabel(r'$k/\pi$')
            if col == 0:
                ax.set_ylabel(r'Energy $\omega/t$ (vacuum = 0)')
    bar = fig.colorbar(mesh, ax=axes, label=r'$t A_{L,\eta}(k,\omega)$', shrink=.88, pad=.02)
    bar.solids.set_rasterized(False)
    fig.suptitle(rf'Single defect: $U=t$, coherent window $L=33$, $\eta={eta:g}t$'+'\n'+
                 r'201 directly measured momenta; white dashed: bare host $\epsilon_k/t=-2\cos k$', fontsize=11)
    name = f'defect-momentum-diagmc-eta{eta:g}'
    save_figure(fig, out/name, meshes+[bar.solids])
    plt.close(fig)
    return dict(name=name, eta=eta, vmax=vmax, vector_color_runs=[m.vector_audit for m in meshes])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fit', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--reference', type=Path, default=ROOT/'results/defect_momentum_ved_report_02')
    parser.add_argument('--ved', type=Path, default=ROOT/'results/defect_momentum_ved_01')
    args = parser.parse_args()
    mf = json.loads((args.fit/'manifest.json').read_text())
    assert mf['complete'] and len(mf['points']) == 402 and not mf['pilot']
    args.out.mkdir(parents=True, exist_ok=False)
    summary = dict(complete=False, method='new bare DiagMC with open endpoints',
                   job_id=os.environ.get('SLURM_JOB_ID'), k_points=201, window=33,
                   energy_points=1601, production_steps=2048000000, chains=16,
                   bootstrap_per_k=mf['bootstrap'], bootstrap_total=402*mf['bootstrap'],
                   spectral_resolution_validated=False, reference_used_for_fitting=False,
                   no_k_interpolation=True, no_energy_shifts=True, no_column_normalization=True,
                   source_manifest_sha256=sha256((args.fit/'manifest.json').read_bytes()).hexdigest(),
                   report_source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                   couplings={}, plots=[])
    bundles = []
    for lam in [.25, .5]:
        spectra, records, quantiles, greens, errors = [], [], [], [], []
        for index in range(201):
            path = args.fit/f'lambda{lam:.2f}_k{index:03d}'
            selected = next(r for r in mf['points'] if r['coupling'] == lam and r['k_index'] == index)
            assert sha256((path/'selection.json').read_bytes()).hexdigest() == selected['selection_sha256']
            assert sha256((path/'spectra.npz').read_bytes()).hexdigest() == selected['spectra_sha256']
            r = json.loads((path/'summary.json').read_text())
            assert r['complete'] and r['bootstrap_draws'] == mf['bootstrap'] and not r['reference_used']
            with np.load(path/'spectra.npz') as f:
                spectra.append(f['selected'])
                quantiles.append(f['bootstrap_quantiles'])
                greens.append(f['green'])
                errors.append(np.max(f['blocking_errors'], axis=0))
                axis, eta, tau = f['energy'], f['eta'], f['tau']
                sensitivity = np.trapezoid(abs(f['candidates']-f['selected']), axis, axis=-1)/np.trapezoid(f['selected'], axis, axis=-1)
            records.append(dict(k_index=index, k_over_pi=index/200, config=r['config'],
                                alpha=r['fit']['alpha'], chi2_per_mode=r['fit']['chi2_per_mode'],
                                moment_residuals=r['fit']['moment_residuals'],
                                candidate_relative_l1=sensitivity,
                                negative_green_bins=r['negative_green_bins'],
                                green_bins_below_three_sigma=r['green_bins_below_three_sigma']))
        A = np.stack(spectra, axis=1)
        bs = np.stack(quantiles, axis=2)
        with np.load(args.reference/f'lambda{lam:.2f}_U1_Nh12_R96_W32.npz') as v:
            V = v['A_W16']
            k = v['k']
            np.testing.assert_array_equal(axis, v['energy'])
            np.testing.assert_array_equal(eta, v['eta'])
        assert np.isfinite(A).all() and A.min() >= 0
        assert bs.shape == (3, 5, 201, 1601)
        relative = np.trapezoid(abs(A-V), axis, axis=2)/np.trapezoid(V, axis, axis=2)
        p = r['parameters']
        Vg = green_reference(lam, tau, 2*p['tau_max']/p['bins'], p['mu'], args.ved)
        G, se = np.array(greens), np.array(errors)
        z = (G-Vg)/se
        output = dict(A=A, A_ved=V, bootstrap_16=bs[0], bootstrap_84=bs[2],
                      energy=axis, eta=eta, k=k, coupling=lam, U=1., t=1., omega0=1.,
                      window=33, method='new nonlocal bare DiagMC',
                      spectral_resolution_validated=False, bootstrap_draws=mf['bootstrap'],
                      tau=tau, green=G, green_ved=Vg, green_standard_error=se,
                      relative_l1_vs_ved=relative)
        np.savez_compressed(args.out/f'defect-momentum-diagmc-lambda{lam:.2f}.npz', **output)
        summary['couplings'][str(lam)] = dict(records=records, relative_l1_vs_ved=relative,
                   minimum_relative_l1_by_eta=relative.min(axis=1),
                   median_relative_l1_by_eta=np.median(relative, axis=1),
                   maximum_relative_l1_by_eta=relative.max(axis=1),
                   maximum_moment_error=max(max(abs(np.array(rr['moment_residuals']))) for rr in records),
                   green_max_abs_standardized_error_tau_le_6=float(abs(z[:, tau <= 6]).max()),
                   green_median_abs_standardized_error_tau_le_6=float(np.median(abs(z[:, tau <= 6]))),
                   green_negative_bins=int((G < 0).sum()))
        bundles.append(output)
    plt.rcParams.update({'font.size':10, 'mathtext.fontset':'stix', 'svg.fonttype':'path',
                         'svg.hashsalt':'holstein-defect-momentum-mc-v1', 'savefig.facecolor':'white'})
    for eta in [.15, .25, 1.]:
        summary['plots'].append(draw(args.out, bundles, eta))
    summary['complete'] = True
    (args.out/'summary.json').write_text(json.dumps(public(summary), indent=2)+'\n')
    print('COMPLETE nonlocal DiagMC comparison', args.out, flush=True)


if __name__ == '__main__':
    main()
