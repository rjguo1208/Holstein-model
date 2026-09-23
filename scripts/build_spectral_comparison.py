"""Publish all 41 refined DiagMC columns alongside independent VED.

This is a presentation step: it does not fit spectra or interpolate momenta.
Run with the packages pinned in research/spectral-fullgrid/requirements.txt.
"""
from hashlib import sha256
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from spectral_plot import vector_map, save_figure
from matplotlib.colors import LogNorm, Normalize
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "site/data"
PLOTS = ROOT / "site/results"
COUPLINGS = [.25, .5]


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def edges(centers):
    middle = (centers[1:] + centers[:-1]) / 2
    return np.r_[centers[0] - (middle[0] - centers[0]), middle,
                 centers[-1] + (centers[-1] - middle[-1])]


def assemble(coupling):
    sources = {
        "baseline": DATA / f"diagmc41-lambda{coupling:.2f}_map.npz",
        "refinement": DATA / f"diagmc-fullgrid-lambda{coupling:.2f}.npz",
        "ved": DATA / f"lambda{coupling:.2f}_map.npz",
    }
    with np.load(sources["baseline"]) as b, np.load(sources["refinement"]) as p, np.load(sources["ved"]) as v:
        indices = p["indices"]
        assert indices.dtype.kind in "iu"
        np.testing.assert_array_equal(indices, np.arange(41))
        assert b["A"].shape == (3, 41, 701)
        assert p["selected"].shape == (3, 41, 701)
        assert p["quantiles"].shape == (4, 3, 41, 701)
        for source in (p, v):
            assert float(source["coupling"]) == coupling
            np.testing.assert_array_equal(b["energy"], source["energy"])
        np.testing.assert_array_equal(b["k"], v["k"])
        np.testing.assert_array_equal(b["k"][indices], p["k"])
        np.testing.assert_array_equal(b["eta"], p["eta"])
        np.testing.assert_array_equal(b["A"][:, indices], p["baseline"])
        reference = np.stack([v["A"][np.flatnonzero(v["eta"] == eta).item()] for eta in b["eta"]])
        np.testing.assert_array_equal(reference[:, indices], p["reference"])
        for key in ("t", "omega0"):
            assert float(b[key]) == float(v[key]) == 1.

        refined = np.ones(len(b["k"]), dtype=bool)
        spectrum = p["selected"].copy()
        np.testing.assert_array_equal(spectrum, p["selected"])
        assert np.isfinite(spectrum).all() and (spectrum >= 0).all()
        assert np.isfinite(reference).all() and (reference >= 0).all()
        bands = {}
        for suffix, probability in [("16", .16), ("84", .84)]:
            iq = np.flatnonzero(p["bootstrap_probabilities"] == probability).item()
            bands[f"bootstrap_{suffix}"] = p["quantiles"][iq].copy()
            assert np.isfinite(bands[f"bootstrap_{suffix}"]).all()
        assert np.all(bands["bootstrap_16"] <= bands["bootstrap_84"])

        bundle = dict(A=spectrum, A_ved=reference, k=b["k"], energy=b["energy"],
                      eta=b["eta"], coupling=coupling, t=1., omega0=1.,
                      refined_mask=refined, refined_indices=indices,
                      source_per_k=np.full(41, "refinement"),
                      bootstrap_replicates=np.full(41, 128),
                      spectral_resolution_validated=False,
                      method="41 measured momenta: all use independent 4x sampling and frozen refined fits",
                      **bands)
        output = DATA / f"diagmc-updated41-lambda{coupling:.2f}-map.npz"
        np.savez_compressed(output, **bundle)
        denominator = np.trapezoid(reference, b["energy"], axis=2)
        error = np.trapezoid(abs(spectrum - reference), b["energy"], axis=2) / denominator
        record = dict(coupling=coupling, shape=list(spectrum.shape),
                      refined_indices=indices.tolist(), unchanged_baseline_columns=0,
                      refined_columns_exact=True, all_columns_refined=True,
                      axes_identical=True, refined_ved_matches_full_reference=True,
                      input_sha256={str(path.relative_to(ROOT)): digest(path) for path in sources.values()},
                      output=str(output.relative_to(ROOT)), output_sha256=digest(output),
                      relative_l1_by_eta_and_k=error.tolist())
    return bundle, record


def plot(bundles, eta, logarithmic=False):
    ie = np.flatnonzero(bundles[0]["eta"] == eta).item()
    upper = max(float(b[key][ie].max()) for b in bundles for key in ["A", "A_ved"])
    norm = LogNorm(vmin=1e-3, vmax=upper) if logarithmic else Normalize(vmin=0., vmax=upper)
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.6), sharex=True, sharey=True,
                             constrained_layout=True)
    meshes=[]
    for row, b in enumerate(bundles):
        for col, (key, title) in enumerate([("A", "Refined DiagMC"), ("A_ved", "VED reference")]):
            ax = axes[row, col]
            mesh = vector_map(ax,b["k"]/np.pi,b["energy"],b[key][ie],"viridis",norm)
            meshes.append(mesh)
            bare_k=np.linspace(0.,np.pi,401)
            ax.plot(bare_k/np.pi,-2*np.cos(bare_k),color="white",linestyle="--",linewidth=1.2,
                    label=r"Bare $\epsilon_k/t=-2\cos k$",
                    path_effects=[path_effects.Stroke(linewidth=2.,foreground="black"),path_effects.Normal()])
            ax.legend(loc="upper left",fontsize=8,framealpha=.85)
            ax.set(xlim=(0., 1.), ylim=(b["energy"][0], b["energy"][-1]))
            ax.set_title(rf"$\lambda={b['coupling']:g}$: {title}", pad=15)
            if col == 0:
                ax.set_ylabel(r"$\omega/t$")
            if row == 1:
                ax.set_xlabel(r"$k/\pi$")
            ax.set_xticks([0., .25, .5, .75, 1.], ["0", "1/4", "1/2", "3/4", "1"])
    bar = fig.colorbar(mesh, ax=axes, label=r"$t A_\eta(k,\omega)$", shrink=.91, pad=.02)
    bar.solids.set_rasterized(False)
    scale = "logarithmic color" if logarithmic else "linear color"
    fig.suptitle(rf"Spectral function comparison: $\eta={eta:g}t$ ({scale})" + "\n"
                 "All 41 momenta: 4x DiagMC sampling and the same refinement protocol", fontsize=12)
    name = f"spectral-comparison-eta{eta:g}-{'log' if logarithmic else 'linear'}"
    save_figure(fig,PLOTS/name,meshes+[bar.solids])
    plt.close(fig)
    return dict(name=name, eta=eta, scale=scale, vmin=float(norm.vmin), vmax=float(norm.vmax),
                one_color_scale_for_all_four_panels=True,
                bare_band="epsilon_k/t=-2 cos(k), analytic 401 points, dashed",
                vector_color_runs=[m.vector_audit for m in meshes],
                files={extension: digest(PLOTS / f"{name}.{extension}") for extension in ["svg", "pdf", "png"]})


def main():
    plt.rcParams.update({"font.size": 11, "svg.fonttype": "none",
                         "svg.hashsalt": "holstein-spectral-comparison-v2"})
    summary_path = DATA / "diagmc-fullgrid-summary.json"
    summary = json.loads(summary_path.read_text())
    assert summary["complete"] and summary["indices"] == list(range(41))
    assert summary["bootstrap_points"] == 82 and summary["bootstrap_draws"] == 10496
    assert summary["bootstrap_failures"] == 0
    bundles, records = zip(*(assemble(coupling) for coupling in COUPLINGS))
    for record in records:
        points = sorted((p for p in summary["points"] if p["coupling"] == record["coupling"]),
                        key=lambda p: p["index"])
        assert [p["index"] for p in points] == list(range(41))
        np.testing.assert_allclose(record["relative_l1_by_eta_and_k"],
                                   np.array([p["errors"]["selected"] for p in points]).T,
                                   rtol=1e-12, atol=1e-12)
    for key in ["k", "energy", "eta", "refined_mask"]:
        np.testing.assert_array_equal(bundles[0][key], bundles[1][key])
    plots = [plot(bundles, .25), plot(bundles, 1.), plot(bundles, .25, logarithmic=True)]
    provenance = dict(complete=True, script_sha256=digest(Path(__file__)),
                      plot_helper_sha256=digest(Path(__file__).with_name("spectral_plot.py")),
                      points=41, refined_points_per_coupling=41, baseline_points_per_coupling=0,
                      scientific_summary_sha256=digest(summary_path),
                      eta=bundles[0]["eta"].tolist(), energy_window=[-3.25, 5.5],
                      no_interpolation=True, no_column_normalization=True, no_energy_shifts=True,
                      display_script_does_not_sample_or_fit=True, ved_used_only_for_comparison=True,
                      spectral_resolution_validated=False, couplings=list(records), plots=plots,
                      bootstrap_note="16-84% empirical ranges: 128 draws at every momentum; representation bias not included")
    (DATA / "spectral-comparison-provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print("Verified all 41 columns refined per coupling; created three paired maps.")


if __name__ == "__main__":
    main()
