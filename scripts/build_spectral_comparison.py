"""Publish the available refined columns alongside VED on the original grid.

This is a presentation step: it does not fit spectra or interpolate momenta.
Run with the packages pinned in research/spectral-refinement/requirements.txt.
"""
from hashlib import sha256
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
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
        "refinement": DATA / f"diagmc-refinement-lambda{coupling:.2f}.npz",
        "ved": DATA / f"lambda{coupling:.2f}_map.npz",
    }
    with np.load(sources["baseline"]) as b, np.load(sources["refinement"]) as p, np.load(sources["ved"]) as v:
        indices = p["indices"]
        assert indices.dtype.kind in "iu"
        assert len(indices) == len(np.unique(indices)) == 10
        assert b["A"].shape == (3, 41, 701)
        assert p["selected"].shape == (3, 10, 701)
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

        refined = np.zeros(len(b["k"]), dtype=bool)
        refined[indices] = True
        spectrum = b["A"].copy()
        spectrum[:, indices] = p["selected"]
        np.testing.assert_array_equal(spectrum[:, ~refined], b["A"][:, ~refined])
        np.testing.assert_array_equal(spectrum[:, refined], p["selected"])
        assert np.isfinite(spectrum).all() and (spectrum >= 0).all()
        assert np.isfinite(reference).all() and (reference >= 0).all()
        bands = {}
        for suffix, probability in [("16", .16), ("84", .84)]:
            quantile = b[f"bootstrap_{suffix}"].copy()
            iq = np.flatnonzero(p["bootstrap_probabilities"] == probability).item()
            quantile[:, indices] = p["quantiles"][iq]
            bands[f"bootstrap_{suffix}"] = quantile
        assert np.all(bands["bootstrap_16"] <= bands["bootstrap_84"])

        bundle = dict(A=spectrum, A_ved=reference, k=b["k"], energy=b["energy"],
                      eta=b["eta"], coupling=coupling, t=1., omega0=1.,
                      refined_mask=refined, refined_indices=indices,
                      source_per_k=np.where(refined, "refinement", "baseline"),
                      bootstrap_replicates=np.where(refined, 128, 8),
                      spectral_resolution_validated=False,
                      method="41 measured momenta: 10 frozen refined fits and 31 unchanged baseline fits",
                      **bands)
        output = DATA / f"diagmc-updated41-lambda{coupling:.2f}-map.npz"
        np.savez_compressed(output, **bundle)
        denominator = np.trapezoid(reference, b["energy"], axis=2)
        error = np.trapezoid(abs(spectrum - reference), b["energy"], axis=2) / denominator
        record = dict(coupling=coupling, shape=list(spectrum.shape),
                      refined_indices=indices.tolist(), unchanged_baseline_columns=int((~refined).sum()),
                      refined_columns_exact=True, other_columns_unchanged=True,
                      axes_identical=True, pilot_ved_matches_full_reference=True,
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
    for row, b in enumerate(bundles):
        for col, (key, title) in enumerate([("A", "Updated DiagMC"), ("A_ved", "VED reference")]):
            ax = axes[row, col]
            mesh = ax.pcolormesh(edges(b["k"] / np.pi), edges(b["energy"]), b[key][ie].T,
                                 cmap="viridis", norm=norm, shading="flat", rasterized=False)
            ax.set(xlim=(0., 1.), ylim=(b["energy"][0], b["energy"][-1]))
            ax.set_title(rf"$\lambda={b['coupling']:g}$: {title}", pad=15)
            if col == 0:
                ax.set_ylabel(r"$\omega/t$")
                ax.scatter(b["k"] / np.pi, np.full(len(b["k"]), 1.013), marker="|", s=42,
                           linewidths=2., color=np.where(b["refined_mask"], "#d66324", "#b7b7b7"),
                           transform=ax.get_xaxis_transform(), clip_on=False)
            if row == 1:
                ax.set_xlabel(r"$k/\pi$")
            ax.set_xticks([0., .25, .5, .75, 1.], ["0", "1/4", "1/2", "3/4", "1"])
    bar = fig.colorbar(mesh, ax=axes, label=r"$t A_\eta(k,\omega)$", shrink=.91, pad=.02)
    bar.solids.set_rasterized(False)
    scale = "logarithmic color" if logarithmic else "linear color"
    fig.suptitle(rf"Spectral function comparison: $\eta={eta:g}t$ ({scale})" + "\n"
                 "Orange ticks: 10 refined momenta; gray ticks: 31 original momenta", fontsize=12)
    name = f"spectral-comparison-eta{eta:g}-{'log' if logarithmic else 'linear'}"
    for extension in ["svg", "pdf", "png"]:
        metadata = {"Date": None} if extension == "svg" else (
            {"CreationDate": None, "ModDate": None} if extension == "pdf" else {})
        fig.savefig(PLOTS / f"{name}.{extension}", dpi=170, bbox_inches="tight", metadata=metadata)
    plt.close(fig)
    return dict(name=name, eta=eta, scale=scale, vmin=float(norm.vmin), vmax=float(norm.vmax),
                one_color_scale_for_all_four_panels=True,
                files={extension: digest(PLOTS / f"{name}.{extension}") for extension in ["svg", "pdf", "png"]})


def main():
    plt.rcParams.update({"font.size": 11, "svg.fonttype": "none",
                         "svg.hashsalt": "holstein-spectral-comparison-v1"})
    bundles, records = zip(*(assemble(coupling) for coupling in COUPLINGS))
    for key in ["k", "energy", "eta", "refined_mask"]:
        np.testing.assert_array_equal(bundles[0][key], bundles[1][key])
    plots = [plot(bundles, .25), plot(bundles, 1.), plot(bundles, .25, logarithmic=True)]
    provenance = dict(complete=True, script_sha256=digest(Path(__file__)),
                      points=41, refined_points_per_coupling=10, baseline_points_per_coupling=31,
                      eta=bundles[0]["eta"].tolist(), energy_window=[-3.25, 5.5],
                      no_interpolation=True, no_column_normalization=True, no_energy_shifts=True,
                      no_new_sampling_or_fitting=True, ved_used_only_for_comparison=True,
                      spectral_resolution_validated=False, couplings=list(records), plots=plots,
                      bootstrap_note="16-84% empirical ranges: 128 draws at refined points, 8 at original points; representation bias not included")
    (DATA / "spectral-comparison-provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print("Verified 10 replaced and 31 unchanged columns per coupling; created three paired maps.")


if __name__ == "__main__":
    main()
