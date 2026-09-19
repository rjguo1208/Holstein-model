# Updated spectral maps alongside VED

The spectral page now opens with a direct DiagMC–VED comparison. Each row
fixes the coupling and each column fixes the method. Linear maps at eta=.25t
and eta=t, plus a logarithmic map at eta=.25t, use identical momentum and
energy axes and one common absolute-intensity color scale per figure.

The available refinement covers ten of the 41 measured momenta for each
coupling: indices 0, 24, 25, 26, 35, 36, 37, 38, 39 and 40. Only those columns
are replaced by the frozen `selected` spectra. The other 31 columns retain
their original values. Orange ticks identify updated points; gray ticks
identify the original points. No interpolation, refitting, energy shift or
column normalization is applied. This display does not extend the scientific
resolution claims from the [refinement report](spectral-refinement-completed.md).

The former standalone maps and their calculation records remain available
in a collapsed history section. Main heatmaps use PNG previews for loading
speed, with full vector SVG and PDF versions available for download.

## Data and reproduction

Run `python scripts/build_spectral_comparison.py` with the packages pinned
in `research/spectral-refinement/requirements.txt`, then `npm run build`.
The plotting script reads the three existing NPZ products per coupling:

- `diagmc41-lambda*_map.npz`: original 41-point DiagMC.
- `diagmc-refinement-lambda*.npz`: frozen ten-point refined spectra.
- `lambda*_map.npz`: independent 41-point VED.

The new `diagmc-updated41-lambda*-map.npz` files contain paired `A` and `A_ved`
arrays, shape `(3,41,701)` with axes `(eta,k,energy)`. The eta axis is
`[.25,.5,1]`. `refined_mask`, `refined_indices` and `source_per_k` preserve the
per-column provenance. The two empirical bootstrap bands retain their source
values: 128 draws at refined points, eight at original points. The per-point
counts are stored in `bootstrap_replicates`.

The script checks exact equality of the coordinate arrays, the pilot's saved
baseline and VED slices, the ten replacement columns, and the 31 unchanged
columns. It also verifies finite nonnegative spectra and ordered uncertainty
bands. Input/output SHA-256 hashes, all plot color limits and per-point VED
comparison errors are saved in
[`spectral-comparison-provenance.json`](../site/data/spectral-comparison-provenance.json).

## Verification

The build and static checks passed for three pages, 276 LaTeX/MathML formulas,
four TikZ figures and 28 scientific plots. All generated SVG files remain
vector graphics. Browser checks passed for all three pages at widths 1280,
390 and 320 pixels with JavaScript disabled and external requests blocked.
They verified the paired figures appear before the refinement report, the
original maps are initially collapsed, the global page directory works,
images and formulas load, and the page has no horizontal overflow, including
when the historical section is opened.
