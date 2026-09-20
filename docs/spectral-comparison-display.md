# Full-grid spectral maps alongside VED

The spectral page opens with a direct DiagMC–VED comparison. Each row fixes
the coupling and each column fixes the method. Linear maps at eta=.25t and
eta=t, plus a logarithmic map at eta=.25t, use identical momentum and energy
axes and one common absolute-intensity color scale per figure.

All 41 measured momenta now use the improved DiagMC results for each coupling.
The ten completed pilot momenta are reused unchanged, and the other 31 are
independently sampled and analyzed with the same protocol. Every point uses
four times the original production sampling and 128 bootstrap draws. Details,
measured errors and limitations are in the
[full-grid report](full-grid-refinement-completed.md).

The plotting script performs no additional fit, interpolation, energy shift
or column normalization. The former mixed 10/31 display has been replaced.
The original standalone maps and ten-point pilot remain available in collapsed
history sections. Main heatmaps use PNG previews for loading speed, with full
vector SVG and PDF versions available for download.

## Data and reproduction

Run `python scripts/build_spectral_comparison.py` with the packages pinned in
`research/spectral-fullgrid/requirements.txt`, then `npm run build`.
The plotting script reads three NPZ products per coupling:

- `diagmc41-lambda*_map.npz`: original 41-point DiagMC, for provenance checks.
- `diagmc-fullgrid-lambda*.npz`: frozen full-grid refined spectra and quantiles.
- `lambda*_map.npz`: independent 41-point VED.

The `diagmc-updated41-lambda*-map.npz` files contain paired `A` and `A_ved`
arrays, shape `(3,41,701)` with axes `(eta,k,energy)`. The eta axis is
`[.25,.5,1]`. Every `refined_mask` entry is true, `refined_indices` is 0..40,
`source_per_k` is uniformly `refinement`, and `bootstrap_replicates` is 128
at every point. The 16–84% bands come from those same 128-draw calculations.

The script verifies exact equality of coordinate arrays, saved baseline and
VED arrays, and all 41 refined columns. It checks finite nonnegative spectra,
ordered uncertainty bands, completion of all 10,496 bootstrap draws and
agreement of every plotted VED comparison error with the scientific report.
Input/output SHA-256 hashes, plot color limits and per-point errors are saved
in [the display provenance](../site/data/spectral-comparison-provenance.json).

All points use the improved protocol, but `spectral_resolution_validated`
remains false. The existing synthetic tests cover six pilot coupling/momentum
points and do not certify every narrow feature across the full map.

## Completed validation

`npm run build` and `npm run check` passed: 3 pages,
279 LaTeX expressions with MathML, 4 TikZ figures and
31 scientific SVG/PDF plots. Chromium checked all three pages at
1280, 390 and 320 pixels with JavaScript disabled and external requests
blocked. Navigation, formulas, images, downloads and both expanded history
sections passed; there was no horizontal page overflow. Every full-grid
column and its 128-draw uncertainty bands were checked against the scientific
report. The complete Pages directory is 991,246,655 bytes.

All 20 archives passed fresh extraction and hash checks, including 2,624 raw
files. The extracted report reproduced every numerical field and NPZ array
exactly. Its 17 Python tests and C++ kernel tests passed; the latter used a
fresh compilation from archived source. Original maps and pilot records
remain byte-for-byte unchanged in the tracked repository.
