# Single on-site defect: implementation and verification

The published extension changes only the electron potential to `-U n0` in the
1D, zero-temperature, single-electron Holstein model. Parameters are
`t=omega=1`, `lambda=g²/(2t omega)=0.25,0.5`, and `U=0.5,1,2`.

## Repository synchronization

Before publication, the repository was fetched and fast-forwarded from
`5070b21` to `5a20e613b48c9bcd672e1ed21fd59aea7d046979`. This includes the worked
sampling example, Kestrel 41-point results, full-grid refinement and bilingual
maintenance rules. These existing results and data are preserved. Changes to
existing pages add navigation and the defect entry; the prior spectral body
and all existing scientific assets retain their remote contents.

The new page is available in Chinese and English. Its 90 new translation keys
cover prose, metadata, table labels, captions, alternative text and explanatory
equation text. Both languages share numerical tables, mathematical expressions,
scientific figures and downloadable data, with the existing language switch.

## Scientific implementation

Reviewable source is in [research/holstein-defect](../research/holstein-defect/README.md).
The continuous-time direct-space bare expansion samples electron hops and
phonon contractions. The defect potential enters exactly through residence
time. Electron paths live on the infinite chain; the reporting radius is not
a physical boundary. The code permits crossed and nested phonon contractions.

Normalization uses the zero-hop, zero-phonon sector. Conditional external-time
integration supplies bin averages; correlated long-time fits yield energy and
injection residue. Midpoint projectors independently measure density, local
phonon-vacuum pole weights, total zero-phonon weight, phonon number and radius.
Binding energies subtract independent clean and defect DiagMC energies. VED
ground energies select the auxiliary importance shift only; they do not fix
the measured energies or weights.

The independent VED basis varies relative-cloud generations and the electron's
absolute-position radius separately. Ground checks use generations 8,10,12
and enlarge the radius by 8 sites. Local spectra compare generations 10→12,
radius 24→40 and Lanczos 200→400 steps at all nine measured positions.

## Completed sampling and validation

All jobs completed on the Anvil `highmem` partition, node `b000`, exit `0:0`:

| Stage | Job | Chains | Production steps | Elapsed |
|---|---:|---:|---:|---:|
| Pilot and analytic checks | 20881037 | 28 | 128,000,000 | 00:00:31 |
| Six ground-state points | 20881090 | 48 | 1,536,000,000 | 00:05:52 |
| Local spectra and references | 20881164 | 72 | 1,152,000,000 | 00:05:55 |

Total: 148 chains and 2.816 billion production steps, excluding warm-up.
Recorded parameters, block counts, unique seeds, source/executable/input hashes
and zero production order-cap attempts pass the report audit. The five sampled
analytic checks cover free hopping, a noninteracting impurity, the atomic limit,
one phonon line and one hop pair. A further finite-time atomic check validates
the midpoint phonon-vacuum probability and phonon number using the pilot data.

C++ tests cover detailed balance, locality, residence times, crossed topology,
time Jacobians and conditional energy quadrature. The 17 Python tests include
independent-reference limits, Hermiticity, local spectral moments,
Hellmann–Feynman derivatives, normalization and synthetic local continuation.

## Results and limits

All six ground energies differ from VED by at most 1.108 estimated statistical
standard errors. At `U=1`:

| lambda | DiagMC Edef | DiagMC binding energy | p0 | Zb |
|---|---|---|---|---|
| 0.25 | -2.50638 ± 0.00123 | 0.27763 ± 0.00124 | 0.5206 ± 0.0033 | 0.8102 ± 0.0014 |
| 0.5 | -2.80218 ± 0.00162 | 0.33234 ± 0.00163 | 0.6121 ± 0.0040 | 0.6161 ± 0.0026 |

Errors are conservative estimates across several block lengths and
independent-chain jackknifing. Three projection windows and chain estimates
are retained. These are not rigorous bounds on all residual projection bias.
Density, local pole weights and total zero-phonon weight are distinct quantities.
Reported localization lengths are finite-window fits; tail-window dependence
prevents treating them as precise asymptotic lengths.

For `U=1`, local spectra are measured at sites 0..8 for each coupling. Reflection
symmetry supplies the remaining eight displayed columns, with no interpolation.
Map arrays have axes `(eta,site,energy)` and shape `(3,17,801)`. Plot labels use
English and mathematical notation; SVG and PDF remain vector graphics.

Positive continuation uses local unbroadened moments, covariance, exact time-bin
kernels and held-out-chain selection. VED spectra do not enter reconstruction.
The eight block bootstrap resamples do not include every support or
regularization uncertainty. Whole-spectrum relative L1 differences from VED:

| lambda | eta=0.25t | eta=t |
|---|---|---|
| 0.25 | 11.4–44.0% | 1.3–5.2% |
| 0.5 | 12.2–36.7% | 1.3–6.7% |

Fine spectral accuracy remains unvalidated. Some narrow peaks lack support in
VED. Coarser broadening only supports broad-envelope comparisons and does not
establish narrow peaks, sideband widths or lifetimes.

## Distribution

Five immutable ZIP archives are hosted in the
[single-defect release](https://github.com/rjguo1208/Holstein-model/releases/tag/single-defect-20260923),
following the existing full-grid distribution pattern. They share the
`holstein-diagmc/` root and contain source, immutable run snapshots, exact input
manifests, references, all 148 raw chains and the propagator block arrays.
The main archive includes clean-system energy summaries used for binding; their
original full raw data remain in earlier published datasets.

Each ZIP passes an integrity check, and byte counts and SHA-256 digests are
listed in [defect-archives.json](../site/data/defect-archives.json).
After extracting all five archives into a fresh directory, both C++ test
executables and all 17 Python tests pass. The reporting script is rerun there;
all numerical fields and audit results reproduce exactly, excluding only the
report timestamp and relocated clean-source paths. Details are in
[defect-archive-validation.json](../site/data/defect-archive-validation.json).
The per-job accounting and numerical audit are provided separately in
[defect-jobs.json](../site/data/defect-jobs.json) and
[defect-summary.json](../site/data/defect-summary.json).

## Website checks

`npm run build` and `npm run check` pass: eight bilingual pages, 652 rendered
LaTeX expressions with matching MathML, four existing TikZ figures and 35
vector scientific figures. The checks also verify paired language navigation,
identical mathematical expressions/numerical tables/resources, translations,
local links, anchors and local fonts. Scientific verification and successful
page generation are reported separately.

Chromium also checks all four pages in both languages at 1280, 390 and 320
pixels: 24 page/viewport combinations, with JavaScript disabled and external
requests blocked. All equations and images load, documents have no horizontal
overflow, and language switching and navigation preserve the
`/Holstein-model/` deployment prefix. Desktop and mobile screenshots of the
new page and the local spectra were visually reviewed. The recorded browser
results are in [defect-browser-check.json](defect-browser-check.json).
