# Full-grid spectral refinement completed on Kestrel

Both couplings now use the same improved protocol at every one of the 41
measured momenta, k=j*pi/40. The public DiagMC–VED maps contain no remaining
baseline columns. All central fits are retained, including cases where the
independent VED comparison worsens. Narrow-peak accuracy remains unvalidated.

## Completed scope

- Couplings .25 and .5, indices 0..40. The ten pilot momenta were reused after
  verifying 825 file hashes; 31 more momenta were computed.
- New: 496 chains and 23.808 billion production steps. Reused: 160 chains and
  7.68 billion steps. Complete grid: 656 chains and 31.488 billion steps.
- Every coupling/momentum uses eight chains, each with 480 blocks of 100,000
  production steps, after 500,000 warmup steps. This is four times the original
  96 million production steps per point.
- Job outcomes: 18680984: TIMEOUT (exit 0:0); 18681467: FAILED (exit 1:0); 18681535: COMPLETED (exit 0:0); 18681588: FAILED (exit 1:0); 18681618: COMPLETED (exit 0:0). Completed point outputs were preserved after the time limit and unfinished point calculations were resumed in the final successful job.
  Production sampling and inference ran on compute nodes under the existing
  standby allocation. Exact scheduler settings, start/end times and resources
  are in the jobs JSON.
- The unchanged sampler executable, source hashes, unique seeds disjoint from
  the original 328 chains, complete block counts and zero order-cap attempts
  were checked. The raw audit contains 2,624 file hashes.
- Every one of the 82 coupling/momentum pairs has 128 successful bootstrap
  draws: 10,496 total. 21 numerical failures were recovered from their
  original random samples and alpha selections, including the pilot's three
  recoveries. Successful draws and central selections remained unchanged.
- 2 draws required the numerical fallback in
  `recovery/stable_maxent.py`: the identical MaxEnt dual is solved with
  extended-precision pivoted Newton systems and centered objective differences.
  Three known analytic optima and converged real-draw controls passed. Its
  source and validation records are separate from the unchanged inference snapshot.
- C++ kernel tests and 17 Python tests passed in the compute-job preflight.
  The packaged source and all numerical report fields are independently checked
  again after fresh extraction; see the archive-validation JSON.

## Whole-spectrum comparison

The following values use all 41 momenta at eta=.25t and energy window
[-3.25,5.5], with identical VED grids and absolute intensities. They include
the ground-state contribution and are not pointwise error bounds.

| lambda | Original median | 4x sampling only median | Refined median | Original maximum | Refined maximum |
|---|---|---|---|---|---|
| 0.25 | 16.3% | 14.8% | 10.2% | 39.9% | 17.9% |
| 0.5 | 23.3% | 17.2% | 10.1% | 53.7% | 24.2% |

70 of 82 cases improve and 12 worsen. Worsened (lambda,index) pairs:
`[(0.25, 3), (0.25, 6), (0.25, 10), (0.25, 13), (0.25, 14), (0.25, 23), (0.25, 24), (0.25, 37), (0.5, 9), (0.5, 11), (0.5, 29), (0.5, 30)]`. None is replaced using VED. The imaginary-time standard-error
median ratios are 0.5007 and 0.4971.
The maximum unbroadened moment residual is 8.85e-08.

The predeclared 14 continuation candidates and one-standard-error heuristics
are unchanged from the pilot. All candidates use common held-out scoring;
reference data never enter sampling or model selection. Source and central
fits are frozen before new comparison. Bootstrap resampling holds the spectral
representation fixed and reselects alpha; representation sensitivity is saved
separately. The new example curves use indices 8, 20 and 32, chosen before the
new comparisons. Full methodological details are in the source README.

## Resolution limits

The original 36 synthetic cases (864 draws) are reused byte for byte, covering
only the two couplings at indices 25, 36 and 39. There are no additional
synthetic draws and no new full-grid resolution certificate. The unresolved
.6t doublets at k/pi=.9, eta=.1t (0/24 for each coupling) and the possible
false peaks remain relevant. `spectral_resolution_validated=false` is retained.
Bootstrap bands do not include every analytic-continuation model bias, and
this extension does not demonstrate momentum-grid convergence.

## Data and reproducibility

- [Full report](../site/data/diagmc-fullgrid-summary.json), [jobs](../site/data/diagmc-fullgrid-jobs.json),
  [raw audit](../site/data/diagmc-fullgrid-raw-audit.json), [frozen fits](../site/data/diagmc-fullgrid-selection-freeze.json).
- [Source snapshot](../site/data/diagmc-fullgrid-inference-source.json),
  [reused-file hashes](../site/data/diagmc-fullgrid-reuse-manifest.json),
  [reviewable scientific source](../research/spectral-fullgrid/README.md).
- [Complete release](https://github.com/rjguo1208/Holstein-model/releases/tag/full-grid-refinement-20260919): three main archives, one reused-synthetic archive, four bootstrap
  archives and 12 raw-chain archives (20 total). Extract all into the same directory. The final report
  runs offline using packaged original/reference arrays and diagnostics.
- [Archive sizes and SHA-256](../site/data/diagmc-fullgrid-archives.json) and
  [independent extraction/replay checks](../site/data/diagmc-fullgrid-archive-validation.json).

The full-grid NPZs contain all analysis stages, empirical quantiles and
representation spread. The paired web NPZs contain A, A_ved and 16–84% bands;
all 41 `refined_mask` entries are true and all bootstrap counts are 128.
The three paired heatmaps share one color scale across four panels. No
momentum interpolation, column normalization or energy shift is applied.
Original maps and the ten-point pilot remain available in collapsed history.

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
