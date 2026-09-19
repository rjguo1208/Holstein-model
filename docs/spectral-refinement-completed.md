# Spectral refinement completed on Kestrel

The independent pilot reduced the typical whole-spectrum error, while the
synthetic tests still reject a general narrow-peak precision claim. The full
41-point map remains a labeled baseline; the new results cover ten measured
momenta per coupling, without interpolating the missing momenta.

## Completed computation

- Indices: [0, 24, 25, 26, 35, 36, 37, 38, 39, 40] of the 41-point grid, for lambda=.25 and .5.
- 160 new independent chains; eight chains per point, 48 million production
  steps each, 500,000 warmup steps. Total: 7.68 billion production steps.
- Identical C++ sampler and conditional measurement executable; disjoint seeds,
  complete raw block counts, source/input hashes and zero order-cap attempts verified.
- Sampling job 18680020 and analysis job 18680024 completed with exit code 0:0.
  The three bootstrap numerical recoveries ran in job 18680098 on debug-stdby, verifying the numerical solver fix.
  Unstarted recovery 18680085 was replaced by a native RHEL 8 submission.
  Synthetic tests ran in step 18680024.0 using idle cores during the two slower
  maximum-entropy bootstrap points. The duplicate unstarted jobs 18680031 and
  18680035 were cancelled. Exact accounting is in the published jobs JSON.
- C++ kernel/conditional checks and 17 Python tests passed (16 in the main analysis, plus the extended-precision recovery comparison). New tests compare
  entropy inference with an independent constrained primal solve, check exact
  synthetic moments/integrals, and verify covariance transformation under rebinning.
- All 2,560 bootstrap draws and 864 synthetic draws completed successfully.
  Three bootstrap fits initially failed convergence. The same random draws
  and selected alpha were recovered with a fresh dual start, extended-precision dual accumulation and up to 4,000
  Newton iterations; the objective, constraints and central fits were unchanged.
  Initial failures, successful-draw invariance checks and the recovery job/source
  hashes are retained in the bootstrap recovery record and main archive.

## Measured changes

All numbers below compare the same ten momenta, energy window [-3.25,5.5] and
eta=.25t against unchanged VED arrays. They include the complete spectrum,
including the ground-state contribution, and are not pointwise error bounds.

| lambda | Original median | 4x sampling only median | Refined median | Original maximum | Refined maximum |
|---|---|---|---|---|---|
| 0.25 | 21.9% | 26.1% | 14.3% | 37.2% | 16.6% |
| 0.5 | 40.4% | 39.6% | 11.0% | 47.1% | 24.2% |

Eighteen of the twenty coupling/momentum cases improve. At lambda=.25,
k/pi=.6 and .925 have small increases; these fits remain unchanged after
evaluation. Sampling alone does not consistently improve the reconstructed
spectrum. The median pointwise standard-error ratio for imaginary-time G is
0.5096 and
0.5109, close to the 1/2 expectation.
The block-16/block-4 median-error ratios range from .95 to 1.03.

The 14 continuation configurations share a fixed held-out scoring convention.
They vary time resolution/range, covariance cutoff, spectral grid/support,
block length and maximum entropy. The one-standard-error rules are selection
heuristics, not significance tests. Full details and equations are in the
[scientific source README](../research/spectral-refinement/README.md).

The 20 central fits and source/input records were frozen and hashed before
the first VED evaluation. The bootstrap reselects alpha within each frozen
representation, using groups of 16 raw blocks within each chain. It does not
resample representation choice. Near-optimal representation spread is saved
separately. Independent old-data prediction scores are also in the report,
and are not used to tune the new fits.

## Resolution limitations

Six problem points have six known synthetic spectra each: one narrow peak,
four doublet separations, and a broad continuum. Each case has 24 independent
Gaussian noise realizations with the measured absolute covariance; covariance
is treated as known in these conditional tests. All spectra obey M0, M1 and M2.

Doublet recovery requires maxima within separation/4 of their true positions,
at least 5% prominence, and at least 10% valley contrast. Figures mark cases
where the broadened truth itself fails the criterion as n/a. Failures remain
in the denominator. At k/pi=.9 and eta=.1t, neither coupling recovers the .6t
doublet in any of 24 draws. Broad continua can acquire spurious peaks.
Thus `spectral_resolution_validated=false` remains appropriate. A 24-draw
empirical recovery rate is not a general resolution or confidence guarantee.

## Reproduction and artifacts

- [Pilot summary](../site/data/diagmc-refinement-summary.json),
  [statistical diagnostics](../site/data/diagmc-refinement-diagnostics.json),
  [Slurm accounting](../site/data/diagmc-refinement-jobs.json).
- [Selection freeze](../site/data/diagmc-refinement-selection-freeze.json) and
  [raw data audit](../site/data/diagmc-refinement-raw-audit.json).
- [Package sizes and SHA-256](../site/data/diagmc-refinement-archives.json).
  Extract the main archive, bootstrap-draw archive and the two raw-chain archives into one directory.
  The main archive includes baseline/reference arrays, so the final evaluation
  can run without an external repository or network connection.
- [Archive reproduction verification](../site/data/diagmc-refinement-archive-validation.json).
  The published source is under `research/spectral-refinement/`.

Each pilot NPZ contains baseline, sampling-only, refined and reference spectra
with shape (3,10,701), axes (eta,k,energy), eta=[.25,.5,1]. Bootstrap quantiles
have shape (4,3,10,701), with probabilities [.025,.16,.84,.975]. The original
41-point downloads and plots retain their names and values.

The page provides five new SVG/PDF/PNG figures. Fresh extraction of all four
archives verified 640 raw-file hashes, 37 reviewable source files, 20 original
inference source files and all 320 duplicate entries. Replaying the report
from packaged inputs reproduced its numerical fields exactly. The 253
previously successful draws at the two repaired points remained identical.
All 17 Python tests passed, as did the C++ kernel tests after recompiling
the archived source with the local compiler. Prebuilt executables depend on
the build system's glibc; use `make -B test` to rebuild on another system.

`npm run build` and `npm run check` passed: three pages, 271 LaTeX/MathML
expressions, four TikZ figures and 25 scientific plots, including local links,
anchors and fonts. Browser checks passed for all three pages at widths 1280,
390 and 320 pixels with JavaScript disabled and external requests blocked.
The checks covered the page directory, formulas, image loading, downloads and
horizontal page overflow.
