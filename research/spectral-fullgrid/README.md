# Spectral refinement across the full 41-point grid

Every measured momentum `k = j*pi/40`, `j=0..40`, for `lambda=.25,.5` uses
384 million production steps and the same continuation protocol. The ten
completed pilot momenta are reused byte for byte; the remaining 31 momenta
receive 496 new independent chains and 23.808 billion production steps.
The full data set contains 656 chains and 31.488 billion production steps.
Each chain has 480 blocks of 100,000 steps, after 500,000 warmup steps.

The C++ sampler, conditional estimator and executable are unchanged from the
original 41-point calculation. All 656 seeds are unique and disjoint from the
328 original chains. The plan retains its legacy name `plans/pilot.json` so
that the unchanged inference modules can consume it; its indices span all 41
momenta. `scripts/full_grid_pipeline.py` orchestrates the extension on Kestrel.

## Frozen protocol

The scientific modules are the completed pilot's original source, including
its numerical recovery. Fourteen continuation configurations compare time
resolution/range, covariance cutoffs, frequency grids/support, block length
and a maximum-entropy alternative. Each uses the same held-out chain data,
time bins, scoring range and covariance cutoff. The one-standard-error rule
selects a regularization parameter; a representation replaces the baseline
only when its paired validation gain exceeds one fold standard error. Shared
training folds make these selection heuristics, not significance tests.

All fits use nonnegative weights, finite time-bin kernels, full covariance and
exact unbroadened M0, M1 and M2 constraints. Maximum entropy solves

```
0.5 * ||W (K w - G)||^2 + alpha * sum(w * log(w / m))
subject to w >= 0, sum(w)=1, sum(E*w)=epsilon_k,
sum(E^2*w)=epsilon_k^2+g^2.
```

The implementation uses a convex dual and a uniform spectral-density default.
Convergence is checked explicitly. The independent primal-solver comparison
and the other numerical tests are preserved under `tests/`.

Every coupling/momentum has 128 bootstrap draws (10,496 total). Groups of 16
raw blocks are resampled within each independent chain, and alpha is selected
again within the frozen representation. If a draw fails numerical convergence,
the original recovery reruns the same draw and selected alpha with extended
precision dual accumulation and a 4,000-iteration limit. Initial failed outputs
are retained; successful draws and central fits are preserved.

For draws that still do not converge, `recovery/stable_maxent.py` keeps that
same objective, frequency grid, covariance and selected alpha. It replaces
the clipped-eigenvalue Newton solve with extended-precision pivoted
elimination and evaluates objective differences using a centered cumulant
to avoid cancellation. Three known analytic optima and real draws already
solved by the original method check the fallback before its outputs are
accepted. Only originally failed draws use these replacement results.
`results/bootstrap_recovery.json` records all recoveries, control errors
and source hashes; `results/numerical_recovery-source.json` supplements the
unchanged original inference snapshot. Diagnostic failures and scheduler
timeouts remain in the complete job history.

The 20 old central fits remain frozen, and the 62 new selections are hashed
before new VED comparisons. No VED array is read by the sampling/inference
pipeline. The comparison is computed afterwards on identical energy, momentum
and broadening arrays, with no energy shifts, interpolation or column
normalization. Fits are retained even where comparison errors increase.

## What the uncertainty does and does not establish

The 16–84% and 2.5–97.5% bands describe bootstrap variation with the selected
representation held fixed. Near-optimal representation spread is saved
separately. Neither captures all analytic-continuation bias.

The six problem coupling/momentum points from the pilot retain their original
36 synthetic cases and 864 noise realizations. These tests were not rerun or
expanded to the other momenta. They include known single peaks, doublets and
a broad continuum, with exact moments and measured noise covariance. The
existing unresolved doublets and false peaks remain evidence against general
narrow-peak precision. `spectral_resolution_validated=false` remains in all
new data products. Using the improved protocol at every point does not prove
momentum-grid convergence or accuracy of every narrow feature.

## Reproduction

The full-grid release contains three main archives, one reused-synthetic archive,
four bootstrap archives and 12 raw-chain archives (20 in total). Extract all archives together: every archive has a
`holstein-diagmc/` root. File sizes and SHA-256 hashes are in the archive
manifest. These release assets are separate from the static web payload.

The three main archives together supply the baseline/VED arrays, covariance diagnostics,
frozen selections, source snapshots, plans, logs, audit records and report.
The reused-synthetic archive contains the unchanged six-point resolution tests. The bootstrap archives include final samples and any
original failed outputs. The raw archives contain all 656 chains, including
the 160 verified pilot chains. No external data path or network connection is
needed to reproduce the final report:

```bash
cd holstein-diagmc
python -m pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python scripts/full_grid_report.py --out results/report-replayed
python -m pytest -q tests
make -B test
```

The report checks 2,624 raw-file hashes, the reused files and frozen source,
regenerates central spectra from their weights, and evaluates all 82 points.
It refuses to overwrite an existing report directory. The replay's numerical
fields and spectral arrays can be compared with `results/report/`.
Prebuilt binaries depend on the original compiler/glibc; recompile tests with
a compatible local compiler when needed.

Fresh sampling needs a new output directory, a documented plan and the old
baseline raw data for independent diagnostics and seed-disjointness checks.
The Slurm script records Kestrel account, partition and environment settings;
`results/jobs.json` records actual scheduler settings and accounting. Heavy
sampling and inference run only inside the Slurm allocation. The full-grid
pipeline resumes completed points without changing the frozen source and
verifies reused files again after inference.

Final compact NPZ products have shape `(3,41,701)` with axes `(eta,k,energy)`;
eta is `[.25,.5,1]`. Bootstrap quantiles have shape `(4,3,41,701)` with
probabilities `[.025,.16,.84,.975]`. The paired web maps are generated by
`scripts/build_spectral_comparison.py` in the website repository.
