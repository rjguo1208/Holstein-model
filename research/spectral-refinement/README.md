# Spectral refinement pilot

This directory contains the reviewable scientific source for the Kestrel pilot
at `k/pi = [0, .6, .625, .65, .875, .9, .925, .95, .975, 1]`, for
`lambda = .25, .5`. The C++ sampler and conditional measurement are unchanged
from the completed 41-point calculation. Results and large raw files are linked
from the [spectral page](https://rjguo1208.github.io/Holstein-model/spectral-map.html#refinement).

Each point has eight new independent chains, each with 480 blocks of 100,000
production steps and 500,000 warmup steps. Thus each coupling/momentum receives
384 million production steps, four times the original 96 million. All seeds
are disjoint from the original 328 chains. The complete pilot uses 7.68 billion
new production steps. `plans/pilot.json` records the original paths, settings
and source/executable hashes.

## What changes in inference

`refinement.py` compares 14 configurations: the original frequency grid and
time range, finer time bins, shorter/longer time ranges, two covariance cutoffs,
coarser/finer/shifted energy grids, two support cutoffs, two block sizes, and
maximum entropy with exact unbroadened moments. All are scored against the
same held-out-chain time bins and covariance cutoff. VED is absent from
sampling, model selection, bootstrap and synthetic inference.

Within each configuration, a one-standard-error heuristic selects the largest
regularization parameter consistent with the minimum validation score. A
different configuration replaces the baseline when its paired validation gain
exceeds one fold standard error. Shared training folds mean these rules are
heuristics, not calibrated significance tests. The final 20 selected fits are
frozen and hashed before VED is loaded for evaluation.

The maximum-entropy implementation minimizes

```
0.5 * ||W (K w - G)||^2 + alpha * sum(w * log(w / m))
subject to w >= 0, sum(w)=1, sum(E*w)=epsilon_k,
sum(E^2*w)=epsilon_k^2+g^2.
```

It uses a small convex dual, normalized exponential weights, exact moment
multipliers, projection of fixed moment directions and a preconditioned Newton
solve. The default `m` is a uniform spectral density. Numerical convergence is
checked; unconverged fits cannot win cross validation. A test compares the
result with a separate constrained primal SLSQP solve.
Kernel singular values below 1e-12 of the largest projected value are treated
as numerically null; reported residuals still use the original bin kernel.

`pilot_bootstrap.py` attempts 128 draws per point, resampling groups of 16 raw
blocks within each independent chain and reselecting alpha. The representation
is held fixed. Near-optimal representation spread is saved separately; neither
band captures every possible analytic-continuation bias.

Three bootstrap draws initially failed numerical convergence. The recovery
stage repeats exactly those random draws and their original cross-validation
rule, keeps the selected alpha, and retries the same entropy objective from a
fresh dual start and extended-precision dual accumulation, with up to 4,000 iterations. It preserves all successful draws
byte for byte at the array level and never changes the frozen central fits.
The failed outputs and recovery provenance are included in the archive.

`resolution.py` constructs known single peaks, doublets and a broad continuum
with exactly the Hamiltonian's first three moments. Six problem points each
have six synthetic cases and 24 independent noise realizations. The eight
synthetic chain means use the measured absolute covariance, with covariance
treated as known. Recovery requires two local maxima near the true positions,
at least 5% prominence and 10% valley contrast. Failures count as non-recovery.
These conditional tests do not certify every physical spectrum or momentum.

## Validation and reproduction

Python 3.12 and the pinned packages in `requirements.txt` were used. For the
small numerical tests, run:

```bash
python -m pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q tests
make test
```

The completed download comprises a main archive, a bootstrap-draw archive and
two raw-chain archives. Extract all four into one directory; each contains a `holstein-diagmc/` root.
The main archive includes the baseline/reference NPZ arrays, diagnostics,
frozen fits, synthetic results, figures and hashes. The separate bootstrap
archive contains the individual resampled spectra. To rerun
the final evaluation entirely from those files, preserve the published output
directory first:

```bash
cd holstein-diagmc
mv results/report results/report-published
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/pilot_report.py
```

The report prefers packaged inputs and raw chains, checks frozen selections
and input hashes, regenerates central spectra from their weights, and compares
all three stages against identical VED grids. Repeating the original
82-point diagnostic additionally requires the previous 41-point raw archives.

For a fresh computation, use an empty results directory and update the paths
in a new run plan. `sample_pilot.py` checks its executable against the recorded
hash; a new build needs its own documented compiler/executable provenance.
The Slurm files record Kestrel's account, standby partition and Python paths.
The actual run consolidated diagnostics, inference and uncertainty into job
18680024; a second step used idle cores for synthetic tests. The two unstarted
duplicate jobs were cancelled. Heavy sampling/inference ran on compute nodes.

Files refuse to overwrite completed suites. `pilot_resolution.py` uses a file
lock so the serial driver can join the synthetic stage without recomputation.

## Background

- [Koch, analytic continuation lecture notes](https://cond-mat.de/events/correl22/manuscripts/koch.pdf)
- [Shao and Sandvik, stochastic analytic continuation](https://arxiv.org/abs/2202.09870)
- [Schumm, Yang and Sandvik, cross validation](https://arxiv.org/abs/2406.06763)

The algorithm above is a maximum-entropy comparison, not an implementation of
the cited stochastic analytic continuation sampler.
