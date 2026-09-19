# 41-point DiagMC calculation completed on Kestrel — 2026-09-19

The pending denser momentum calculation has completed on Kestrel. All three
final Slurm jobs are COMPLETED with exit code 0:0, and all 328 raw chains and
82 coupling/momentum reconstructions have passed the output-integrity audit.
The spectral page now displays the measured 41-point maps. Frequency-resolution
validation remains open; completing the grid does not establish convergence of
the fine real-frequency spectrum or of momentum discretization itself.

## Completed jobs

| Job | Stage | Elapsed | Requested / allocated CPUs | Node |
|---|---|---|---|---|
| 18674188 | preflight | 00:00:14 | 4 / 104 | x1005c0s0b0n0 |
| 18674189 | sampling | 00:00:20 | 104 / 104 | x1005c0s0b0n0 |
| 18674190 | continuation_and_audit | 00:04:00 | 1 / 104 | x1005c0s0b0n0 |

Submission host: kl5 (RHEL 9); account sipv, QoS standby, partition short-stdby.
The partition allocates whole 104-core nodes. Sampling used up to 13 independent
momentum folders concurrently, with eight single-threaded chains per folder.
The original continuation remains single-threaded. The JSON accounting records
include exact start/end times in UTC and total CPU time.

An initial environment preflight (18673763) failed before compilation because
the RHEL 9 module hierarchy could not load the RHEL 8 Python module name.
Submission and runtime setup were corrected, and preflight 18674188 passed.
The unstarted dependent jobs 18673765/18673766 were cancelled and replaced.
Historical Anvil jobs 20829123/20829124 were neither queried nor cancelled here;
the old maintenance message remains only in the historical submission record.

Working directory:
`/scratch/rjguo/holstein-diagmc-kestrel-20260919/holstein-diagmc`.

## Sampling and validation

- Grid: k=j*pi/40, j=0,...,40, exactly the stored VED coordinates.
- t=omega0=1, lambda=g^2/(2*t*omega0)=0.25 and 0.5.
- Four chains per coupling and momentum; 240 blocks of 100,000 production
  steps per chain, with 500,000 warmup steps for each newly sampled chain.
- mu=-2.8, tau_max=12, 96 time bins, maximum order 96, rb_thin=100.
- Nine matching momenta reuse 72 chains byte for byte; the other 32 momenta
  contain 256 new chains with seed offset 300,000,000. New seeds are disjoint
  from all 136 old map seeds. No old point was interpolated or relabeled.
- New production steps: 6,144,000,000; final total: 7,872,000,000.
- C++ diagram/conditional-estimator checks and 12 Python tests passed.
  Free-particle, atomic-limit and first-order-truncation comparisons passed for
  both the direct and conditional measurements (six analytic comparisons).
- All 328 completion flags, unique seeds, model parameters, raw block counts,
  source snapshots and continuation input hashes passed the audit. Reused file
  hashes match their originals. Order-cap attempts: 0.
- All 32 new points used the same tested executable. Core scientific sources,
  the continuation grid/method, regularization scan and eight block-bootstrap
  replicates retain the submitted design; only cluster startup and scheduling
  were adapted.

Python 3.12.5, GCC 12.3.0, NumPy 2.3.5, SciPy 1.16.3, Matplotlib 3.10.6 and
pytest 8.4.2 were used. Full package versions and executable hashes are in
`site/data/diagmc-map41-environment.json` and the archived preflight log.

## Full-spectrum comparison at all 41 momenta

The metric is the relative integral of |A_DiagMC-A_VED| over [-3.25,5.5], divided
by the VED integral in that same window, with identical k, energy and broadening
arrays. VED spectra enter only this comparison, after continuation and model
selection. This is the full spectrum, including the ground-state contribution,
not the separate incoherent-spectrum error reported on the k=0 results page.

| lambda | eta=.25t range | eta=.5t range | eta=t range | Largest held-out score |
|---|---|---|---|---|
| 0.25 | 4.43%–39.88% | 1.68%–15.28% | 0.39%–3.50% | 107.66 |
| 0.5 | 7.25%–53.74% | 2.95%–23.43% | 0.93%–5.93% | 47.92 |

Maximum unbroadened moment residual: 9.353e-10. Arrays are finite and
nonnegative; bootstrap quantile arrays are finite and correctly ordered. The
selected regularization parameter lies at a scan boundary at eight points for
lambda=.25 and ten points for lambda=.5; this is retained as an open diagnostic.
The central spectra at the nine reused momenta match the earlier reconstructions
to maximum relative L1 differences of 2.79e-10 and 2.65e-9, respectively. Their
bootstrap replicas use the new grid-index seeds from the unchanged algorithm.
The detailed reuse comparison is in `diagmc-map41-reused-spectra.json`.
The eight-replicate 16–84% ranges are empirical block-bootstrap ranges and do not
include all spectral-parametrization bias. Large held-out scores and the fine
spectrum differences remain visible; `spectral_resolution_validated=false`.

## Data and reproduction

New NPZ arrays have shape (3,41,701), axes (eta,k,energy), and eta=[.25,.5,1].
Downloads use the `diagmc41-`/`diagmc-map41-` prefixes; the earlier 17-point
downloads and Anvil submission snapshots are retained separately.

For a complete local audit, unpack `diagmc-spectral-map-01.zip` (the unchanged
VED reference) first, followed by `diagmc-map41-kestrel-completed-01.zip` and
the four `diagmc-map41-raw-lambda*-k*.zip` archives into one directory. These
five Kestrel archives contain all 328 chains, including the reused chains;
the earlier raw-chain archives are not needed for this audit. ZIP sizes,
file counts and SHA-256 values are in `diagmc-map41-kestrel-archives.json`.

From the extracted `holstein-diagmc` root, set HOLSTEIN_PYTHON to a compatible
Python environment with the recorded numerical-library versions, then run:

```bash
bash scripts/python.sh scripts/audit_dense_map.py \
  --data results/map_data_41_01 --map results/map_diagmc_41_01 \
  --out results/map_audit_41_recheck
```

To rerun the spectral continuation, use a fresh output directory:

```bash
bash scripts/python.sh scripts/continue_map.py \
  --data results/map_data_41_01 --out results/map_diagmc_41_repeat \
  --points 41 --bootstrap 8
```

Sampling from scratch uses the legacy inputs listed in
`results/kestrel_input_archives.json` and the current source files, in a fresh
working directory without the completed `map_data_41_01` output. The Kestrel
preflight prepares the nine reused points, then sampling fills the 32 remaining
folders and the dependent finish job performs continuation and audit. Submit
from a login node matching the compute OS. Each stage refuses to overwrite its
output, so preserve existing results and choose a new directory for a rerun.

## Publication checks

- All five new ZIP files passed integrity and SHA-256 checks. A fresh extraction
  with the VED reference reran the full 328-chain, 41-point audit and reproduced
  every numerical audit value exactly. See the
  [archive verification record](../site/data/diagmc-map41-archive-validation.json).
- Static build and site checks passed: three pages, 266 LaTeX expressions with
  MathML, four TikZ figures and 20 scientific SVG/PDF plots; local links, anchors
  and fonts are valid.
- Chromium checks passed at 1280, 390 and 320px with JavaScript disabled and
  external requests blocked: no missing images, math errors or body overflow.
  Global navigation, current-page markers and the worked sampling example
  remain functional. The fine and coarse PNG exports were also inspected.
