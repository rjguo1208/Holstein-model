# Dense local spectra for a single on-site defect

This extends the completed 1D, zero-temperature, one-electron calculation with
`t = omega = U = 1`, `lambda = 0.25, 0.5`. The original ground-state results and
published archives remain historical records. The new run samples 17 injection
sites `0,...,16` independently. Reflection gives 33 displayed sites `-16,...,16`.
More sites extend spatial coverage; the lattice spacing is still one.

## Acquisition and performance

`defect_dense_sample.py` launches 16 independent C++ processes on highmem.
Each of the 34 parameters has eight chains, each with 500,000 warm-up steps,
320 blocks of 100,000 production steps, and thinning 100. There are 8.704 billion
new production steps. The maximum imaginary time remains 12; raw bins increase
from 96 to 192. Each position has four times the previous production statistics.
Seeds, input parameters, source files and executable hashes are saved. Production
must finish without attempted arc/hop cap violations.

The sampler reuses state and position-vector storage; hopping removal samples
ranks without allocating temporary lists. Conditional time-bin integrals reuse
log-gamma evaluation and incomplete-gamma recurrences. `--observables green`
omits unused projector/energy measurements and profile output. The default
`--observables all` preserves the previous interface. Four fixed-seed cases
verify bitwise-identical Green-function blocks and Monte Carlo trajectories
between baseline and optimized code. Other observables agree to rounding error.
The measured green-only speedup is about 1.64–1.68 on one highmem CPU; it is a
benchmark of these cases, not a universal performance guarantee.

## Continuation

The positive local spectral measure obeys
`M0 = 1`, `M1 = V_i`, `M2 = V_i^2 + 2 t^2 + g^2`, where `V_i = -U delta_i0`.
Its improved analytic lower support bound is `-sqrt(4 t^2 + U^2) - g^2/omega`.
Completing the phonon square leaves the free electron plus defect Hamiltonian,
whose exact lowest energy is `-sqrt(4 t^2 + U^2)`; the remaining squares are
positive. This tighter bound improves conditioning without a VED energy prior.
The legacy Tikhonov candidate retains its older bound `-2t-U-g^2/omega`. There is no imposed
pole energy, residue or VED spectral prior. The auxiliary sampling shift mu is
inherited from the preceding ground calculation and appears in the kernel only.

`defect_dense.py` compares nine predeclared configurations using common validation
bins (`dt = 0.25`, `tau <= 8`) and full retained covariance. Eight chains form four
folds, each holding out two chains. Within each configuration, the largest alpha
within one fold-standard-error of the minimum mean prediction score is selected.
The representation with the lowest selected score is then frozen. Because folds
share training observations, this is a selection heuristic, not a confidence
interval or a resolution guarantee.

The primary maximum-entropy representation uses 641 energies, upper support 10,
`dt = 0.125`, `tau <= 8`, relative-covariance cutoff 1e-8 and 400,000-step groups.
Variants use `dt = 0.0625`, `tau <= 12`, 321 energies, upper support 14,
800,000-step groups, or covariance cutoffs 1e-6 / 1e-10. The previous 321-energy
Tikhonov representation is included with the new common selection rule. The
selected grid need not be the finest grid. Candidate spectra and scores are saved.

The exact-moment MaxEnt dual projects out moment directions and performs a small
kernel SVD once per alpha scan. Cached bin-integral kernels and grid matrices are
reused. Long-double Newton solves, centered objective differences, adaptive
linear-solve damping and a smooth-to-sharp homotopy stabilize the long-time
kernel. Damping changes the search direction, not the objective. Every retained
fit must satisfy explicit gradient and moment-residual convergence criteria.
This addresses numerical failures discovered during the pilot; failed trials
are not reported as physical spectra.

64 independent block-bootstrap draws per position resample inside each chain
and reselect alpha with the representation fixed. The pointwise 16–84% ranges
are conditional empirical ranges. Finite grids, support, covariance truncation,
and representation selection remain systematic uncertainties. Blocking and
whole-chain error estimates are saved separately. Process-level parallelism is
used across positions; BLAS/OpenMP threads are restricted to one per process.

## References and resolution diagnostics

Only after all continuation choices have been frozen is VED run at all 17 sites.
The reference independently varies phonon-cloud generations 10/12, electron
radius 40/64, and Lanczos recursion 200/400. The finest combination is used for
comparison. These are convergence differences, not rigorous infinite-basis
error bounds. Local spectral moments are checked independently.

`defect_dense_resolution.py` uses prescribed single/doublet spectra with the
exact local moments and relative covariance measured in the old/new Monte Carlo
data. Both noise levels use the same 641-energy MaxEnt representation, common
time bins, and chain-held-out alpha selection. Sixteen Gaussian-noise experiments
per case measure peak recovery and false splitting. The finite-sample diagnostic
does not certify that a particular physical sideband is resolved.

All displayed spectra use 1601 energies from -4.5 to 5.5 and Lorentzian half-widths
0.15, 0.25, 0.5, 1.0. Smaller broadening is exploratory. Plotting a denser grid and
using narrower eta do not by themselves improve inferential resolution.
Absolute weights and a common vacuum energy origin are retained.

## Plotting and bare-electron guides

The local map has site on the horizontal axis. Dashed horizontal lines indicate
the `g=0` continuum edges `+-2t` and the impurity bound energy
`-sqrt(4 t^2 + U^2)`. A cosine dispersion is meaningful only on the clean-system
momentum map, where the dashed curve is `epsilon_k = -2 t cos(k)`.

`spectral_plot.vector_map` merges only adjacent energy cells with exactly equal
RGBA under the selected colormap. It verifies that expanding the color runs
reproduces the original map. The underlying spectral arrays are neither averaged
nor interpolated. SVG and PDF remain vectors; PNG is also provided.

## Reproduction

Use the pinned scientific environment in `requirements.txt`. `scripts/python.sh`
accepts `HOLSTEIN_PYTHON=/path/to/python` and sets numerical library thread counts.
From the source root:

```bash
make test
./scripts/python.sh -m pytest -q tests
sbatch --account=che190065 scripts/defect_dense_sample.slurm
# After acquisition and pilot verification:
sbatch --account=che190065 scripts/defect_dense_continue.slurm
# After all continuation fits have been frozen:
sbatch --account=che190065 scripts/defect_dense_reference.slurm
sbatch --account=che190065 scripts/defect_dense_resolution.slurm
./scripts/python.sh scripts/report_defect_dense.py
```

Scripts refuse to overwrite output directories. Choose fresh output names for a
new run. The acquisition checks the verified optimization benchmark and existing
ground-state validation. On another platform, compile and run the benchmark
there before acquiring new data; executable hashes are platform-specific.
