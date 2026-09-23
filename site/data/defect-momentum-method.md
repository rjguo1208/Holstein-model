# Momentum spectra of a single Holstein defect

The model remains a one-dimensional infinite chain at zero temperature with
one added electron, t = omega0 = U = 1 and lambda = g^2/(2 t omega0) = 0.25, 0.5.
Only the electronic on-site potential at site zero is changed. A single defect
breaks translation symmetry; k labels the incident probe, not an exact conserved
quantum number of the defect Hamiltonian.

## Observable and resolution

Define the normalized coherent injection

    |v(k,L)> = (1/sqrt(L)) sum_{j=-W}^W exp(i k j) c_j^dagger |vac>, L=2W+1.
    A_L(k,w) = -(1/pi) Im <v(k,L)|(w+i eta-H)^-1|v(k,L)>.

The main observable uses W=16 (33 sites centered on the defect). A 65-site
window is also evaluated by VED to show the dependence on the observation
region. The window restricts the injected and removed electron only. It is not
a hard-wall chain, a periodic defect array or a truncation of intermediate
Monte Carlo paths. Enlarging the window changes the observable and dilutes
the relative weight of an isolated impurity; it is distinct from increasing
the numerical boundary R used in VED.

All 201 k points from 0 to pi (delta k=pi/200) are evaluated by an explicit
Fourier contraction of off-diagonal propagation. No momentum interpolation or
Fourier transform of local spectral intensities is used. Nearby k values share
the same real-space information and, for Monte Carlo, correlated noise. The
intrinsic momentum width of the rectangular probe is of order 2 pi/L; adding
plot points does not overcome that finite-window width. The 1601 output energies
span -4.5 to 5.5 with spacing 0.00625. Broadening eta is a Lorentzian half-width,
not an error bar or a certified inference resolution.

The dashed curve is the clean bare-host reference epsilon_k=-2t cos(k).
Absolute spectral weights are retained, with a common vacuum energy zero and
no normalization of individual energy columns. The unbroadened measure has
M0=1. For L>=3 its exact moments are

    M1 = -2t(1-1/L) cos(k) - U/L,
    M2 = 2t^2 + 2t^2(1-2/L) cos(2k) + 4tU cos(k)/L + U^2/L + g^2.

These are moments before Lorentzian broadening. In particular the second moment
of the infinitely extended Lorentzian-broadened curve does not converge.

## VED first

The existing electron-position times relative-phonon-cloud basis is reused.
The Hamiltonian is assembled into CSR and checked against the independently
written matrix-free operator. Processes cache a basis/operator while advancing
through source sites. BLAS/OpenMP use one thread per process.

For each injection site j>=0, an 800-step Lanczos recursion keeps only the
components of each vector on the zero-phonon target sites. These small projected
vectors and the tridiagonal matrix give all nonlocal resolvent columns. Reflection
supplies j<0. Their full double Fourier sum gives every requested k; the expensive
Hamiltonian multiplication is shared across the entire grid. This reduces the
number of large recursions compared with separate complex recursions at all k.

The main basis has cloud generation Nh=12 and electron boundary R=96. Independent
checks compare Nh=10/12 at R=64, R=64/96 at Nh=12, and 400/800 recursion steps.
Both U=0 and U=1, both couplings and every k are included. Five additional directly
injected momenta per parameter pair check the off-diagonal reconstruction using
separate cosine and sine recursions. Exact small-matrix diagonalization and the
infinite-chain g=0 impurity resolvent independently test the method.

At eta=0.25 the largest full-spectrum relative L1 change from Nh=10 to 12 is
0.846% over clean/defect cases. At eta=0.15 it is 3.021%; at eta=0.1 it is 6.749%.
Thus the narrower maps expose cutoff sensitivity and are not described as fully
converged fine structure. Numerical radius and recursion errors are much smaller;
full per-parameter diagnostics and positive-spectrum/moment checks are archived.
Cutoff differences are not rigorous bounds on error in the infinite phonon space.

## Fresh nonlocal bare DiagMC, after VED completion

The real-space continuous-time expansion now permits different electron
endpoints. Single-hop insertion/removal changes the displacement and includes
odd hopping orders. Balanced opposite-hop pairs, phonon-line insertion/removal,
vertex-time moves and external-time rescaling remain available. An exact heat
bath samples the initial site among all translations for which both endpoints
lie in the observation window. All intermediate positions are unrestricted.

For a path with H hopping vertices and N phonon lines, the nonnegative weight
is proportional to t^H g^(2N) exp[-omega0 sum(line durations)+U residence_at_0+mu tau].
The single-hop insertion ratio includes 2t tau/(H+1). An opposite-hop pair includes
t^2 tau^2/[(N_plus+1)(N_minus+1)]. These proposal factors, their inverse moves,
phonon locality, origin heat bath and excursion beyond the observation window
are tested. Conditional integration over external time uses the exact truncated
Gamma law with shape H+2N+1.

The acquisition retains the joint histogram of |j-i| and time, not only local
returns. A cosine factor cos[k(j-i)] projects each histogram into momentum
space. Those projected observations may be noisy and negative; no bins are
clipped. The zero-vertex normalization, divided by L, is known analytically:

    Z_ref/L = [(L-1) F(-mu)+F(-mu-U)]/L,
    F(r)=(1-exp(-r T))/r.

The sampler writes compact float64 blocks rather than large decimal CSV files,
reuses working vectors and parallelizes independent chains. Free hopping, the
g=0 attractive impurity and the atomic electron-phonon limit are checked with
24 independent chains before production. Their largest standardized discrepancy
over the predeclared momenta/time bins is 2.211, below the predeclared threshold 8.

Production uses eight independent chains per coupling, 128 million steps per
chain, one million warm-up steps, thinning 50, 640 blocks of 200,000 steps,
192 time bins and T=12. Total production is 2.048 billion steps. No arc/hop cap
is reached. The auxiliary shift is the analytic spectral lower bound minus 0.1:

    mu = -sqrt(4t^2+U^2) - g^2/omega0 - 0.1t.

It does not use a VED pole energy or spectral shape.

## Continuation and comparison

Positive spectra are reconstructed using exact coherent-window moments and the
analytic lower bound -sqrt(4t^2+U^2)-g^2/omega0. Variance-scaled full covariance
whitening remains defined for signed noisy measurements. The main MaxEnt support
has 641 points up to 10t, dt=0.125, tau<=6 and 1.6-million-step statistical blocks.
Variants use dt=0.25, tau<=8, 321 energies, 3.2-million-step blocks, or covariance
eigenvalue cutoff 1e-6 rather than 1e-8.

Four folds each hold out two complete independent chains. Every configuration
uses the same validation bins dt=0.25, tau<=6. Within a configuration choose
the largest alpha within one fold-standard-error of the minimum mean score,
then choose the configuration with the smallest selected score. All selections
are saved before any VED spectral comparison. This heuristic does not establish
the true resolution of the inverse Laplace transform.

Each k has 32 block-bootstrap repeats, resampling within each chain and reselecting
alpha with the representation fixed. Their 16-84% pointwise empirical ranges are
conditional on the representation; they exclude systematic representation bias
and do not describe independent errors between adjacent k. Spectral pole weights
are stored so broadened bootstrap curves can be reproduced without storing many
redundant output grids. Negative raw Green-function bins, multiple blocking
estimates, candidate sensitivity and comparisons of uncontinued Green functions
are retained. Narrow features are not declared validated by the dense grid.

In the completed run, the largest uncontinued-Green discrepancy from VED for
tau<=6 is 2.471 conservative standard errors, using the maximum error from
several block lengths and whole-chain jackknife estimates. At eta=0.25 the
full-spectrum relative L1 differences span 3.97-12.13% for lambda=0.25 and
6.23-15.67% for lambda=0.5. At eta=1 the corresponding maxima are 1.685% and
1.749%. All 402 selected fits converged, with maximum exact-moment residual
9.295e-8. The retained covariance ranks range from 10 to 16. These checks support
the acquisition and broad spectra; they do not certify narrow peak resolution.

## Reproduction

Use requirements.txt and scripts/python.sh (or set HOLSTEIN_PYTHON explicitly).
Run make all, make test and scripts/python.sh -m pytest -q tests. Production and
postprocessing run on highmem, account che190065, with process-level parallelism.
All acquisition scripts refuse existing output directories. VED sampling,
analytic validation, production Monte Carlo, continuation and reporting are
separate steps with recorded Slurm dependencies and source/input hashes. See
the downloadable job manifest for exact commands and output directories.

The VED report_02 only repairs a panel-title math delimiter in report_01; its
numerical NPZ arrays are identical. Final plots retain the raw color cells, with
only exactly identical neighboring colors merged in vector exports. Web pages
embed PNG previews and link to full SVG/PDF figures and NPZ data in the release.

The distributed archives share a holstein-diagmc/ root. Extract every ZIP into
the same parent directory; the top-level manifest lists the hash of every file.
The direct-download figures and main map NPZ files are also included in the
archives. Numerical reproduction from these archived inputs needs no new Monte
Carlo sampling. For example, from the extracted source directory:

```bash
make all test
bash scripts/python.sh -m pytest -q tests
bash scripts/python.sh scripts/report_defect_momentum.py \
  --data results/defect_momentum_ved_01 --out results/reproduced_ved --jobs 8
bash scripts/python.sh scripts/report_defect_momentum_mc.py \
  --fit results/defect_momentum_fit_01 --out results/reproduced_mc \
  --reference results/defect_momentum_ved_report_02 \
  --ved results/defect_momentum_ved_01
```

For a new acquisition, use the included Slurm scripts with fresh output names.
The default wrappers refer to this cluster's scientific Python environment;
set HOLSTEIN_PYTHON to a compatible environment on another machine. Regenerating
the HTML additionally requires a checkout of the bilingual website repository.
