# Completed runs and publication checks — 2026-09-18

This update replaces pending-job messages with completed data, while retaining
the independent VED reference and all earlier source/data archives.

## Confirmed Slurm completion

All jobs completed with exit code 0:0:

| Job | Task | Partition/node | Runtime |
|---|---|---|---|
| 20826452 | 17-momentum map sampling | highmem / b000 | 3m02s |
| 20826455 | Map continuation | highmem / b000 | 2m29s |
| 20825973 | Additional k=0 spectral data | debug / a004 | 56s |
| 20825985 | Additional k=0 long-time data | debug / a004 | 37s |
| 20826180 | Merged k=0 analysis | shared / a240 | 13m32s |

The map has 136/136 completed chains and 3,264,000,000 production steps. Raw
block counts, unique seeds, completion flags, input and source hashes were
checked. No order-cap rejection occurred. Map arrays are finite and nonnegative;
all bootstrap interval arrays are finite.

## Spectral accuracy remains limited

Comparison uses only the nine exactly shared k values, with no interpolation
of the reference, after the actual-data continuation and model selection.
The metric is relative L1 error of the **full** broadened spectrum over
[-3.25,5.5], not the incoherent component or a pointwise error.

| lambda | eta=.25t error range | eta=t error range | largest held-out score |
|---|---|---|---|
| .25 | 5.0–37.2% | 0.44–2.57% | 45.1 |
| .5 | 11.8–43.6% | 1.37–5.11% | 31.7 |

These checks do not validate the other eight k values or fine sideband peaks.
Large held-out scores at some momenta require further statistics/covariance
checks. Outputs retain spectral_resolution_validated=false.

The merged k=0 mass estimates are 1.152053 ± .000396 and 1.350906 ± .000743.
The mass checks pass, with a separate window-dependence assessment; the first
fit still has chi²/dof≈2.25. The selected k=0 continuation methods have
incoherent-spectrum errors 24.1% and 24.5%, so the 10% fine-spectrum goal fails
and the stronger-coupling expansion remains unstarted.

## Reproducible downloads

The main map archive contains conditional blocks, source snapshots, final maps
and the audit. Raw chains are split by coupling to keep individual files below
the repository upload limit. The completed k=0 archive contains the new raw
chains, analysis, references and reports; merged reproduction also uses the
retained original archive for older chains.

Exact sizes, file counts and SHA-256 values are in
site/data/completed-archives.json. All four archives passed ZIP integrity checks:

- diagmc-map-completed-01.zip — 22,062,694 bytes.
- diagmc-map-raw-lambda0.25.zip — 37,542,239 bytes.
- diagmc-map-raw-lambda0.50.zip — 42,955,690 bytes.
- diagmc-completed-03.zip — 48,724,635 bytes.

Actual map NPZ axes are (eta,k,energy), shape (3,17,701), with eta=[.25,.5,1].
The prior reference NPZ remains distinct, shape (4,41,701).

## Rendering

- Build and site checks pass: 141 LaTeX/MathML expressions, four TikZ diagrams,
  and 16 scientific SVG/PDF plots.
- Updated map and results pages passed Chromium checks at 1280, 390 and 320px,
  with JavaScript disabled and external resources blocked.
- No missing images/fonts, math errors or body overflow; long figures and
  formulas scroll within their containers on narrow displays.
- Mobile status text and actual map figures were visually checked. Vector-cell
  seam fixes from the previous publication are retained.
- Current scientific results are labeled preliminary where spectral precision
  has not passed; the VED reference is not relabeled as a DiagMC result.
