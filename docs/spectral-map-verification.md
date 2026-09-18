# Spectral map verification — 2026-09-18

## Scope and provenance

The new spectral-map.html displays **VED/Lanczos reference data**, not a completed
DiagMC map. Both the page lead and image title state the method. Existing
DiagMC spectral data cover k=0; the previous small-k dispersion data are not
interpolated or relabeled as full-zone spectra.

The reference covers lambda=.25,.5, t=omega0=1, zero temperature, one electron,
41 independent momenta from 0 to pi and 701 energy values [-3.25,5.5].
Both panels use the same vacuum energy zero and absolute t*A color scale.
No per-column normalization, energy shift, k interpolation or image smoothing
is performed. The second image only changes the color mapping to logarithmic.
The main Lorentzian half width is eta=.25t.

Source code, model routines and hashes, spectral poles, validation records,
full maps and plotting source are in data/diagmc-spectral-map-01.zip:

- Size: 7,824,276 bytes; 240 files.
- SHA-256: 769935d258e8727df254f7b1cef7321cbc37410c9c70a1d5febfdd8947f419be.
- The archived reference routines include a package-layout copy for standalone
  reproduction without the separate NNQS checkout.
- ZIP integrity and source hashes were verified.

## Numerical checks

- Main reference: Nh16, dimension 178617; all 41 momenta also at Nh14.
- Nh18 at k=0, pi/4, pi/2, 3pi/4, pi.
- Maximum whole-spectrum relative L1 changes at eta=.25t:
  14→16: 0.1020% / 0.0956%; 16→18 at five check points: 0.0253% / 0.0366%.
- 400→800 Lanczos steps checked at the five momenta in each basis:
  largest relative L1 change <3.26e-11.
- Unbroadened M0/M1/M2 maximum absolute residual <1.20e-14.
- Arrays have shape (4 etas, 41 momenta, 701 energies), finite positive values,
  correct axis ordering and endpoints.
- A nonzero-k column of each coupling was independently regenerated directly
  from saved poles and compared with the corresponding plotted data.
- Main-window integrated weight ranges .9187–.9624 (lambda=.25) and
  .9068–.9606 (lambda=.5). No rescaling conceals omitted Lorentzian tails.

The comparisons test finite-space and recursion stability at stated broadening,
not a rigorous bound on the exact solution or all narrow spectral features.

## Submitted DiagMC calculation

Sampling job 20826452: 17 momenta × 2 couplings × 4 chains = 136 chains,
24 million production steps per chain, four CPUs. The initial array submission
exceeded a submit-count limit; a single batch processing successive momenta was
accepted. Publication-time state: PENDING (Priority).

Analysis job 20826455 depends on successful sampling: PENDING (Dependency).
It saves full-spectrum nonnegative continuation, held-out-chain selection,
block bootstrap and map figures; it does not publish automatically.
The pipeline does not impose E/Z priors or use a reference spectrum as input.
All output is marked preliminary, with spectral_resolution_validated=false.

The shared scientific Python tests pass (12 total). New tests cover the free
particle limit without a pole prior, absolute spectral weight and array axes,
and integrated-G preservation under time rebinning.

A three-momentum, four-chain, small-statistics end-to-end run exercised
acquisition, conditional summaries, held-out selection, one bootstrap replicate,
serialization and plotting. These temporary smoke-test spectra are not included
as scientific results. An initial smaller test correctly stopped when a bootstrap
replicate had insufficient covariance rank; the complete test used 80 blocks.

## Display checks

- npm run build / npm run check pass: three pages, 136 LaTeX/MathML expressions,
  four genuine TikZ diagrams and 12 scientific SVG/PDF figures.
- New page: 19 LaTeX expressions; all have MathML.
- Chromium at 1280, 390 and 320px, with page JavaScript disabled and external
  requests blocked: no missing images/fonts, errors or body-level overflow.
- Long equations and figures scroll within their containers on narrow screens.
- SVG contains vector cells, labels and colorbars; no embedded raster image.
- PNG previews and standalone vector PDFs are available.
- Desktop plots and mobile introductory layout were visually inspected.

SVG mesh-only crispEdges and same-color PDF cell borders remove viewer hairline artifacts; the revised browser and PDF exports were visually checked. No spectral array was changed.
