# 41-point DiagMC submission — 2026-09-18

Update (2026-09-19): the 41-point calculation has now completed on Kestrel. See [the completion record](dense-map-kestrel-completed.md). The Anvil queue state below is the historical submission snapshot.

The actual DiagMC grid is now configured to match the existing VED grid exactly:
k_j=j*pi/40 for j=0,...,40. The current published DiagMC figure still contains
17 measured points and remains correctly labeled until the new jobs complete.

Nine old momenta coincide with the target grid. Their 72 raw chains and conditional
blocks were copied without changing data, names, seeds or metadata. The other
32 target momenta require 256 new independent chains at the existing 24 million
production steps per chain. The final map will use 328 chains, totaling
7.872 billion production steps, of which 6.144 billion are newly requested.

The preparation checked the VED grid against both couplings, matched physical
parameters, required a one-to-one reuse map, and checked that all 256 planned
new seeds are unique and disjoint from the entire previous 136-chain map.
The manifest includes hashes of all copied files. The production executable
and numerical continuation algorithm are unchanged.

Sampling job 20829123 and dependent continuation/audit job 20829124 were accepted
by highmem. Their resource requests are 4 CPUs / 8 GB / 30 minutes and
1 CPU / 4 GB / 2 hours. At submission, the sampling job was PENDING with reason
ReqNodeNotAvail, Reserved for maintenance; analysis was PENDING on
afterok:20829123. The full-cluster maintenance reservation runs from
2026-09-18 23:30 to 2026-09-20 21:00 America/Indianapolis (UTC-4). The scheduler
had not assigned a start time.

The dependent job will automatically produce the map and audit all 41 matching
reference momenta. Acquisition completion is separate from spectral accuracy;
finer momentum sampling does not validate narrow frequency features.

The website adds a concise status note, the new grid formula and downloadable
job/grid records plus the reproducible submission sources. Existing 17-point
results, numerical data and reference figures remain intact. Source syntax,
Slurm syntax, the prepared reuse data and the site's math/link checks are
verified for this update; the new 256-chain calculation has not yet run.

The site build passes with 142 rendered LaTeX/MathML expressions. Chromium
checks pass at 1280, 390 and 320px with JavaScript disabled: 25 formulas on the
map page, no missing images/fonts, no math errors and no page-width overflow.
Long formulas and plots remain scrollable on narrow screens. The new mobile
status note and grid formula were visually checked.

The 58-file submission archive passes ZIP integrity checking; its size is
132,627 bytes and SHA-256 is
e1e48d1b3c430dd39dd062b7d50f0091d9fc9658038cd86258038455c16b139c.
