# Chinese and English website — 2026-09-21

All three existing pages have complete English versions under `en/`: theory
notes, numerical results, and spectral maps including the historical records.
The top-right Chinese/English links switch to the corresponding page; page
navigation stays in the selected language. Links work without JavaScript.

## Shared content and future updates

`src/*.html` retains the shared markup and Chinese content. The 464 reviewed
entries in `src/locales/en.json` translate prose, titles, navigation, captions,
table labels, download descriptions, metadata, alternative text, inline SVG
labels and explanatory text inside equations. The English pages use the same
scientific plots and downloadable files. No sampling, inference or numerical
data changed in this update.

`scripts/localize.py` uses normalized source text as translation keys. A changed
Chinese passage requires a corresponding English entry. Missing, empty,
obsolete or untranslated entries fail the build. Mathematical expressions
must agree across languages, allowing translated `\text{...}` explanations and
different ordering of inline expressions within English prose. HTML entities
and LaTeX alignment ampersands retain their original spelling.

The bilingual maintenance requirement is recorded in `AGENTS.md` and README.
Both the CI and Pages workflows run the bilingual build and site checks, and
verify that committed publication files match the generated output.

## Verification

- `npm run build` renders six pages: 201 theory expressions, 41 results
  expressions and 37 spectral-map expressions per language, 558 in total.
- `npm run check` passes: complete language pairs and language metadata,
  corresponding-page switches, same-language page navigation, identical anchors
  and numerical table entries, identical scientific resource targets, valid
  local links/fonts, four TikZ figures and 31 scientific SVG/PDF plot pairs.
- The checks pass with both the local Python 3.6 and Python 3.12 environments.
- Five deliberate invalid catalogs are rejected: a missing entry, an empty
  entry, an obsolete entry, Chinese fallback text and a changed equation.
- Chromium/Playwright passes all 18 combinations of six pages and viewport
  widths 1280, 390 and 320, with JavaScript disabled and external requests
  blocked. Checks cover actual language-switch clicks in both directions,
  page navigation, anchors, completed images, MathML, shared downloads, opening
  historical details and no horizontal overflow of the document body.
- Language links have at least 44 × 44 pixel targets and appear at the top
  right. The site directory stays on the left on desktop and above the content
  on narrow screens. Long tables, equations and figures keep their own scroll
  containers. The language switch is hidden in print; English theory PDF output
  was checked.
- Browser checks use the `/Holstein-model/` URL prefix. Screenshots of English
  desktop/mobile pages and the translated worked-weight example were reviewed.
- `git diff --check` passes. The Chinese scientific source pages and files under
  `site/data/`, `site/results/` and `site/figures/` are unchanged.

These checks concern translation, navigation and publication. They do not
change the documented limits on analytic continuation or narrow-peak accuracy.
