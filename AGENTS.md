# Website maintenance

The user requires Chinese and English versions to be updated together.

- Keep the Chinese routes at the site root and their English equivalents under
  `en/`. Every website page needs both versions and the top-right language switch.
- Edit shared content and markup in `src/*.html`, and update all affected English
  translations in `src/locales/en.json` in the same change. Translate headings,
  body text, captions, tables, download descriptions, metadata, alternative text,
  inline SVG labels and explanatory labels inside equations. Preserve the full
  scientific meaning, numerical qualifications and historical status in both.
- Translation keys are the Chinese source strings with whitespace normalized.
  `python3 scripts/localize.py --extract` lists the current keys. Missing, empty,
  obsolete or untranslated entries fail the build. Do not bypass these checks or
  add Chinese fallback text to make an English page build.
- Both languages share equations, anchors, scientific plots and downloadable
  data. Keep scientific plot labels in English or mathematical notation so the
  same assets are readable from both language versions. Translate their page
  captions and alternative text through the catalog.
- Run `npm run build` and `npm run check` and commit the regenerated `site/`
  pages/assets with source changes. Check language switching and navigation on
  desktop and mobile, including the GitHub Pages `/Holstein-model/` prefix.
- Publish both language versions together through the existing Pages workflow.
  Do not mark calculations or narrow spectral features as validated merely
  because page generation, translation or deployment succeeded.
