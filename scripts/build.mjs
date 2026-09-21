import { readFile, writeFile, mkdir, cp, copyFile, readdir } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';
import path from 'node:path';
import { execFileSync } from 'node:child_process';
import katex from 'katex';

const root = fileURLToPath(new URL('../', import.meta.url));
// Localize before writing any pages: a missing/stale translation blocks the
// entire bilingual build instead of silently publishing one language.
const english = JSON.parse(execFileSync('python3', [path.join(root, 'scripts/localize.py')], { encoding: 'utf8' }));
await mkdir(path.join(root, 'site/en'), { recursive: true });
for (const name of (await readdir(path.join(root, 'src'))).filter(name => name.endsWith('.html')).sort()) {
  const source = await readFile(path.join(root, 'src', name), 'utf8');
  for (const lang of ['zh-CN', 'en']) {
    const translated = lang === 'en' ? english[name] : source;
    const prefix = lang === 'en' ? '../' : '';
    const zh = `${prefix}${name}`;
    const en = lang === 'en' ? name : `en/${name}`;
    const alternates = `  <link rel="alternate" hreflang="zh-CN" href="${zh}">\n  <link rel="alternate" hreflang="en" href="${en}">\n`;
    const switcher = `<nav class="language-switch" aria-label="${lang === 'en' ? 'Language' : '语言'}">\n` +
      `  <a href="${zh}" lang="zh-CN" hreflang="zh-CN"${lang === 'zh-CN' ? ' aria-current="true"' : ''}>中文</a>\n` +
      `  <a href="${en}" lang="en" hreflang="en"${lang === 'en' ? ' aria-current="true"' : ''}>English</a>\n</nav>`;
    const bilingual = translated.replace('</head>', `${alternates}</head>`)
      .replace(/(<a class="skip-link"[^>]*>[\s\S]*?<\/a>)/, `$1\n${switcher}`);
    let count = 0;
    const html = bilingual.replace(/\\\[([\s\S]*?)\\\]|\\\(([\s\S]*?)\\\)/g, (_, display, inline) => {
      count++;
      const rendered = katex.renderToString((display ?? inline).trim(), {
        displayMode: display !== undefined,
        output: 'htmlAndMathml',
        throwOnError: true,
        strict: 'error',
        trust: false,
      });
      return display === undefined ? rendered : `<div class="equation">${rendered}</div>`;
    });
    const output = lang === 'en' ? `en/${name}` : name;
    await writeFile(path.join(root, 'site', output), `<!-- Generated from src/${name}${lang === 'en' ? ' and src/locales/en.json' : ''}. Run npm run build. -->\n${html}`);
    console.log(`${output}: rendered ${count} LaTeX expressions with HTML and MathML.`);
  }
}
const vendor = path.join(root, 'site/assets/katex');
await mkdir(vendor, { recursive: true });
await copyFile(path.join(root, 'node_modules/katex/dist/katex.min.css'), path.join(vendor, 'katex.min.css'));
await copyFile(path.join(root, 'node_modules/katex/LICENSE'), path.join(vendor, 'LICENSE'));
await cp(path.join(root, 'node_modules/katex/dist/fonts'), path.join(vendor, 'fonts'), { recursive: true });
await copyFile(path.join(root, 'src/style.css'), path.join(root, 'site/assets/style.css'));
console.log('CSS and fonts are local; no client JavaScript is required.');
