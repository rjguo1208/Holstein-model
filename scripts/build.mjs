import { readFile, writeFile, mkdir, cp, copyFile } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';
import path from 'node:path';
import katex from 'katex';

const root = fileURLToPath(new URL('../', import.meta.url));
const source = await readFile(path.join(root, 'src/index.html'), 'utf8');
let count = 0;
const html = source.replace(/\\\[([\s\S]*?)\\\]|\\\(([\s\S]*?)\\\)/g, (_, display, inline) => {
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
const vendor = path.join(root, 'site/assets/katex');
await mkdir(vendor, { recursive: true });
await copyFile(path.join(root, 'node_modules/katex/dist/katex.min.css'), path.join(vendor, 'katex.min.css'));
await copyFile(path.join(root, 'node_modules/katex/LICENSE'), path.join(vendor, 'LICENSE'));
await cp(path.join(root, 'node_modules/katex/dist/fonts'), path.join(vendor, 'fonts'), { recursive: true });
await copyFile(path.join(root, 'src/style.css'), path.join(root, 'site/assets/style.css'));
await writeFile(path.join(root, 'site/index.html'), `<!-- Generated from src/index.html. Run npm run build. -->\n${html}`);
console.log(`Rendered ${count} LaTeX expressions. HTML, MathML, CSS and fonts are local; no client JavaScript is required.`);
