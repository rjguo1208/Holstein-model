"""Check published links, rendered math, fonts and LaTeX figure provenance."""
from hashlib import sha256
from html.parser import HTMLParser
import json
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit
import xml.etree.ElementTree as ET
from localize import build as localize, HAN

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / "site"


class Page(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []
        self.ids = set()
        self.math_count = 0
        self.mathml_count = 0
        self.images = 0
        self.lang = None
        self.nav = None
        self.language_links = []
        self.page_links = []
        self.alternates = []
        self.untranslated = []
        self.cell = None
        self.numeric_cells = []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        assert tag not in {"script", "merror"}, f"Unexpected element: {tag}"
        classes = attrs.get("class", "").split()
        if tag == "html":
            self.lang = attrs.get("lang")
        if tag == "nav":
            self.nav = classes
        if tag == "a" and self.nav:
            if "language-switch" in self.nav:
                self.language_links.append(attrs)
            elif "site-nav" in self.nav:
                self.page_links.append(attrs)
        if tag == "link" and attrs.get("rel") == "alternate":
            self.alternates.append(attrs)
        if tag == "td":
            self.cell = []
        if self.lang == "en":
            for name in ("alt", "title", "aria-label", "content"):
                if HAN.search(attrs.get(name, "")):
                    self.untranslated.append(attrs[name])
        assert "katex-error" not in classes, "Unrendered formula"
        self.math_count += "katex" in classes
        self.mathml_count += tag == "math"
        if tag == "img":
            assert attrs.get("alt"), "A diagram needs alternative text"
            self.images += 1
        if "id" in attrs:
            assert attrs["id"] not in self.ids, f"Duplicate id: {attrs['id']}"
            self.ids.add(attrs["id"])
        for name in ("src", "href"):
            if name in attrs:
                self.links.append(attrs[name])

    def handle_endtag(self, tag):
        if tag == "nav":
            self.nav = None
        if tag == "td" and self.cell is not None:
            self.numeric_cells.append(re.findall(r"\d+(?:\.\d+)?", "".join(self.cell)))
            self.cell = None

    def handle_data(self, data):
        if self.lang == "en" and HAN.search(data):
            if not (self.nav and "language-switch" in self.nav and data.strip() == "中文"):
                self.untranslated.append(data)
        if self.cell is not None:
            self.cell.append(data)


def check_link(path, link):
    url = urlsplit(link)
    if url.scheme or url.netloc:
        return
    target = (path.parent / unquote(url.path)).resolve() if url.path else path
    assert SITE in target.parents, f"Link escapes site: {link}"
    assert target.is_file(), f"Missing local resource: {link}"
    if url.fragment and target.suffix == ".html":
        parser = Page()
        parser.feed(target.read_text())
        assert unquote(url.fragment) in parser.ids, f"Missing anchor: {link}"


pages = {}
catalog = json.loads((ROOT / "src/locales/en.json").read_text())
english, _ = localize(catalog)
names = sorted(english)
expected_pages = set(names) | {"en/" + name for name in names}
math_count = 0
for path in sorted(SITE.rglob("*.html")):
    relative = path.relative_to(SITE).as_posix()
    assert relative in expected_pages, f"Unexpected or unpaired page: {relative}"
    page = Page()
    page.feed(path.read_text())
    pages[relative] = page
    for link in page.links:
        check_link(path, link)
    is_english = relative.startswith("en/")
    source = english[path.name] if is_english else (ROOT / "src" / path.name).read_text()
    expressions = re.findall(r"\\\[([\s\S]*?)\\\]|\\\(([\s\S]*?)\\\)", source)
    assert page.math_count == len(expressions) == page.mathml_count, "Incomplete LaTeX/MathML rendering"
    assert page.math_count > 0
    assert page.lang == ("en" if is_english else "zh-CN"), f"Incorrect language: {relative}"
    assert not page.untranslated, f"Untranslated English content in {relative}: {page.untranslated}"
    targets = {"zh-CN": ("../" if is_english else "") + path.name,
               "en": ("" if is_english else "en/") + path.name}
    assert len(page.language_links) == 2, f"Missing language switch: {relative}"
    assert {link.get("hreflang"): link.get("href") for link in page.language_links} == targets
    assert [link["hreflang"] for link in page.language_links if link.get("aria-current") == "true"] == [page.lang]
    assert {link.get("hreflang"): link.get("href") for link in page.alternates} == targets
    assert sorted(link["href"] for link in page.page_links) == names, f"Incomplete page navigation: {relative}"
    assert [link["href"] for link in page.page_links if link.get("aria-current") == "page"] == [path.name]
    math_count += page.math_count
assert set(pages) == expected_pages, "Every source page must have Chinese and English outputs"
for name in names:
    zh, en = pages[name], pages["en/" + name]
    assert zh.ids == en.ids, f"Section anchors differ between languages: {name}"
    assert zh.numeric_cells == en.numeric_cells, f"Numerical table data differ between languages: {name}"
    def resources(relative, page):
        return [(SITE / relative).parent.joinpath(urlsplit(link).path).resolve()
                for link in page.links if urlsplit(link).path
                and not urlsplit(link).scheme and not urlsplit(link).netloc
                and not urlsplit(link).path.endswith(".html")]
    assert resources(name, zh) == resources("en/" + name, en), f"Scientific resources differ between languages: {name}"
for css in (SITE / "assets").rglob("*.css"):
    for link in re.findall(r"url\(['\"]?([^)'\"]+)['\"]?\)", css.read_text()):
        check_link(css, link)

figures = list((SITE / "figures").glob("*.tex"))
assert pages['index.html'].images == len(figures) == 4
for source in figures:
    svg = source.with_suffix(".svg")
    digest = sha256(source.read_bytes()).hexdigest()
    assert f"source-sha256: {digest}" in svg.read_text(), f"Recompile changed LaTeX: {source.name}"
    tree = ET.parse(svg)
    assert tree.getroot().get("viewBox"), f"Missing scalable bounds: {svg.name}"
    assert not tree.findall(".//{http://www.w3.org/2000/svg}image"), f"Raster image in {svg.name}"
    assert not tree.findall(".//{http://www.w3.org/2000/svg}script"), f"Script in {svg.name}"
    assert source.with_suffix(".pdf").is_file(), f"Missing PDF: {source.name}"

plots = list((SITE / "results").glob("*.svg"))
for plot in plots:
    tree = ET.parse(plot)
    assert tree.getroot().get("viewBox"), f"Missing scalable bounds: {plot.name}"
    assert not tree.findall(".//{http://www.w3.org/2000/svg}image"), f"Raster image in {plot.name}"
    assert plot.with_suffix('.pdf').is_file(), f"Missing PDF: {plot.name}"
print(f"OK: {len(pages)} bilingual pages; {math_count} LaTeX expressions with MathML; {len(figures)} TikZ figures; {len(plots)} scientific plots; translations, paired language navigation, identical numerical tables/resources, local links, anchors and fonts.")
