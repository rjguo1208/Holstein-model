"""Check published links, rendered math, fonts and LaTeX figure provenance."""
from hashlib import sha256
from html.parser import HTMLParser
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit
import xml.etree.ElementTree as ET

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

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        assert tag not in {"script", "merror"}, f"Unexpected element: {tag}"
        classes = attrs.get("class", "").split()
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


page = Page()
page.feed((SITE / "index.html").read_text())
for link in page.links:
    check_link(SITE / "index.html", link)
for css in (SITE / "assets").rglob("*.css"):
    for link in re.findall(r"url\(['\"]?([^)'\"]+)['\"]?\)", css.read_text()):
        check_link(css, link)

source = (ROOT / "src/index.html").read_text()
expressions = re.findall(r"\\\[([\s\S]*?)\\\]|\\\(([\s\S]*?)\\\)", source)
assert page.math_count == len(expressions) == page.mathml_count, "Incomplete LaTeX/MathML rendering"
assert page.math_count > 0
figures = list((SITE / "figures").glob("*.tex"))
assert page.images == len(figures) == 4
for source in figures:
    svg = source.with_suffix(".svg")
    digest = sha256(source.read_bytes()).hexdigest()
    assert f"source-sha256: {digest}" in svg.read_text(), f"Recompile changed LaTeX: {source.name}"
    tree = ET.parse(svg)
    assert tree.getroot().get("viewBox"), f"Missing scalable bounds: {svg.name}"
    assert not tree.findall(".//{http://www.w3.org/2000/svg}image"), f"Raster image in {svg.name}"
    assert not tree.findall(".//{http://www.w3.org/2000/svg}script"), f"Script in {svg.name}"
    assert source.with_suffix(".pdf").is_file(), f"Missing PDF: {source.name}"

print(f"OK: {page.math_count} LaTeX expressions with MathML; {len(figures)} TikZ vector figures; local links, anchors and fonts.")
