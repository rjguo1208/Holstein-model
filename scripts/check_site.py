"""Validate site links and the exported scientific data using the standard library."""
from html.parser import HTMLParser
import json
import math
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / "site"


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []
        self.ids = set()

    def handle_starttag(self, tag, attrs):
        for name, value in attrs:
            if name in {"href", "src"}:
                self.links.append(value)
            if name == "id":
                assert value not in self.ids, f"Duplicate HTML id: {value}"
                self.ids.add(value)


def check_link(path, link):
    url = urlsplit(link)
    if url.scheme or url.netloc:
        return
    target = (path.parent / unquote(url.path)).resolve() if url.path else path
    assert SITE in target.parents, f"Link escapes site: {path} -> {link}"
    assert target.is_file(), f"Missing link target: {path} -> {link}"
    if url.fragment and target.suffix == ".html":
        parser = Links()
        parser.feed(target.read_text())
        assert unquote(url.fragment) in parser.ids, f"Missing anchor: {link}"


for path in SITE.rglob("*.html"):
    parser = Links()
    parser.feed(path.read_text())
    for link in parser.links:
        check_link(path, link)
for path in SITE.rglob("*.md"):
    for link in re.findall(r"!?\[[^\]]+\]\(([^)]+)\)", path.read_text()):
        check_link(path, link)

count = 0
max_moment_error = 0.
for case in "ABCD":
    data = json.loads((SITE / "data" / f"{case}.json").read_text())
    assert data["case"] == case
    assert [p["k_pi"] for p in data["momenta"]] == [0, .5, 1]
    omega = data["omega"]
    assert len(omega) >= 1701
    assert all(math.isfinite(w) for w in omega)
    assert all(b > a for a, b in zip(omega[:-1], omega[1:]))
    for point in data["momenta"]:
        assert point["reference"]["generations"] == 20
        assert point["reference"]["dimension"] == 2975103
        ref = point["reference"]["poles"]
        source_energy = -2 * math.cos(math.pi * point["k_pi"])
        exact = [1., source_energy, source_energy ** 2 + data["model"]["g"] ** 2]
        for order in range(3):
            moment = math.fsum(w * e ** order for w, e in zip(ref["weights"], ref["energies"]))
            max_moment_error = max(max_moment_error, abs(moment - exact[order]))
            assert abs(moment - exact[order]) < 1e-10
        assert {(r["samples"], r["rcond"]) for r in point["lr"]} == {
            (n, c) for n in (65536, 262144) for c in (.01, .001, .0001)}
        for record in [point["reference"], *point["lr"]]:
            e, w = record["poles"]["energies"], record["poles"]["weights"]
            assert len(e) == len(w) and len(e) > 0
            assert all(math.isfinite(x) for x in e + w)
            assert all(x >= 0 for x in w)
            assert all(b >= a for a, b in zip(e[:-1], e[1:]))
            if "total_weight" in record:
                assert abs(math.fsum(w) - record["total_weight"]) < 1e-11
        assert len(point["convergence"]) == 4
        for r in point["convergence"]:
            if r["eta"] in (.025, .05):
                assert r["reference_certified"] is False
        count += len(point["lr"])
assert count == 72
assert len(json.loads((SITE / "data/calibration.json").read_text())) == 96
assert json.loads((SITE / "data/provenance.json").read_text())["physical_convergence_certified"] is False
print(f"Links and exported data OK: 4 cases, 12 momenta, {count} LR spectra; "
      f"largest VED moment error = {max_moment_error:.3e}")
