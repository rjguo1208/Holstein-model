"""Compile the self-contained LaTeX/TikZ figures to PDF and vector SVG."""
from hashlib import sha256
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "site" / "figures"
tectonic = os.environ.get("TECTONIC", "tectonic")
if not shutil.which(tectonic) or not shutil.which("pdftocairo"):
    raise SystemExit("Install Tectonic and Poppler (pdftocairo), or set TECTONIC to its executable path.")

for source in sorted(FIGURES.glob("*.tex")):
    with tempfile.TemporaryDirectory(prefix="holstein-tikz-") as output:
        subprocess.run([tectonic, "--keep-logs", "--outdir", output, str(source)], check=True)
        pdf = Path(output) / source.with_suffix(".pdf").name
        svg = source.with_suffix(".svg")
        subprocess.run(["pdftocairo", "-svg", str(pdf), str(svg)], check=True)
        data = svg.read_text()
        digest = sha256(source.read_bytes()).hexdigest()
        end = data.index("?>") + 2 if data.startswith("<?xml") else 0
        data = data[:end] + f"\n<!-- source-sha256: {digest} -->\n" + data[end:]
        svg.write_text(data)
        shutil.copyfile(pdf, source.with_suffix(".pdf"))
        print(f"Compiled {source.name} -> PDF + SVG")
