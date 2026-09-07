"""Export completed calculations for the website; never launch new training.

Usage: python scripts/export_holstein.py /path/to/holstein_lrvmc
Requires NumPy. The website and its checks have no package dependencies.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / "site"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    args = parser.parse_args()
    source = args.source.resolve()
    results = source / "results"
    output = SITE / "data"
    output.mkdir(parents=True, exist_ok=True)
    hashes = {}

    def track(path):
        hashes[str(path.relative_to(source))] = hashlib.sha256(path.read_bytes()).hexdigest()
        return path

    def read(path):
        return json.loads(track(path).read_text())

    def write(path, obj):
        path.write_text(json.dumps(obj, ensure_ascii=False, allow_nan=False,
                                   separators=(",", ":")) + "\n")

    def poles(path):
        with np.load(track(path), allow_pickle=False) as data:
            energy, weight = data["energies"], data["weights"]
            assert len(energy) == len(weight)
            assert np.isfinite(energy).all() and np.isfinite(weight).all()
            assert (weight >= 0).all()
            return {"energies": energy.tolist(), "weights": weight.tolist()}

    rows = read(results / "comparison_vmc_centered/summary.json")
    resampling = read(results / "audit/resampling.json")
    convergence = read(results / "audit/ved_changes.json")
    fullsum = read(results / "comparison_fullsum/summary.json")
    names = {"A": "弱耦合", "B": "中等耦合", "C": "较强耦合", "D": "强耦合"}
    table = []
    for case in "ABCD":
        base = next(r for r in rows if r["case"] == case)
        frequencies = np.load(track(results / "comparison_vmc_centered" / base["omega_file"]),
                              allow_pickle=False)
        data = {"case": case, "label": names[case], "model": base["model"],
                "omega": frequencies.tolist(), "momenta": []}
        for k in [0.0, 0.5, 1.0]:
            display = next(r for r in rows if r["case"] == case and r["k_pi"] == k
                           and r["metric_rcond"] == .001 and r["eta"] == .05)
            ref = read(results / "ved_high" / f"{case}_k{k:.6f}_nh20.json")
            point = {"k_pi": k, "reference": {"E": ref["energy"], "Z": ref["Z"],
                     "generations": ref["generations"], "dimension": ref["dimension"],
                     "residual": ref["residual"],
                     "poles": poles(results / "comparison_vmc_centered" / display["reference_file"])},
                     "vmc": display["measurement"], "lr": [],
                     "convergence": [r for r in convergence if r["case"] == case
                                     and r["k_pi"] == k and r["nh_to"] == 20]}
            for row in resampling:
                if row["case"] != case or row["k_pi"] != k or row["eta"] != .05:
                    continue
                # Use exactly the finalized run selected by the numerical report.
                if row["source_run"] != display["source_run"]:
                    continue
                metrics = [r for r in resampling if r["source_run"] == row["source_run"]
                           and r["draws"] == row["draws"]
                           and r["lr"]["metric_rcond"] == row["lr"]["metric_rcond"]]
                point["lr"].append({"samples": row["lr"]["sample_count"],
                    "rcond": row["lr"]["metric_rcond"], "E": row["E_candidate"],
                    "Z": row["Z_candidate"], "total_weight": row["weight_candidate"],
                    "diagnostics": row["lr"],
                    "poles": poles(source / row["source_run"] / row["candidate_file"]),
                    "metrics": [{key: r[key] for key in ["eta", "delta_E", "delta_Z",
                        "spectral_l1", "rest_spectral_l1", "weight_candidate"]} for r in metrics]})
            assert len(point["lr"]) == 6, (case, k, len(point["lr"]))
            point["lr"].sort(key=lambda r: (r["samples"], -r["rcond"]))
            data["momenta"].append(point)
            table.append({key: display[key] for key in ["case", "k_pi", "E_reference",
                "E_candidate", "Z_reference", "Z_candidate", "delta_E", "delta_Z",
                "spectral_l1", "rest_spectral_l1", "weight_candidate"]})
        write(output / f"{case}.json", data)

    write(output / "calibration.json", [{key: r[key] for key in
          ["case", "k_pi", "hidden", "seed", "delta_E", "delta_Z", "spectral_l1"]}
          for r in fullsum if r["eta"] == .05 and r["lr"]["metric_rcond"] == 1e-10])
    write(output / "ved-convergence.json", convergence)
    write(output / "validation.json", read(results / "audit/artifact_checks.json"))
    downloads = SITE / "downloads"
    downloads.mkdir(exist_ok=True)
    with (downloads / "benchmark.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(table[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(table)

    # Preserve report-relative links and figures; no unrelated project files or logs.
    visited = set()

    def copy_report(relative):
        relative = Path(relative)
        if relative in visited:
            return
        visited.add(relative)
        src = (source / relative).resolve()
        if not src.is_relative_to(source) or not src.is_file():
            raise ValueError(f"Invalid report dependency: {relative}")
        dest = SITE / "reports" / relative
        dest.parent.mkdir(parents=True, exist_ok=True)
        if src.suffix == ".csv":
            dest.write_text(track(src).read_text())
            return
        if src.suffix != ".md":
            shutil.copy2(track(src), dest)
            return
        text = track(src).read_text()

        def link(match):
            prefix, label, target = match.groups()
            if target.startswith(("http:", "https:", "#", "mailto:")):
                return match.group(0)
            linked = (relative.parent / target).as_posix()
            if linked == "README.md":
                return f"{prefix}[{label}](../index.html#methods)"
            if (source / linked).is_file():
                copy_report(linked)
                return match.group(0)
            return label + "（见原始计算项目）"

        text = re.sub(r"(!?)\[([^\]]+)\]\(([^)]+)\)", link, text)
        dest.write_text(text)

    copy_report("RESULTS.md")
    for item in ["results/audit/capacity.pdf", "results/audit/ved_changes.json",
                 "results/audit/coarse_broadening.json"]:
        copy_report(item)
    write(output / "provenance.json", {
        "schema_version": 1, "dataset_date": "2026-09-07",
        "observable": "zero-temperature single-electron addition from vacuum",
        "energy_zero": "vacuum", "t": 1, "coupling": "lambda = g^2 / (2 t omega0)",
        "reference": "infinite-lattice VED, Nh=20, 1600 Lanczos steps",
        "candidate": "L=32, h=40, SR, training seed=0, LR sampling seed=2026",
        "sample_note": "The two sample counts are prefixes of one measurement record.",
        "physical_convergence_certified": False,
        "source_files_sha256": dict(sorted(hashes.items())),
        "exporter_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    print(f"Exported 4 cases, 12 momenta, 72 LR spectra to {output}")


if __name__ == "__main__":
    main()
