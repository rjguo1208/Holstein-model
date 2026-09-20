"""Extend the frozen spectral refinement protocol to every measured momentum.

All sampling and inference run inside the Slurm allocation. No VED spectrum
is read here. Completed pilot points are reused byte for byte.
"""
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import traceback

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pilot_compare import run as compare
from pilot_bootstrap import run as bootstrap
from pilot_diagnostics import diagnose, public, statistic
from refinement import candidates

PLAN = json.loads((ROOT / "plans/pilot.json").read_text())
CPUS = int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))
JOB = os.environ.get("SLURM_JOB_ID", "local")


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def write(path, value):
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(public(value), indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def status(phase, **extra):
    record = dict(phase=phase, job_id=JOB, utc=datetime.now(timezone.utc).isoformat(), **extra)
    write(ROOT / "results/pipeline_status.json", record)
    print(json.dumps(record), flush=True)


def verify_reuse():
    reused = json.loads((ROOT / "results/reuse_manifest.json").read_text())
    for name, expected in reused["sha256"].items():
        assert digest(ROOT / name) == expected, name
    return reused


def source_snapshot():
    destination = ROOT / "results/inference_source"
    sources = list(ROOT.glob("*.py")) + list((ROOT / "scripts").glob("*"))
    sources += list((ROOT / "src").glob("*")) + [ROOT / "plans/pilot.json", ROOT / "Makefile"]
    hashes = {str(p.relative_to(ROOT)): digest(p) for p in sources if p.is_file()}
    manifest = ROOT / "results/inference_source.json"
    if manifest.exists():
        assert json.loads(manifest.read_text())["sha256"] == hashes
    else:
        for name in hashes:
            target = destination / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT / name, target)
        write(manifest, dict(sha256=hashes, reference_used=False,
                             utc=datetime.now(timezone.utc).isoformat()))


def sample(ik):
    output = ROOT / f"results/pilot_data/ik{ik:03d}"
    existing = list(output.glob("lambda*/run.json"))
    if len(existing) == 16 and (output / "summary.json").exists():
        if all(json.loads(p.read_text()).get("complete") for p in existing):
            return ik
    if output.exists():
        backup = ROOT / "results/interrupted_sampling"
        backup.mkdir(exist_ok=True)
        output.rename(backup / f"ik{ik:03d}-job{JOB}")
    command = [sys.executable, str(ROOT / "scripts/run_suite.py"), "spectral_map",
               "--map-points", "41", "--map-index", str(ik), "--out", str(output),
               "--couplings", ".25", ".5", "--chains", "8", "--blocks", "480",
               "--steps-per-block", "100000", "--seed-offset", str(PLAN["seed_offset"]),
               "--rb-thin", "100", "--jobs", str(min(8, CPUS))]
    with (ROOT / f"logs/sample-ik{ik:03d}-job{JOB}.out").open("w") as log:
        subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
    print(f"Completed new sampling at k-index {ik}", flush=True)
    return ik


def raw_audit():
    rows = []
    hashes = {}
    for ik in PLAN["indices"]:
        folder = ROOT / f"results/pilot_data/ik{ik:03d}"
        manifest = json.loads((folder / "manifest.json").read_text())
        assert manifest["executable_sha256"] == PLAN["baseline_executable_sha256"]
        for name, expected in manifest["source_hashes"].items():
            assert digest(folder / "source" / name) == expected
        for name in ["src/main.cpp", "src/diagmc.hpp", "src/conditional.hpp"]:
            assert digest(ROOT / name) == manifest["source_hashes"][name]
        for lam in PLAN["couplings"]:
            paths = sorted(folder.glob(f"lambda{lam:.2f}_*/run.json"))
            assert len(paths) == 8
            for path in paths:
                meta = json.loads(path.read_text())
                assert meta["complete"] and meta["blocks"] == 480 and meta["steps_per_block"] == 100000
                assert meta["warmup"] == 500000 and meta["rb_thin"] == 100
                assert meta["bins"] == 96 and meta["tau_max"] == 12 and meta["max_order"] == 96
                assert meta["mu"] == -2.8 and meta["t"] == meta["omega"] == 1
                assert np.isclose(meta["g"] ** 2, 2 * lam, rtol=1e-14)
                assert np.isclose(meta["k"], ik * np.pi / 40, rtol=1e-14, atol=1e-14)
                assert sum(meta["attempted"]) == 48000000 and meta["cap_attempts"] == 0
                for name in ["run.json", "blocks.csv", "rb_blocks.csv", "orders.csv"]:
                    raw = path.parent / name
                    hashes[str(raw.relative_to(ROOT))] = digest(raw)
                    if name in ["blocks.csv", "rb_blocks.csv"]:
                        with raw.open() as handle:
                            header = handle.readline().strip().split(",")
                            count = 0
                            for count, line in enumerate(handle, start=1):
                                fields = line.strip().split(",")
                                assert len(fields) == len(header) and int(fields[0]) == count - 1
                            assert count == 480
                rows.append(meta)
    baseline_seeds = {json.loads(p.read_text())["seed"] for p in Path(PLAN["baseline"]).glob("ik*/lambda*/run.json")}
    assert len(baseline_seeds) == 328, "The complete original baseline must be available for seed checks"
    seeds = {r["seed"] for r in rows}
    assert len(rows) == len(seeds) == 656 and not seeds.intersection(baseline_seeds)
    record = dict(chains=656, complete_chains=656, cap_attempts=0,
                  new_chains=496, reused_chains=160, raw_block_counts_verified=True,
                  source_hashes_verified=True, executable_unchanged=True,
                  unique_disjoint_seeds=True, baseline_seed_count=328, production_steps=PLAN["total_refined_production_steps"],
                  new_production_steps=PLAN["new_production_steps"], raw_sha256=hashes)
    write(ROOT / "results/raw_audit.json", record)


def diagnostics():
    out = ROOT / "results/diagnostics/summary.json"
    if out.exists():
        assert len(json.loads(out.read_text())["points"]) == 164
        return
    prior = json.loads((ROOT / "results/reused/diagnostics.json").read_text())
    tasks = [("pilot", str(ROOT / "results/pilot_data"), lam, ik)
             for lam in PLAN["couplings"] for ik in PLAN["new_indices"]]
    with ProcessPoolExecutor(max_workers=min(32, CPUS)) as pool:
        points = prior["points"] + list(pool.map(diagnose, tasks))
    comparisons = []
    for lam in PLAN["couplings"]:
        for ik in PLAN["indices"]:
            a = next(r for r in points if (r["dataset"], r["coupling"], r["index"]) == ("baseline", lam, ik))
            b = next(r for r in points if (r["dataset"], r["coupling"], r["index"]) == ("pilot", lam, ik))
            ta, ga, ca = (np.asarray(a[key]) for key in ["tau", "G", "covariance"])
            tb, gb, cb = (np.asarray(b[key]) for key in ["tau", "G", "covariance"])
            np.testing.assert_array_equal(ta, tb)
            use = ta <= 8
            record = statistic(ga[use], gb[use], ca[np.ix_(use, use)], cb[np.ix_(use, use)])
            record.update(coupling=lam, index=ik, sigma_ratio_pilot_over_baseline=np.sqrt(np.diag(cb) / np.diag(ca)))
            comparisons.append(record)
    out.parent.mkdir(exist_ok=True)
    write(out, dict(complete=True, points=points, comparisons=comparisons))


def has_point(folder, lam, ik):
    stem = ROOT / f"results/{folder}/lambda{lam:.2f}_ik{ik:03d}"
    try:
        meta = json.loads(Path(str(stem) + ".json").read_text())
        with np.load(Path(str(stem) + ".npz")) as data:
            key = "primary" if folder == "comparison" else "samples"
            assert np.isfinite(data[key]).all()
        return True
    except (OSError, ValueError, KeyError, AssertionError):
        return False


def inference():
    tasks = [(lam, ik) for lam in PLAN["couplings"] for ik in PLAN["new_indices"]
             if not has_point("comparison", lam, ik)]
    with ProcessPoolExecutor(max_workers=min(64, CPUS)) as pool:
        for future in as_completed([pool.submit(compare, task) for task in tasks]):
            print("Completed selection: " + json.dumps(public(future.result())), flush=True)
    records = []
    for lam in PLAN["couplings"]:
        for ik in PLAN["indices"]:
            path = ROOT / f"results/comparison/lambda{lam:.2f}_ik{ik:03d}.json"
            r = json.loads(path.read_text())
            assert has_point("comparison", lam, ik)
            records.append(dict(coupling=lam, index=ik, primary=r["primary"],
                                score=r["primary_cv"]["score"], failures=len(r["failures"])))
    write(ROOT / "results/comparison/provenance.json", dict(complete=True, job_id=JOB,
          reference_used_in_inference=False, candidates=candidates(), points=records,
          reused_indices=PLAN["reused_indices"], new_indices=PLAN["new_indices"],
          source_manifest="../inference_source.json", reused_source_manifest="../reused/comparison-provenance.json"))
    verify_reuse()
    hashes = {p.name: digest(p) for p in sorted((ROOT / "results/comparison").glob("*.json"))}
    write(ROOT / "results/evaluation_freeze.json", dict(comparison_hashes=hashes,
          frozen_utc=datetime.now(timezone.utc).isoformat(), all_82_points_frozen=True,
          reference_used_in_selection=False, method_frozen_in_completed_pilot=True))


def uncertainty():
    tasks = [(lam, ik) for lam in PLAN["couplings"] for ik in PLAN["new_indices"]
             if not has_point("bootstrap", lam, ik)]
    with ProcessPoolExecutor(max_workers=min(64, CPUS)) as pool:
        for future in as_completed([pool.submit(bootstrap, task) for task in tasks]):
            r = future.result()
            print(f"Completed bootstrap: lambda={r['coupling']} index={r['index']} draws={r['successful']}", flush=True)
    points = [json.loads((ROOT / f"results/bootstrap/lambda{lam:.2f}_ik{ik:03d}.json").read_text())
              for lam in PLAN["couplings"] for ik in PLAN["indices"]]
    write(ROOT / "results/bootstrap/summary.json", dict(complete=True, points=points, job_id=JOB))
    if any(p["failures"] for p in points):
        status("repairing_original_failed_bootstrap_draws")
        subprocess.run([sys.executable, str(ROOT / "scripts/pilot_repair_bootstrap.py")], cwd=ROOT, check=True)
    final = json.loads((ROOT / "results/bootstrap/summary.json").read_text())
    assert len(final["points"]) == 82
    assert all(p["successful"] == 128 and not p["failures"] for p in final["points"])


def main():
    assert "SLURM_JOB_ID" in os.environ, "Run production work in a compute allocation"
    assert digest(ROOT / "build/diagmc") == PLAN["baseline_executable_sha256"]
    source_snapshot()
    status("verifying_reused_pilot")
    verify_reuse()
    status("sampling_remaining_31_momenta")
    with ThreadPoolExecutor(max_workers=max(1, CPUS // min(8, CPUS))) as pool:
        list(pool.map(sample, PLAN["new_indices"]))
    status("auditing_all_656_chains")
    raw_audit()
    status("covariance_diagnostics")
    diagnostics()
    status("selecting_remaining_62_spectra")
    inference()
    status("bootstrap_remaining_62_spectra")
    uncertainty()
    verify_reuse()
    status("complete", complete=True, points=82, bootstrap_draws=82 * 128,
           new_chains=496, reused_chains=160)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        status("failed", complete=False, error=str(error), traceback=traceback.format_exc())
        raise
