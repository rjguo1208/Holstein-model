"""Fresh extraction, checksum audit, numerical regeneration and fixed-seed replay."""
from __future__ import annotations
import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys
import zipfile
import xml.etree.ElementTree as ET
import numpy as np


def digest(path):
    h = sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''):
            h.update(b)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--release', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    mf = json.loads((args.release/'manifest.json').read_text())
    assert mf['complete']
    args.out.mkdir(parents=True, exist_ok=False)
    for item in mf['archives']+mf['assets']:
        file = args.release/item['name']
        assert file.stat().st_size == item['bytes'] and digest(file) == item['sha256']
        if file.suffix == '.svg':
            tree = ET.parse(file)
            assert tree.getroot().get('viewBox')
            assert not tree.findall('.//{http://www.w3.org/2000/svg}image')
            assert not tree.findall('.//{http://www.w3.org/2000/svg}script')
    for item in mf['archives']:
        with zipfile.ZipFile(args.release/item['name']) as z:
            for name in z.namelist():
                assert name.startswith('holstein-diagmc/') and '..' not in Path(name).parts
            assert z.testzip() is None
            z.extractall(args.out)
    for item in mf['files']:
        file = args.out/item['path']
        assert file.stat().st_size == item['bytes'] and digest(file) == item['sha256']
    root = args.out/'holstein-diagmc'
    with (args.out/'tests.log').open('w') as log:
        subprocess.run(['make', 'all', 'test'], cwd=root, stdout=log, stderr=subprocess.STDOUT, check=True)
        subprocess.run([sys.executable, '-m', 'pytest', '-q', 'tests'], cwd=root,
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    sys.path.insert(0, str(root))
    from defect_dense import curves
    recovered = []
    for lam in [.25, .5]:
        for index in [0, 100, 200]:
            path = root/'results/defect_momentum_fit_01'/f'lambda{lam:.2f}_k{index:03d}'
            with np.load(path/'spectra.npz') as f:
                selected = curves(dict(energies=f['energies'], weights=f['weights']), tuple(f['eta']))
                np.testing.assert_allclose(selected, f['selected'], atol=2e-13, rtol=2e-13)
                with np.load(path/'bootstrap.npz') as bs:
                    replicas = np.array([curves(dict(energies=bs['energies'], weights=w), tuple(f['eta'])) for w in bs['weights']])
                np.testing.assert_allclose(np.quantile(replicas, [.16, .5, .84], axis=0), f['bootstrap_quantiles'], atol=2e-13, rtol=2e-13)
            recovered.append(dict(coupling=lam, k_index=index, selected_and_bootstrap_reproduced=True))
    pilot = root/'results/defect_momentum_mc_pilot_01'
    p = next(t for t in json.loads((pilot/'manifest.json').read_text())['tasks'] if t['name'] == 'free_s00')
    command = [str(root/'build/defect_momentum_diagmc')]
    for key, value in p['parameters'].items():
        command += ['--'+key, str(value)]
    replay = args.out/'sampler_replay'
    command += ['--out', str(replay)]
    subprocess.run(command, check=True)
    assert digest(replay/'blocks.f64') == digest(pilot/'chains/free_s00/blocks.f64')
    result = dict(passed=True, extracted_files=len(mf['files']), archives=len(mf['archives']),
                  verification_job_id=os.environ.get('SLURM_JOB_ID'),
                  direct_assets=len(mf['assets']), build_and_tests_passed=True,
                  python_tests=26, regenerated_momenta=recovered,
                  analytic_sampler_replay_bitwise_identical=True,
                  release_manifest_sha256=digest(args.release/'manifest.json'))
    (args.release/'verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
