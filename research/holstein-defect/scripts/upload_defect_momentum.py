"""Upload a verified release, audit server hashes, then publish its assets."""
from hashlib import sha256
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent/'Holstein-model'
OUT = ROOT/'results/defect_momentum_release_02'
TAG = 'defect-momentum-20260923'


def run(args):
    return subprocess.run(args, cwd=REPO, text=True, capture_output=True, check=True).stdout


def digest(path):
    h = sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''):
            h.update(b)
    return h.hexdigest()


def main():
    assert json.loads((OUT/'verification.json').read_text())['passed']
    assert not run(['git', 'status', '--porcelain']).strip()
    commit = run(['git', 'rev-parse', 'HEAD']).strip()
    # A fresh fetch before upload is required by the user's instruction.
    run(['git', 'fetch', 'origin'])
    assert commit == run(['git', 'rev-parse', 'origin/main']).strip()
    notes = ROOT/'results/defect_momentum_release_notes.md'
    notes.write_text('''Dense momentum spectra for one attractive on-site defect in the 1D, zero-temperature, one-electron Holstein model.

Parameters: t = omega0 = U = 1; lambda = 0.25 and 0.5. The coherent probe covers 33 sites centered on the defect, with 201 momenta from 0 to pi and 1601 displayed energies. Every figure includes the bare-host cosine guide. VED also provides a clean-system comparison, defect-minus-clean maps and a 65-site probe.

VED was completed first. Fresh bare DiagMC then sampled open electron endpoints using 16 independent chains and 2.048 billion production steps on highmem. The archive includes all raw displacement/time blocks, 12864 conditional bootstrap repeats, frozen reference-free continuation choices, VED recursions, convergence checks and source code.

Narrow spectral peaks are not certified as resolved. A denser k grid does not remove finite-window momentum width, phonon-cloud cutoff error or inverse-Laplace uncertainty. See DEFECT_MOMENTUM.md and the bilingual page for quantitative limitations.

Extract all part ZIPs into the same parent directory. manifest.json lists every archived file and its SHA-256. Individual SVG/PDF/PNG figures and the principal NPZ maps are also supplied for direct download. verification.json records fresh-extraction checks, compiled/Python tests, representative spectral/bootstrap regeneration and bitwise replay of a fixed-seed analytic sampler.

中文：新增单缺陷动量谱使用201个 k 点与1601个能量点，并以33格点相干窗口定义观测量。先完成 VED，再以非局域 bare DiagMC 全新采样20.48亿步；保留原始链、12864次条件 bootstrap、无 VED 谱形先验的选参记录及收敛检查。窄峰仍未通过定量分辨率验证。所有分卷解压到同一父目录，完整哈希清单与全新解压复现检查随包提供。
''')
    result = subprocess.run(['gh', 'release', 'view', TAG, '--json', 'isDraft,targetCommitish'],
                            cwd=REPO, text=True, capture_output=True)
    if result.returncode:
        run(['gh', 'release', 'create', TAG, '--target', commit, '--draft', '--latest=false',
             '--title', 'Single-defect momentum spectra: dense VED and nonlocal DiagMC',
             '--notes-file', str(notes)])
    else:
        existing = json.loads(result.stdout)
        assert existing['isDraft'], 'Do not overwrite an already public immutable release'
    mf = json.loads((OUT/'manifest.json').read_text())
    names = [r['name'] for r in mf['archives']+mf['assets']]
    names += ['manifest.json', 'defect-momentum-archives.json', 'defect-momentum-jobs.json', 'verification.json']
    releases = json.loads(run(['gh', 'api', 'repos/rjguo1208/Holstein-model/releases?per_page=100']))
    detail = next(r for r in releases if r['tag_name'] == TAG)
    assert detail['target_commitish'] == commit
    present = {a['name']:a for a in detail['assets']}
    wanted = {name:dict(bytes=(OUT/name).stat().st_size, sha256=digest(OUT/name)) for name in names}
    missing = []
    for name in names:
        if name in present:
            assert present[name]['size'] == wanted[name]['bytes']
            assert present[name].get('digest') == 'sha256:'+wanted[name]['sha256']
        else:
            missing.append(str(OUT/name))
    if missing:
        # gh handles upload concurrency within this single release operation.
        subprocess.run(['gh', 'release', 'upload', TAG, *missing], cwd=REPO, check=True)
    detail = json.loads(run(['gh', 'api', f'repos/rjguo1208/Holstein-model/releases/{detail["id"]}']))
    present = {a['name']:a for a in detail['assets']}
    assert set(present) == set(names)
    for name, local in wanted.items():
        assert present[name]['size'] == local['bytes']
        assert present[name].get('digest') == 'sha256:'+local['sha256']
    record = dict(passed=True, release_id=detail['id'], tag=TAG, commit=commit,
                  assets=len(names), uploaded_assets=present,
                  all_server_sizes_and_sha256_match=True)
    (OUT/'upload-verification.json').write_text(json.dumps(record, indent=2)+'\n')
    run(['gh', 'release', 'edit', TAG, '--draft=false', '--latest=false'])
    print('Published verified release', TAG, 'with', len(names), 'assets at', commit, flush=True)


if __name__ == '__main__':
    main()
