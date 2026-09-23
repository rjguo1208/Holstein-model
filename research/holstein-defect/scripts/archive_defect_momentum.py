"""Create bounded reproducible archives and direct-download map/figure assets."""
from __future__ import annotations
import argparse
import ast
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parents[1]
TAG = 'defect-momentum-20260923'
URL = 'https://github.com/rjguo1208/Holstein-model/releases/download/'+TAG+'/'


def digest(path):
    h = sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def sources():
    paths = {ROOT/'Makefile', ROOT/'requirements.txt', ROOT/'DEFECT_MOMENTUM.md', ROOT/'scripts/python.sh'}
    paths.update((ROOT/'src').glob('*'))
    paths.update((ROOT/'tests').glob('*.cpp'))
    paths.update((ROOT/'tests').glob('*.py'))
    paths.update((ROOT/'tests/data').glob('*'))
    paths.update((ROOT/'scripts').glob('*defect_momentum*'))
    todo = [p for p in paths if p.suffix == '.py']
    while todo:
        p = todo.pop()
        for node in ast.walk(ast.parse(p.read_text())):
            names = [node.module] if isinstance(node, ast.ImportFrom) else [a.name for a in node.names] if isinstance(node, ast.Import) else []
            for name in names:
                if not name:
                    continue
                for parent in [ROOT, ROOT/'scripts']:
                    local = parent/(name.split('.')[0]+'.py')
                    if local.is_file() and local not in paths:
                        paths.add(local)
                        todo.append(local)
    return {p for p in paths if p.is_file()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out = args.out.resolve()
    for directory, name in [('defect_momentum_mc_report_01', 'summary.json'),
                             ('defect_momentum_ved_report_02', 'summary.json'),
                             ('defect_momentum_fit_01', 'manifest.json')]:
        assert json.loads((ROOT/'results'/directory/name).read_text())['complete']
    args.out.mkdir(parents=True, exist_ok=False)
    job_ids = ['20883897', '20883985', '20884137', '20884064', '20884069', '20884108', '20884144', '20884149']
    result = subprocess.run(['sacct', '-P', '-n', '-j', ','.join(job_ids),
                             '--format=JobID,JobName,Partition,State,ExitCode,Elapsed,NCPUS,MaxRSS'],
                            text=True, capture_output=True, check=True)
    rows = [line.split('|') for line in result.stdout.splitlines() if line.strip()]
    main_rows = [r for r in rows if r[0] in job_ids]
    assert len(main_rows) == len(job_ids) and all(r[3] == 'COMPLETED' and r[4] == '0:0' for r in main_rows)
    jobs = dict(complete=True, jobs=[dict(zip(['job_id','name','partition','state','exit_code','elapsed','cpus','max_rss'], r)) for r in rows],
                ordering='VED report complete, then analytic nonlocal DiagMC validation, production sampling, continuation and comparison',
                archive_job_id=os.environ.get('SLURM_JOB_ID'))
    (args.out/'defect-momentum-jobs.json').write_text(json.dumps(jobs, indent=2)+'\n')
    paths = sources()
    roots = ['defect_momentum_ved_01', 'defect_momentum_ved_report_02',
             'defect_momentum_mc_pilot_01', 'defect_momentum_mc_data_01',
             'defect_momentum_fit_pilot_01', 'defect_momentum_fit_01', 'defect_momentum_mc_report_01']
    for name in roots:
        paths.update(p for p in (ROOT/'results'/name).rglob('*')
                     if p.is_file() and '__pycache__' not in p.parts and p.suffix != '.pyc')
    paths.add(ROOT/'results/defect_momentum_mc_validation_01.json')
    # Keep the original summary/source needed by acquisition's VED-first guard.
    paths.add(ROOT/'results/defect_momentum_ved_report_01/summary.json')
    paths.update(p for p in (ROOT/'results/defect_momentum_ved_report_01/source').rglob('*') if p.is_file())
    for job in job_ids:
        paths.update((ROOT/'logs').glob('*-'+job+'.out'))
    paths.add(args.out/'defect-momentum-jobs.json')
    ordered = sorted(paths, key=lambda p: str(p.relative_to(ROOT)))
    manifest = dict(complete=False, files=[], archives=[], assets=[], shared_archive_root='holstein-diagmc/',
                    archive_job_id=os.environ.get('SLURM_JOB_ID'))
    batches, current, size = [], [], 0
    limit = 80*1024**2
    for path in ordered:
        n = path.stat().st_size
        if n > limit:
            raise ValueError(f'One file exceeds the archive input cap: {path}, {n}')
        if current and size+n > limit:
            batches.append(current)
            current, size = [], 0
        current.append(path)
        size += n
    if current:
        batches.append(current)
    for i, batch in enumerate(batches):
        dest = args.out/f'defect-momentum-part{i+1:02d}.zip'
        with zipfile.ZipFile(dest, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
            for path in batch:
                relative = 'holstein-diagmc/'+path.relative_to(ROOT).as_posix()
                z.write(path, relative)
                manifest['files'].append(dict(path=relative, bytes=path.stat().st_size,
                                              sha256=digest(path), archive=dest.name))
        assert dest.stat().st_size < 95*1024**2
        manifest['archives'].append(dict(name=dest.name, bytes=dest.stat().st_size,
                                         sha256=digest(dest), url=URL+dest.name))
        print('ARCHIVE', dest.name, dest.stat().st_size, flush=True)
    # Direct assets avoid extracting archives simply to inspect a plot/map.
    direct = []
    for name in ['defect_momentum_ved_report_02', 'defect_momentum_mc_report_01']:
        directory = ROOT/'results'/name
        direct.extend(p for p in directory.iterdir() if p.suffix in ['.png', '.svg', '.pdf'])
    direct += list((ROOT/'results/defect_momentum_ved_report_02').glob('*Nh12_R96_W32.npz'))
    direct += list((ROOT/'results/defect_momentum_mc_report_01').glob('*.npz'))
    for path in direct:
        dest = args.out/path.name
        assert not dest.exists()
        shutil.copyfile(path, dest)
        manifest['assets'].append(dict(name=dest.name, bytes=dest.stat().st_size,
                                       sha256=digest(dest), url=URL+dest.name))
    manifest['complete'] = True
    (args.out/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    index = dict(complete=True, release='https://github.com/rjguo1208/Holstein-model/releases/tag/'+TAG,
                 archives=manifest['archives'], assets=manifest['assets'], files=len(manifest['files']),
                 manifest_sha256=digest(args.out/'manifest.json'),
                 total_bytes=sum(r['bytes'] for r in manifest['archives']+manifest['assets']))
    (args.out/'defect-momentum-archives.json').write_text(json.dumps(index, indent=2)+'\n')
    print('COMPLETE', len(manifest['archives']), 'archives;', len(manifest['assets']), 'direct assets', flush=True)


if __name__ == '__main__':
    main()
