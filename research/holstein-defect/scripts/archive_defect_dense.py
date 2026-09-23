"""Create immutable, size-limited scientific release archives with file hashes."""
from __future__ import annotations
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
import sys
from zipfile import ZipFile,ZIP_DEFLATED
ROOT=Path(__file__).resolve().parents[1]


def archive(destination,files):
    files=sorted(set(Path(p) for p in files));contents={}
    with ZipFile(destination,'x',compression=ZIP_DEFLATED,compresslevel=6) as z:
        for p in files:
            if not p.is_file():continue
            name='holstein-diagmc/'+str(p.relative_to(ROOT));z.write(p,name)
            contents[name]=dict(bytes=p.stat().st_size,sha256=sha256(p.read_bytes()).hexdigest())
    with ZipFile(destination) as z:
        assert z.testzip() is None
        assert set(z.namelist())==set(contents)
    size=destination.stat().st_size
    assert size<95_000_000,(destination,size)
    return dict(name=destination.name,bytes=size,sha256=sha256(destination.read_bytes()).hexdigest(),files=contents)


def allfiles(path):return [p for p in path.rglob('*') if p.is_file() and '__pycache__' not in str(p) and '.pytest_cache' not in str(p)]


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',type=Path,default=ROOT/'results/defect_dense_release_01');a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    report=json.loads((ROOT/'results/defect_dense_report_01/summary.json').read_text());assert report['complete']
    records=[]
    source=list(ROOT.glob('*.py'))+[ROOT/p for p in ['Makefile','requirements.txt','pytest.ini','README.md','DEFECT.md','DEFECT_DENSE.md']]
    for folder in ['src','tests','scripts']:source+=allfiles(ROOT/folder)
    source+=allfiles(ROOT/'results/defect_optimization_01')
    source+=[ROOT/'results/defect_dense_job_records_01.json',ROOT/'results/defect_scan_01/summary.json',ROOT/'results/defect_spectra_01/summary.json']
    for folder in ['defect_dense_data_01','defect_dense_fit_01']:
        d=ROOT/'results'/folder;source+=[p for p in d.glob('*') if p.is_file()];source+=allfiles(d/'source')
    for folder in ['defect_dense_reference_01','defect_dense_report_01','defect_dense_resolution_01']:
        source+=allfiles(ROOT/'results'/folder)
    source+=list((ROOT/'logs').glob('defect-*20882*.out'))
    records.append(archive(a.out/'defect-dense-core-02.zip',source))
    records.append(archive(a.out/'defect-dense-previous-green-02.zip',allfiles(ROOT/'results/defect_spectra_01/analysis')))
    for lam in [.25,.5]:
        chains=ROOT/'results/defect_dense_data_01/chains'
        files=[p for d in chains.glob(f'lambda{lam:.2f}_*') for p in (allfiles(d) if d.is_dir() else [d])]
        records.append(archive(a.out/f'defect-dense-raw-lambda{lam:.2f}-02.zip',files))
        fit=ROOT/'results/defect_dense_fit_01'
        dirs=sorted(fit.glob(f'lambda{lam:.2f}_i*'));assert len(dirs)==17
        records.append(archive(a.out/f'defect-dense-green-lambda{lam:.2f}-02.zip',[d/'green.npz' for d in dirs]))
        records.append(archive(a.out/f'defect-dense-fits-lambda{lam:.2f}-02.zip',[p for d in dirs for p in d.iterdir() if p.name!='green.npz']))
    url='https://github.com/rjguo1208/Holstein-model/releases/download/defect-spectra-dense-20260923/'
    for r in records:r['url']=url+r['name']
    index=dict(root='holstein-diagmc',complete=True,archives=records)
    (a.out/'manifest.json').write_text(json.dumps(index,indent=2)+'\n')
    public=dict(index,archives=[{k:v for k,v in r.items() if k!='files'} for r in records])
    (a.out/'defect-dense-archives.json').write_text(json.dumps(public,indent=2)+'\n')
    print(json.dumps(public,indent=2))


if __name__=='__main__':main()
