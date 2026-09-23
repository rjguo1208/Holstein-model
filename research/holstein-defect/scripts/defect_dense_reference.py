"""Independent local VED references after continuation has been frozen."""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
from functools import lru_cache
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import sys
import numpy as np
from scipy.linalg import eigh_tridiagonal
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from defect_reference import Cloud,DefectBasis,Parameters,lanczos
from defect_entropy import local_moments
from map_continuation import lorentz_map
from defect_analysis import public


@lru_cache(maxsize=2)
def cloud(nh):return Cloud(nh)


def point(task):
    out,lam,nh,radius,site=task
    basis=DefectBasis(nh,radius,cloud=cloud(nh));p=Parameters(g=np.sqrt(2*lam),U=1.)
    en,w,a,b=lanczos(basis.hamiltonian(p),basis.source(site),400)
    expected=local_moments(dict(t=1.,g=p.g,U=1.,origin=site))
    residual=np.array([w.sum(),en@w,en**2@w])-expected
    if max(abs(residual))>1e-8:raise ValueError('VED moment check failed')
    e2,v2=eigh_tridiagonal(a[:200],b[:199]);w2=v2[0]**2
    axis=np.linspace(-4.5,5.5,1601);etas=np.array([.15,.25,.5,1.])
    A=np.array([lorentz_map(axis,en,w,eta) for eta in etas]);A2=np.array([lorentz_map(axis,e2,w2,eta) for eta in etas])
    path=out/f'lambda{lam:.2f}_nh{nh}_R{radius}_i{site:02d}.npz'
    np.savez_compressed(path,energies=en,weights=w,alpha=a,beta=b,energy=axis,eta=etas,A=A,A_200=A2)
    return dict(coupling=lam,generations=nh,radius=radius,site=site,dimension=basis.dim,
        moment_residuals=residual,steps_relative_l1=np.trapezoid(abs(A-A2),axis,axis=1)/np.trapezoid(A,axis,axis=1),
        file=path.name,sha256=sha256(path.read_bytes()).hexdigest())


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--fit',required=True,type=Path)
    ap.add_argument('--out',required=True,type=Path);ap.add_argument('--jobs',type=int,default=8);a=ap.parse_args()
    mf=json.loads((a.fit/'manifest.json').read_text());assert mf['complete'] and len(mf['points'])==34
    a.out.mkdir(parents=True,exist_ok=False);source=a.out/'source';source.mkdir()
    files=[ROOT/'defect_reference.py',ROOT/'defect_entropy.py',Path(__file__)]
    for f in files:shutil.copyfile(f,source/f.name)
    record=dict(complete=False,job_id=os.environ.get('SLURM_JOB_ID'),records=[],
        frozen_manifest_sha256=sha256((a.fit/'manifest.json').read_bytes()).hexdigest(),
        frozen_selection_sha256={p['selected_sha256']:f"lambda{p['coupling']:.2f}_i{p['site']:02d}" for p in mf['points']},
        source_sha256={f.name:sha256(f.read_bytes()).hexdigest() for f in files},
        settings=[dict(generations=nh,radius=r,steps=400) for nh,r in [(10,40),(12,40),(12,64)]])
    def save(): (a.out/'manifest.json').write_text(json.dumps(public(record),indent=2)+'\n')
    save()
    tasks=[(a.out,lam,nh,r,site) for nh,r in [(10,40),(12,40),(12,64)] for site in range(17) for lam in [.25,.5]]
    with ProcessPoolExecutor(max_workers=a.jobs) as pool:
        for future in as_completed([pool.submit(point,t) for t in tasks]):
            r=future.result();record['records'].append(r);save();print('COMPLETE',len(record['records']),r['file'],flush=True)
    record['complete']=True;save()


if __name__=='__main__':main()
