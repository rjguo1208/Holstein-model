"""Parallel selection and conditional bootstrap for the new single-defect data."""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
from datetime import datetime,timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import sys
import time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from defect_analysis import public
from defect_dense import read_green,configurations,select,bootstrap,curves,prepare,estimate


def write_json(path,obj):
    path=Path(path);temp=path.with_suffix('.tmp');temp.write_text(json.dumps(public(obj),indent=2)+'\n');temp.replace(path)


def point(task):
    data,out,lam,site,nboot=task;start=time.monotonic();label=f'lambda{lam:.2f}_i{site:02d}'
    dest=out/label;dest.mkdir(exist_ok=False)
    paths=[data/'chains'/f'{label}_s{i:02d}' for i in range(8)]
    rows,owners,metas,hashes=read_green(paths);p=metas[0]
    tau,blocks,who,pp=prepare(rows,owners,p,dict(factor=1,block=4));G,C=estimate(blocks,pp)
    np.savez_compressed(dest/'green.npz',tau=tau,G=G,covariance=C,blocks=rows,owners=owners)
    choices=[];spectra=[];selections=[]
    for c in configurations():
        fit,cv=select(rows,owners,p,c)
        choices.append(cv);spectra.append(curves(fit));selections.append(fit)
    chosen=int(np.argmin([c['score'] for c in choices]));c=configurations()[chosen];fit=selections[chosen]
    # The chosen representation and alpha are frozen before any VED reference.
    record=dict(coupling=lam,site=site,parameters=p,input_sha256=hashes,selected=chosen,
        config=c,fit=fit,cv=choices,production_steps=sum(m['blocks']*m['steps_per_block'] for m in metas),
        seeds=[m['seed'] for m in metas],order_caps=sum(m['arc_cap_attempts']+m['hop_cap_attempts'] for m in metas),
        maximum_arcs=max(m['maximum_arcs'] for m in metas),maximum_hops=max(m['maximum_hops'] for m in metas),
        maximum_extent=max(m['maximum_extent'] for m in metas),reference_used=False)
    write_json(dest/'selection.json',record)
    boot=[];bootinfo=[]
    for i in range(nboot):
        bf,cv=bootstrap(rows,owners,p,c,983000001+int(lam*100)*100000+site*1000+i)
        boot.append(curves(bf));bootinfo.append(dict(alpha=bf['alpha'],cv=cv['score']))
    errors=[]
    for factor in [1,2,4,8,16]:
        _,b,_,bp=prepare(rows,owners,p,dict(factor=1,block=factor));_,cov=estimate(b,bp)
        errors.append(np.sqrt(np.diag(cov)))
    chain_blocks=np.array([rows[owners==i].sum(axis=0) for i in range(8)])
    _,chain_cov=estimate(chain_blocks,p);errors.append(np.sqrt(np.diag(chain_cov)))
    np.savez_compressed(dest/'spectra.npz',energy=np.linspace(-4.5,5.5,1601),eta=np.array([.15,.25,.5,1.]),
        selected=curves(fit),candidates=np.array(spectra),bootstrap=np.array(boot),
        bootstrap_quantiles=np.quantile(boot,[.16,.5,.84],axis=0),G=G,tau=tau,
        blocking_errors=np.array(errors),covariance=C)
    record.update(complete=True,bootstrap_draws=nboot,bootstrap=bootinfo,elapsed_seconds=time.monotonic()-start,
        bootstrap_note='16-84% pointwise empirical range, conditional on selected representation; alpha reselected in each draw; does not include representation bias')
    write_json(dest/'summary.json',record)
    return dict(coupling=lam,site=site,configuration=c['name'],score=choices[chosen]['score'],
        alpha=fit['alpha'],seconds=record['elapsed_seconds'],selected_sha256=sha256((dest/'selection.json').read_bytes()).hexdigest())


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--data',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True);ap.add_argument('--jobs',type=int,default=16)
    ap.add_argument('--bootstrap',type=int,default=64);ap.add_argument('--sites',type=int,nargs='*',default=list(range(17)))
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    manifest=json.loads((a.data/'manifest.json').read_text())
    if not manifest['complete']:raise ValueError('Sampling is not complete')
    source=a.out/'source';source.mkdir()
    files=list(ROOT.glob('*.py'))+list((ROOT/'scripts').glob('defect_dense*'))+[ROOT/'scripts/python.sh']
    hashes={}
    for p in files:
        target=source/p.relative_to(ROOT);target.parent.mkdir(exist_ok=True,parents=True);shutil.copyfile(p,target)
        hashes[str(p.relative_to(ROOT))]=sha256(p.read_bytes()).hexdigest()
    record=dict(complete=False,job_id=os.environ.get('SLURM_JOB_ID'),started_utc=datetime.now(timezone.utc).isoformat(),
        source_hashes=hashes,input_manifest_sha256=sha256((a.data/'manifest.json').read_bytes()).hexdigest(),
        points=[],bootstrap_draws_per_point=a.bootstrap,sites=a.sites,workers=a.jobs,reference_used=False,
        spectral_resolution_validated=False)
    write_json(a.out/'manifest.json',record)
    with ProcessPoolExecutor(max_workers=a.jobs) as pool:
        futures=[pool.submit(point,(a.data.resolve(),a.out.resolve(),lam,site,a.bootstrap)) for site in a.sites for lam in [.25,.5]]
        for future in as_completed(futures):
            r=future.result();record['points'].append(r);write_json(a.out/'manifest.json',record)
            print('COMPLETE',json.dumps(r),flush=True)
    record.update(complete=True,finished_utc=datetime.now(timezone.utc).isoformat());write_json(a.out/'manifest.json',record)


if __name__=='__main__':main()
