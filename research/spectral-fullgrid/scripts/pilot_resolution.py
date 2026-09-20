"""Noise-matched synthetic recovery tests after real-data model selection."""
import concurrent.futures as cf
import fcntl
import json
import os
from pathlib import Path
import sys
import numpy as np
from scipy.signal import find_peaks

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from analysis import read_runs,combine_blocks,estimate_green
from continuation import bin_kernel
from map_continuation import lorentz_map
from resolution import synthetic_spectrum,broaden_truth,select_gaussian
from pilot_diagnostics import public

AXIS=np.linspace(-3.25,5.5,1401)
ETAS=[.1,.25,.5,1.]
CASES=[('single',0.,.06),('double_015',.15,.06),('double_030',.30,.06),
       ('double_060',.60,.06),('double_090',.90,.06),('continuum',0.,1.5)]

def peak_test(axis,curve,truth):
    center=truth['main_center'];delta=truth['separation']
    window=abs(axis-center)<max(.75,delta)
    peaks,_=find_peaks(curve,prominence=max(curve[window])*.05)
    peaks=[i for i in peaks if window[i]]
    if delta==0:return dict(peaks=len(peaks),false_split=len(peaks)>1)
    targets=[center-delta/2,center+delta/2]
    match=[]
    for target in targets:
        near=[i for i in peaks if abs(axis[i]-target)<delta/4]
        match.append(max(near,key=lambda i:curve[i]) if near else None)
    recovered=False;contrast=0.
    if all(i is not None for i in match) and match[0]!=match[1]:
        i,j=sorted(match);contrast=float(1-min(curve[i:j+1])/min(curve[i],curve[j]))
        recovered=contrast>=.1
    return dict(peaks=len(peaks),doublet_recovered=recovered,valley_contrast=contrast)

def run(task):
    lam,ik,name,delta,width=task
    selected=json.loads((ROOT/f'results/comparison/lambda{lam:.2f}_ik{ik:03d}.json').read_text())
    config=selected['primary_cv']['config']
    paths=sorted((ROOT/'results/pilot_data').glob(f'ik{ik:03d}/lambda{lam:.2f}_*/run.json'))
    rows,owners,metas=read_runs([p.parent for p in paths],'rb_blocks.csv');p=metas[0]
    G,C,_=estimate_green(combine_blocks(rows,owners,4),p)
    ev,vec=np.linalg.eigh(C);L=vec*np.sqrt(np.maximum(ev,0))
    truth=synthetic_spectrum(p,delta,width)
    tau=(np.arange(p['bins'])+.5)*p['tau_max']/p['bins']
    trueG=bin_kernel(tau,p['tau_max']/p['bins'],truth['energies'],p['mu'])@truth['quadrature']
    reference=np.array([broaden_truth(AXIS,truth,eta) for eta in ETAS])
    rng=np.random.default_rng(593781+int(lam*100)*10000+ik*100+CASES.index((name,delta,width)))
    spectra=[];records=[];failures=[]
    for rep in range(24):
        # Eight independent synthetic chain means have covariance 8*C each,
        # so their average has exactly the measured absolute covariance C.
        chains=trueG+np.sqrt(8.)*rng.normal(size=(8,len(G)))@L.T
        try:
            fit=select_gaussian(chains,C,p,config)
            curves=np.array([lorentz_map(AXIS,fit['energies'],fit['weights'],eta) for eta in ETAS])
            errors=np.trapezoid(abs(curves-reference),AXIS,axis=1)/np.trapezoid(reference,AXIS,axis=1)
            records.append(dict(replicate=rep,alpha=fit['alpha'],relative_l1=errors,
                peaks=[peak_test(AXIS,curve,truth) for curve in curves]))
            spectra.append(curves)
        except (ValueError,RuntimeError,np.linalg.LinAlgError) as error:
            failures.append(dict(replicate=rep,error=str(error)))
    out=ROOT/'results/resolution'
    result=dict(coupling=lam,index=ik,case=name,truth=truth,configuration=config,
        noise='measured absolute covariance; 8 independent Gaussian chain means; covariance treated as known',
        replicates=24,successful=len(records),failures=failures,records=records,eta=ETAS,
        detection='two maxima within separation/4 of true centers, prominence >=5% and valley >=10%',
        note='conditional synthetic benchmark, not a universal resolution certificate')
    (out/f'lambda{lam:.2f}_ik{ik:03d}_{name}.json').write_text(json.dumps(public(result),indent=2)+'\n')
    np.savez_compressed(out/f'lambda{lam:.2f}_ik{ik:03d}_{name}.npz',energy=AXIS,eta=ETAS,
        truth=reference,reconstructed=np.asarray(spectra),tau=tau,G_true=trueG,covariance=C)
    print(f'synthetic lambda={lam} ik={ik:03d} {name}: {len(records)}/24 completed',flush=True)
    return dict(coupling=lam,index=ik,case=name,successful=len(records),failures=len(failures))

if __name__=='__main__':
    # A second step in the same allocation can use idle cores while the two
    # MaxEnt bootstrap points finish. The serial driver joins this file lock
    # and verifies completion rather than starting a duplicate calculation.
    with (ROOT/'results/resolution.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        out=ROOT/'results/resolution'
        if (out/'summary.json').exists():
            complete=json.loads((out/'summary.json').read_text())
            assert complete['complete'] and len(complete['points'])==36
            print('Verified already completed synthetic stage',flush=True)
        else:
            out.mkdir(exist_ok=False)
            tasks=[(lam,ik,*case) for lam in [.25,.5] for ik in [25,36,39] for case in CASES]
            with cf.ProcessPoolExecutor(max_workers=min(36,int(os.environ.get('SLURM_CPUS_PER_TASK','1')))) as pool:
                results=list(pool.map(run,tasks))
            (out/'summary.json').write_text(json.dumps(dict(complete=True,points=results,
                job_id=os.environ.get('SLURM_JOB_ID')),indent=2)+'\n')
