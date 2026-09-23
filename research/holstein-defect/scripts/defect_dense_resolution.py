"""Controlled single/doublet recovery with old/new measured relative covariances.

This is a diagnostic on prescribed atomic spectra, not proof that any physical
sideband is resolved. It holds representation and selection fixed between
noise levels, uses exact local moments, and never uses VED spectral features.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import sys
import numpy as np
from scipy.signal import find_peaks
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from defect_analysis import public
from defect_dense import prepare,estimate,ALPHAS
from defect_entropy import entropy_scan,local_moments
from continuation import bin_kernel
from spectral_cv import whitening
from map_continuation import lorentz_map


def trial(task):
    lam,site,version,separation,rep,tau,pp,C,truth,seed=task
    rng=np.random.default_rng(seed);ev,vec=np.linalg.eigh(C);ev=np.maximum(ev,0)
    draws=truth+(rng.normal(size=(8,len(tau)))@((vec*np.sqrt(ev*8)).T))
    options=dict(size=641,upper=10.,tau_limit=8.,tolerance=1e-8)
    use=tau<=8;W=whitening(truth[use],(4*C)[np.ix_(use,use)],1e-8);folds=[]
    who=np.arange(8)//2
    for fold in range(4):
        g=draws[who!=fold].mean(axis=0);test=draws[who==fold].mean(axis=0)
        if np.any(g<=0):raise ValueError('Synthetic draw was not positive')
        fs=entropy_scan(tau,g,C*(8/6),pp,ALPHAS,**options)
        if not all(f['converged'] for f in fs):raise ValueError('Synthetic solver did not converge')
        rr=np.array([W@(f['predicted'][use]-test[use]) for f in fs]);folds.append(np.sum(rr**2,axis=1)/len(W))
    folds=np.array(folds);mean=folds.mean(axis=0);best=mean.argmin();se=folds[:,best].std(ddof=1)/2
    index=np.flatnonzero(mean<=mean[best]+se)[-1]
    fs=entropy_scan(tau,draws.mean(axis=0),C,pp,ALPHAS,**options);f=fs[index]
    if not f['converged']:raise ValueError('Synthetic central solver did not converge')
    axis=np.linspace(-2.2,-.4,1801);A=lorentz_map(axis,f['energies'],f['weights'],.075)
    peaks,_=find_peaks(A,prominence=.05*(.35/(np.pi*.075)))
    detected=axis[peaks]
    if separation:
        actual=np.array([-1.3-separation/2,-1.3+separation/2])
        success=len(detected)==2 and np.max(abs(detected-actual))<=.10
    else:success=len(detected)==1 and abs(detected[0]+1.3)<=.10
    return dict(coupling=lam,site=site,noise=version,separation=separation,replicate=rep,
        alpha=f['alpha'],success=bool(success),peaks=detected,false_split=bool(separation==0 and len(detected)>=2))


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--fit',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True);ap.add_argument('--jobs',type=int,default=16)
    ap.add_argument('--draws',type=int,default=16);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    assert json.loads((a.fit/'manifest.json').read_text())['complete']
    tasks=[];inputs={}
    for il,lam in enumerate([.25,.5]):
        for site in [0,8]:
            label=f'lambda{lam:.2f}_i{site:02d}';p=json.loads((a.fit/label/'summary.json').read_text())['parameters']
            for iv,version in enumerate(['old','new']):
                path=(ROOT/'results/defect_spectra_01/analysis'/label/'green.npz') if version=='old' else a.fit/label/'green.npz'
                inputs[str(path)]=sha256(path.read_bytes()).hexdigest();s=np.load(path)
                bins=96 if version=='old' else 192
                tau,b,_,pp=prepare(s['blocks'][:,:2+bins],s['owners'],dict(p,bins=bins),dict(factor=bins//48,block=4))
                G,C=estimate(b,pp);relative=C/np.outer(G,G)
                for si,separation in enumerate([0.,.3,.6,1.]):
                    double=np.array([-1.3-separation/2,-1.3+separation/2]);wp=np.full(2,.175)
                    other=np.array([-2.7,0.,3.]);M=np.vstack([np.ones(3),other,other**2])
                    target=local_moments(p)-np.array([wp.sum(),double@wp,double**2@wp])
                    weights=np.linalg.solve(M,target);assert min(weights)>0
                    en=np.r_[other,double];w=np.r_[weights,wp]
                    truth=bin_kernel(tau,pp['tau_max']/pp['bins'],en,pp['mu'])@w
                    covariance=relative*np.outer(truth,truth)
                    for rep in range(a.draws):tasks.append((lam,site,version,separation,rep,tau,pp,covariance,truth,994000001+il*100000+site*1000+iv*10000+si*100+rep))
    records=[]
    with ProcessPoolExecutor(max_workers=a.jobs) as pool:
        for f in as_completed([pool.submit(trial,t) for t in tasks]):
            records.append(f.result())
            if len(records)%32==0:print('COMPLETE',len(records),'/',len(tasks),flush=True)
    report=dict(complete=True,job_id=os.environ.get('SLURM_JOB_ID'),records=records,inputs=inputs,
        protocol='16 independent Gaussian-noise experiments per spectrum and covariance; 8 simulated chains; 4-fold alpha selection; common dt=0.25, 641 energies, MaxEnt local moments; no VED',
        detector='Broaden at eta=0.075; count peaks in [-2.2,-0.4] with absolute prominence 0.05*0.35/(pi*0.075); doublet succeeds only with exactly two peaks within 0.10 of each prescribed pole; single peak tests false splits.',
        limitations='Finite diagnostic sample, fixed measured covariance, artificial spectra; these recovery fractions are not a certificate of physical peak resolution.',
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    (a.out/'summary.json').write_text(json.dumps(public(report),indent=2)+'\n');shutil.copyfile(Path(__file__),a.out/'resolution_source.py')


if __name__=='__main__':main()
