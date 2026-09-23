"""Acquire local G, continue it without VED priors, then compare local spectra."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from defect_analysis import public, read, reblock, jackknife, normalize, summarize
from defect_reference import DefectBasis, Parameters, lanczos
from map_continuation import scan, lorentz_map
from spectral_cv import whitening


def select(tau,blocks,owners,p,alphas):
    def estimate(b):return jackknife(b,lambda total:normalize(total,p))[:2]
    G,C=estimate(blocks)
    V=-p['U'] if p['origin']==0 else 0.
    options=dict(size=321,lower_bound=-2*p['t']-p['U']-p['g']**2/p['omega'],
                 moments=[1.,V,V*V+2*p['t']**2+p['g']**2],tolerance=1e-7)
    use=tau<=8
    scores=[]
    for owner in np.unique(owners):
        train,TC=estimate(blocks[owners!=owner]);test,VC=estimate(blocks[owners==owner])
        W=whitening(test[use],VC[np.ix_(use,use)],1e-7)
        if len(W)<3:raise ValueError('Insufficient validation covariance rank')
        fits=scan(tau,train,TC,p,alphas,**options)
        scores.append([float(np.sum((W@(f['predicted'][use]-test[use]))**2)/len(W)) for f in fits])
    index=int(np.argmin(np.mean(scores,axis=0)))
    fits=scan(tau,G,C,p,alphas,**options)
    return fits[index],dict(alphas=alphas,folds=scores,selected=index,boundary=index in [0,len(alphas)-1],
        rule='minimum mean independent held-out-chain score; no VED spectral or pole prior')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',required=True,type=Path)
    parser.add_argument('--ground',required=True,type=Path)
    parser.add_argument('--jobs',type=int,default=4)
    parser.add_argument('--blocks',type=int,default=160)
    parser.add_argument('--bootstrap',type=int,default=8)
    args=parser.parse_args()
    ground=json.loads((args.ground/'summary.json').read_text())
    if not ground['validation_passed'] or not ground['first_point_passed']:raise ValueError('Ground-state validation gate has not passed')
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    sources=[p for folder in ['src','scripts','tests'] for p in (ROOT/folder).glob('*') if p.is_file()]+list(ROOT.glob('*.py'))+[ROOT/'Makefile']
    manifest=dict(complete=False,job_id=os.environ.get('SLURM_JOB_ID'),sources={},tasks=[],inputs={},
        spectral_resolution_validated=False,model='H_clean - n_0',measured_sites=list(range(9)),reflection_used=True)
    for p in sources:
        dest=out/'source'/p.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
        manifest['sources'][str(p.relative_to(ROOT))]=sha256(p.read_bytes()).hexdigest()
    executable=out/'source/defect_diagmc';shutil.copyfile(ROOT/'build/defect_diagmc',executable);executable.chmod(0o755)
    manifest['executable_sha256']=sha256(executable.read_bytes()).hexdigest()
    def save(): (out/'manifest.json').write_text(json.dumps(public(manifest),indent=2)+'\n')
    save()
    def run(task):
        name,p=task;path=out/'chains'/name;path.parent.mkdir(exist_ok=True)
        command=[str(executable)]
        for key,value in p.items():command+=['--'+key,str(value)]
        command+=['--out',str(path)]
        with (path.parent/(path.name+'.log')).open('w') as stream:subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,check=True)
        return path
    axis=np.linspace(-4.5,5.5,801);etas=np.array([.25,.5,1.]);x=np.arange(-8,9)
    alphas=np.r_[0.,np.logspace(-9,3,13)];summary={}
    for il,lam in enumerate([.25,.5]):
        baseline=ground['cases'][f'lambda{lam:.2f}_U1']['reference']
        maps,low,high,records=[],[],[],[]
        for site in range(9):
            tasks=[]
            for chain in range(4):
                p=dict(t=1,omega=1,g=np.sqrt(2*lam),U=1,mu=baseline['E']-.10,origin=site,
                    **{'tau-max':12,'bins':96,'blocks':args.blocks,'steps-per-block':100000,'warmup':300000,'thin':100,'radius':24,'seed':810_000_001+10_000_000*il+100_000*site+1009*chain})
                tasks.append((f'lambda{lam:.2f}_i{site:02d}_s{chain:02d}',p))
            manifest['tasks'] += [dict(name=n,parameters=p) for n,p in tasks];save()
            with ThreadPoolExecutor(max_workers=args.jobs) as pool: paths=list(pool.map(run,tasks))
            label=f'lambda{lam:.2f}_i{site:02d}'
            info=summarize(paths,out/'analysis'/label,tail=False)
            if info['arc_cap_attempts'] or info['hop_cap_attempts']:raise ValueError('Order cap encountered')
            raw,_,owners,metas=read(paths);p=metas[0]
            # Reblock within each chain and combine pairs of adjacent time bins.
            compact=np.column_stack([raw[:,:2],raw[:,2:2+p['bins']].reshape(len(raw),-1,2).sum(axis=2)])
            blocks=reblock(compact,owners,4);owners=np.repeat(np.arange(4),len(blocks)//4)
            p=dict(p,bins=p['bins']//2);tau=(np.arange(p['bins'])+.5)*p['tau_max']/p['bins']
            fit,cv=select(tau,blocks,owners,p,alphas)
            A=np.array([lorentz_map(axis,fit['energies'],fit['weights'],eta) for eta in etas])
            rng=np.random.default_rng(981_011+il*100+site);boot=[]
            for _ in range(args.bootstrap):
                sampled=np.concatenate([b[rng.integers(len(b),size=len(b))] for b in [blocks[owners==c] for c in range(4)]])
                bf,_=select(tau,sampled,owners,p,alphas)
                boot.append([lorentz_map(axis,bf['energies'],bf['weights'],eta) for eta in etas])
            lower,upper=np.quantile(boot,[.16,.84],axis=0)
            maps.append(A);low.append(lower);high.append(upper)
            inputfile=out/'analysis'/label/'green.npz'
            manifest['inputs'][str(inputfile.relative_to(out))]=sha256(inputfile.read_bytes()).hexdigest()
            records.append(dict(site=site,input=info,fit=fit,cv=cv))
            print(f'CONTINUATION lambda={lam} site={site} CV={np.mean(cv["folds"],axis=0)[cv["selected"]]:.3g}',flush=True)
        reflection=abs(x)
        np.savez_compressed(out/f'lambda{lam:.2f}_map.npz',x=x,measured_sites=np.arange(9),energy=axis,eta=etas,
            A=np.stack(maps,axis=1)[:,reflection],bootstrap_16=np.stack(low,axis=1)[:,reflection],
            bootstrap_84=np.stack(high,axis=1)[:,reflection],coupling=lam,U=1.,reflection_used=True)
        # Only after selection: generate the independent VED LDOS reference.
        reference=[];convergence=[];moments=[]
        for nh,radius in [(10,24),(12,24),(12,40)]:
            basis=DefectBasis(nh,radius);model=Parameters(g=np.sqrt(2*lam),U=1)
            h=basis.hamiltonian(model);curves=[]
            for site in range(9):
                e,w,a,b=lanczos(h,basis.source(site),400)
                V=-1. if site==0 else 0.
                residual=np.array([w.sum(),e@w,e**2@w])-np.array([1.,V,V*V+2+2*lam])
                if max(abs(residual))>1e-9:raise ValueError('Reference local moments failed')
                moments.append(dict(generations=nh,radius=radius,site=site,residual=residual))
                curves.append([lorentz_map(axis,e,w,eta) for eta in etas])
                np.savez_compressed(out/f'reference_lambda{lam:.2f}_nh{nh}_R{radius}_i{site:02d}.npz',energies=e,weights=w,alpha=a,beta=b)
                # Compare 200/400 recursion steps at both cutoff levels.
                from scipy.linalg import eigh_tridiagonal
                e2,v2=eigh_tridiagonal(a[:200],b[:199]);w2=v2[0]**2
                convergence.append(dict(generations=nh,radius=radius,site=site,step_error=float(np.trapezoid(abs(curves[-1][0]-lorentz_map(axis,e2,w2,.25)),axis)/np.trapezoid(curves[-1][0],axis))))
            reference.append(np.stack(curves,axis=1))
        ref=reference[-1]
        np.savez_compressed(out/f'reference_lambda{lam:.2f}_map.npz',x=x,energy=axis,eta=etas,A=ref[:,reflection],coupling=lam,U=1.)
        errors=[]
        for j,eta in enumerate(etas):
            errors.append(dict(eta=float(eta),relative_l1=(np.trapezoid(abs(ref[j]-np.stack(maps,axis=1)[j]),axis,axis=1)/np.trapezoid(ref[j],axis,axis=1)).tolist()))
        cloud_error=np.trapezoid(abs(reference[1]-reference[0]),axis,axis=2)/np.trapezoid(reference[1],axis,axis=2)
        spatial_error=np.trapezoid(abs(reference[2]-reference[1]),axis,axis=2)/np.trapezoid(reference[2],axis,axis=2)
        summary[str(lam)]=dict(records=records,comparison=errors,reference_cloud_relative_l1=cloud_error,
            reference_spatial_relative_l1=spatial_error,reference_steps=convergence,
            reference_moments=moments,reference_radius=40,
            reference_checks_passed=bool(np.max(cloud_error)<.01 and np.max(spatial_error)<.01 and max(c['step_error'] for c in convergence)<.001),
            spectral_resolution_validated=False)
        (out/'summary.json').write_text(json.dumps(public(summary),indent=2)+'\n')
    manifest['complete']=True;save()


if __name__=='__main__':main()
