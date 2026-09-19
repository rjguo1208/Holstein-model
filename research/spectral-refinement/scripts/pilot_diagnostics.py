"""Check raw block statistics before changing the continuation model."""
import concurrent.futures as cf
import json
import os
from pathlib import Path
import sys
import subprocess
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from analysis import read_runs, combine_blocks, estimate_green
from map_continuation import rebin_conditional
from spectral_cv import whitening

def statistic(a,b,Ca,Cb,tolerance=1e-8):
    W=whitening((a+b)/2,Ca+Cb,tolerance)
    r=W@(a-b)
    return dict(chi2_per_mode=float(r@r/len(W)),rank=len(W),
                max_marginal_z=float(np.max(abs(a-b)/np.sqrt(np.diag(Ca+Cb)))))

def diagnose(task):
    label,data,lam,ik=task
    paths=sorted(Path(data).glob(f'ik{ik:03d}/lambda{lam:.2f}_*/run.json'))
    rows,owners,metas=read_runs([p.parent for p in paths],'rb_blocks.csv')
    p=metas[0]
    tau,rows,p=rebin_conditional(rows,p,2)
    use=tau<=8
    estimates={}
    for factor in [1,2,4,8,16]:
        grouped=combine_blocks(rows,owners,factor)
        g,c,_=estimate_green(grouped,p)
        ev=np.linalg.eigvalsh(c[np.ix_(use,use)]/np.outer(g[use],g[use]))
        estimates[factor]=(g,c)
    g,c=estimates[4]
    chain=[]; chain_cov=[]; halves=[]
    for owner in np.unique(owners):
        raw=rows[owners==owner]
        group=combine_blocks(raw,np.zeros(len(raw),int),4)
        cg,cc,_=estimate_green(group,p)
        chain.append(cg);chain_cov.append(cc)
        n=len(group)//2
        a,ac,_=estimate_green(group[:n],p);b,bc,_=estimate_green(group[n:],p)
        halves.append(statistic(a[use],b[use],ac[np.ix_(use,use)],bc[np.ix_(use,use)]))
    chain=np.asarray(chain);chain_cov=np.asarray(chain_cov)
    chain_cov_of_mean=np.cov(chain,rowvar=False,ddof=1)/len(chain)
    block=[]
    for factor,(gg,cc) in estimates.items():
        ev=np.linalg.eigvalsh(cc[np.ix_(use,use)]/np.outer(gg[use],gg[use]))
        ratio=np.sqrt(np.diag(cc)[use]/np.diag(c)[use])
        block.append(dict(factor=factor,blocks=len(rows)//factor,
            median_sigma_over_factor4=float(np.median(ratio)),max_sigma_over_factor4=float(ratio.max()),
            relative_error=np.sqrt(np.diag(cc))/gg,
            eigenvalues=ev,ranks={str(t):int(np.sum(ev>ev[-1]*t)) for t in [1e-6,1e-8,1e-10]}))
    folds=[]
    for owner in np.unique(owners):
        train=combine_blocks(rows[owners!=owner],owners[owners!=owner],4)
        a,ac,_=estimate_green(train,p)
        b,bc=chain[owner],chain_cov[owner]
        folds.append(statistic(a[use],b[use],ac[np.ix_(use,use)],bc[np.ix_(use,use)]))
    # Integrated autocorrelation of the ratio's linear influence, in raw blocks.
    influences=rows[:,3:]-rows[:,0,None]*(rows[:,3:].sum(axis=0)/rows[:,0].sum())
    autocorrelation=[]
    for target in [.125,.625,2.125,6.125,10.125]:
        j=int(np.argmin(abs(tau-target)));times=[]
        for owner in np.unique(owners):
            x=influences[owners==owner,j];x=x-x.mean()
            ac=np.correlate(x,x,mode='full')[len(x)-1:]/np.dot(x,x)
            total=0.
            for q in range(1,min(len(x)-1,100),2):
                pair=ac[q]+ac[q+1]
                if pair<=0:break
                total+=pair
            times.append(float(max(1.,1+2*total)))
        autocorrelation.append(dict(tau=float(tau[j]),inefficiency_in_raw_blocks=times))
    result=dict(dataset=label,coupling=lam,index=ik,k=p['k'],chains=len(paths),
        production_steps=sum(m['blocks']*m['steps_per_block'] for m in metas),
        tau=tau,G=g,covariance=c,reblocking=block,
        chain_means=chain,chain_covariance_of_mean=chain_cov_of_mean,
        chain_sigma_over_block_sigma=np.sqrt(np.diag(chain_cov_of_mean)/np.diag(c)),
        chain_agreement=folds,half_chain_drift=halves,autocorrelation=autocorrelation)
    print(f'{label} lambda={lam} ik={ik:03d}: chain score max={max(f["chi2_per_mode"] for f in folds):.3g}; '
          f'block16/4 sigma median={block[-1]["median_sigma_over_factor4"]:.3g}',flush=True)
    return result

def public(x):
    if isinstance(x,np.ndarray):return public(x.tolist())
    if isinstance(x,np.generic):return public(x.item())
    if isinstance(x,float) and not np.isfinite(x):return None
    if isinstance(x,dict):return {k:public(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [public(v) for v in x]
    return x

if __name__=='__main__':
    plan=json.loads((ROOT/'plans/pilot.json').read_text())
    tasks=[('baseline',plan['baseline'],lam,ik) for lam in plan['couplings'] for ik in range(41)]
    tasks += [('pilot',str(ROOT/'results/pilot_data'),lam,ik) for lam in plan['couplings'] for ik in plan['indices']]
    with cf.ProcessPoolExecutor(max_workers=min(32,int(os.environ.get('SLURM_CPUS_PER_TASK','1')))) as pool:
        results=list(pool.map(diagnose,tasks))
    comparisons=[]
    for lam in plan['couplings']:
        for ik in plan['indices']:
            a=next(r for r in results if (r['dataset'],r['coupling'],r['index'])==('baseline',lam,ik))
            b=next(r for r in results if (r['dataset'],r['coupling'],r['index'])==('pilot',lam,ik))
            use=a['tau']<=8
            s=statistic(a['G'][use],b['G'][use],a['covariance'][np.ix_(use,use)],b['covariance'][np.ix_(use,use)])
            s.update(coupling=lam,index=ik,sigma_ratio_pilot_over_baseline=np.sqrt(np.diag(b['covariance'])/np.diag(a['covariance'])))
            comparisons.append(s)
    out=ROOT/'results/diagnostics';out.mkdir(exist_ok=False)
    (out/'summary.json').write_text(json.dumps(public(dict(points=results,comparisons=comparisons)),indent=2)+'\n')
    # The single allocation carries the remaining independent pilot stages.
    # No reference spectrum is opened by any of these programs.
    subprocess.run([sys.executable,'-m','pytest','-q','tests'],cwd=ROOT,check=True)
    for stage in ['pilot_compare.py','pilot_bootstrap.py','pilot_repair_bootstrap.py','pilot_resolution.py']:
        print('Starting '+stage,flush=True)
        subprocess.run([sys.executable,str(ROOT/'scripts'/stage)],cwd=ROOT,check=True)
