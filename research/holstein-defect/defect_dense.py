"""Continuation of densely sampled local G; all selection is reference-free.

Four folds each hold out two independent chains. A one-fold-standard-error
heuristic selects alpha, then the lowest selected mean score picks the fixed
representation. All representations use identical validation bins/covariances.
These scores and conditional bootstrap ranges are not resolution certificates.
"""
from __future__ import annotations
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
from defect_analysis import jackknife,normalize,reblock
from defect_entropy import entropy_scan,local_moments
from map_continuation import scan
from continuation import bin_kernel
from spectral_cv import whitening

ALPHAS=np.logspace(-3,4,8)
TIKHONOV=np.r_[0.,np.logspace(-9,3,13)]


def configurations():
    base=dict(method='maxent',factor=2,block=4,size=641,upper=10.,tau_limit=8.,tolerance=1e-8)
    return [dict(base,name=name,**{})|change for name,change in [
        ('dense',{}),('time_fine',dict(factor=1)),('time_long',dict(tau_limit=12.)),
        ('energy_coarse',dict(size=321)),('upper_14',dict(upper=14.)),
        ('block_8',dict(block=8)),('cov_cut_1e-6',dict(tolerance=1e-6)),
        ('cov_cut_1e-10',dict(tolerance=1e-10)),
        ('old_representation',dict(method='tikhonov',factor=4,size=321,tolerance=1e-7))]]


def read_green(paths):
    rows=[];owners=[];metas=[];hashes={}
    keys=['t','omega','g','U','mu','origin','tau_max','bins','observables']
    for i,path in enumerate(map(Path,paths)):
        m=json.loads((path/'run.json').read_text())
        if not m['complete'] or m.get('observables')!='green':raise ValueError('Need completed compact green data')
        if metas and any(m[k]!=metas[0][k] for k in keys):raise ValueError('Incompatible chains')
        if m['arc_cap_attempts'] or m['hop_cap_attempts']:raise ValueError('Order cap encountered')
        b=np.loadtxt(path/'blocks.csv',delimiter=',',skiprows=1)
        if b.shape!=(m['blocks'],m['bins']+3) or not np.isfinite(b).all():raise ValueError('Malformed blocks')
        if not np.array_equal(b[:,0],np.arange(m['blocks'])):raise ValueError('Missing or reordered blocks')
        rows.append(b[:,1:]);owners.extend([i]*len(b));metas.append(m)
        for f in ['run.json','blocks.csv']:hashes[str(path.name+'/'+f)]=sha256((path/f).read_bytes()).hexdigest()
    if len(set(m['seed'] for m in metas))!=len(metas):raise ValueError('Duplicate seeds')
    return np.concatenate(rows),np.array(owners),metas,hashes


def prepare(rows,owners,p,c):
    n=p['bins'];factor=c['factor']
    if n%factor or rows.shape[1]!=n+2:raise ValueError('Unexpected time bins')
    compact=np.column_stack([rows[:,:2],rows[:,2:].reshape(len(rows),n//factor,factor).sum(axis=2)])
    b=reblock(compact,owners,c['block'])
    who=np.concatenate([np.full(np.sum(owners==i)//c['block'],i) for i in np.unique(owners)])
    pp=dict(p,bins=n//factor);tau=(np.arange(pp['bins'])+.5)*pp['tau_max']/pp['bins']
    return tau,b,who,pp


def estimate(b,p):return jackknife(b,lambda total:normalize(total,p))[:2]


def fits(tau,G,C,p,c):
    options={k:c[k] for k in ['size','upper','tau_limit','tolerance']}
    if c['method']=='tikhonov':
        return scan(tau,G,C,p,TIKHONOV,lower_bound=-2*p['t']-p['U']-p['g']**2/p['omega'],moments=local_moments(p),**options)
    result=entropy_scan(tau,G,C,p,ALPHAS,**options)
    # A numerical failure is repaired for the identical objective, not hidden
    # by discarding the parameter value or changing a model-selection score.
    for i,f in enumerate(result):
        if f['converged']:continue
        from defect_entropy_stable import solve
        try:
            recovered=solve(tau,G,C,p,f['alpha'],**options,max_iterations=1200)
        except Exception as error:
            raise ValueError(f"{c['name']} g={p['g']} site={p['origin']} alpha={f['alpha']}: {error}; initial gradient={f['dual_gradient']}") from error
        recovered['initial_solver_failure']=dict(gradient=f['dual_gradient'],iterations=f['iterations'])
        if not recovered['converged']:raise ValueError('Same-objective solver recovery failed')
        result[i]=recovered
    return result


def select_prepared(tau,blocks,who,p,c,validation):
    tt,tb,tw,tp=validation;use=tt<=8.;folds=[]
    if len(np.unique(who))!=8:raise ValueError('Expected eight independent chains')
    for fold in range(4):
        train,TC=estimate(blocks[who//2!=fold],p)
        test,VC=estimate(tb[tw//2==fold],tp)
        W=whitening(test[use],VC[np.ix_(use,use)],1e-8)
        fs=fits(tau,train,TC,p,c)
        K=bin_kernel(tt[use],tp['tau_max']/tp['bins'],fs[0]['energies'],p['mu'])
        pred=np.array([f['weights'] for f in fs])@K.T
        residual=(pred-test[use])@W.T
        folds.append(np.sum(residual**2,axis=1)/len(W))
    folds=np.array(folds);means=folds.mean(axis=0);best=int(np.argmin(means))
    se=float(folds[:,best].std(ddof=1)/2)
    selected=int(np.flatnonzero(means<=means[best]+se)[-1])
    G,C=estimate(blocks,p);fs=fits(tau,G,C,p,c)
    return fs[selected],dict(config=c,alphas=[f['alpha'] for f in fs],folds=folds,means=means,
        selected=selected,minimum=best,fold_standard_error=se,score=float(means[selected]),
        rule='largest alpha within one fold-standard-error of minimum; then lowest selected score across configurations')


def select(rows,owners,p,c):
    validation=prepare(rows,owners,p,dict(c,factor=4,block=4))
    return select_prepared(*prepare(rows,owners,p,c),c,validation)


def bootstrap(rows,owners,p,c,seed):
    # Resample 400k-step groups inside each independent chain. block_8 uses
    # 800k-step groups; validation is rebinned from the same resampled groups.
    _,base,who,pp=prepare(rows,owners,p,dict(c,factor=1))
    rng=np.random.default_rng(seed)
    sampled=np.concatenate([b[rng.integers(len(b),size=len(b))] for b in [base[who==i] for i in range(8)]])
    cc=dict(c,block=1)
    # Samples are already blocked; never block them a second time.
    validation=prepare(sampled,who,pp,dict(cc,factor=4))
    return select_prepared(*prepare(sampled,who,pp,cc),cc,validation)


@lru_cache(maxsize=64)
def broadening_kernel(lower,upper,size,eta):
    en=np.linspace(lower,upper,size);axis=np.linspace(-4.5,5.5,1601)
    return eta/np.pi/((axis[:,None]-en)**2+eta**2)


def curves(fit,etas=(.15,.25,.5,1.)):
    e=fit['energies'];w=fit['weights']
    return np.array([w@broadening_kernel(float(e[0]),float(e[-1]),len(e),eta).T for eta in etas])
