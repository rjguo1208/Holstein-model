"""Normalization, covariance and projector observables for real-space diagrams."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import i0e
from analysis import bin_log_factor


def public(value):
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, dict): return {k:public(v) for k,v in value.items()}
    if isinstance(value, (list,tuple)): return [public(v) for v in value]
    return value


def read(paths):
    rows, profiles, owners, metas = [], [], [], []
    keys = ['t','omega','g','U','mu','origin','tau_max','bins','radius','max_arcs','max_hop_pairs']
    for owner, path in enumerate(map(Path, paths)):
        m = json.loads((path/'run.json').read_text())
        if not m['complete']: raise ValueError(f'Incomplete chain: {path}')
        if metas and any(m[k] != metas[0][k] for k in keys): raise ValueError('Incompatible chains')
        b = np.loadtxt(path/'blocks.csv', delimiter=',', skiprows=1)
        p = np.loadtxt(path/'profiles.csv', delimiter=',', skiprows=1)
        if len(b) != m['blocks'] or len(p) != m['blocks']: raise ValueError('Incomplete blocks')
        if not np.isfinite(b).all() or not np.isfinite(p).all(): raise ValueError('Nonfinite data')
        rows.append(b[:,1:]); profiles.append(p[:,1:]); owners.extend([owner]*len(b)); metas.append(m)
    if len({m['seed'] for m in metas}) != len(metas): raise ValueError('Duplicate random seed')
    return np.concatenate(rows), np.concatenate(profiles), np.array(owners), metas


def reblock(rows, owners, factor):
    return np.concatenate([rows[owners==i].reshape(-1,factor,rows.shape[1]).sum(axis=1) for i in np.unique(owners)])


def reference_integral(m):
    rate=-m['mu']-(m['U'] if m['origin']==0 else 0)
    return -np.expm1(-rate*m['tau_max'])/rate


def normalize(total, m):
    bins=m['bins']; dt=m['tau_max']/bins
    if np.any(total[...,0] <= 0): raise ValueError('Insufficient zero-vertex normalization statistics')
    return total[...,2:2+bins]/total[...,0,None]*reference_integral(m)/dt


def jackknife(rows, function):
    total=rows.sum(axis=0)
    values=function(total[None,:]-rows)
    centered=values-values.mean(axis=0)
    covariance=(len(rows)-1)/len(rows)*centered.T@centered
    return function(total),covariance,values


def conservative_ratio(rows, owners, numerator, denominator):
    function=lambda total: numerator(total)/denominator(total)
    point=function(rows.sum(axis=0)); errors=[]
    for factor in [1,2,4,8]:
        if any(np.sum(owners==i)%factor for i in np.unique(owners)): continue
        grouped=reblock(rows,owners,factor)
        leave=function(grouped.sum(axis=0)[None,:]-grouped)
        errors.append(np.sqrt((len(grouped)-1)/len(grouped)*np.sum((leave-leave.mean(axis=0))**2,axis=0)))
    if len(np.unique(owners))>1:
        grouped=np.array([rows[owners==i].sum(axis=0) for i in np.unique(owners)])
        leave=function(grouped.sum(axis=0)[None,:]-grouped)
        errors.append(np.sqrt((len(grouped)-1)/len(grouped)*np.sum((leave-leave.mean(axis=0))**2,axis=0)))
    return point,np.max(errors,axis=0)


def green_estimate(rows, owners, m):
    grouped=reblock(rows,owners,4)
    return jackknife(grouped,lambda total:normalize(total,m))


def fit_window(rows, owners, m, lo, hi):
    bins=m['bins']; dt=m['tau_max']/bins; tau=(np.arange(bins)+.5)*dt
    mask=(tau>=lo)&(tau<hi)
    grouped=reblock(rows,owners,4)
    G,_,leave=green_estimate(rows,owners,m)
    if mask.sum()<4 or np.any(leave[:,mask]<=0): raise ValueError('Insufficient tail data')
    yj=np.log(leave[:,mask]); c=yj-yj.mean(axis=0)
    covariance=(len(grouped)-1)/len(grouped)*c.T@c
    values,vectors=np.linalg.eigh(covariance)
    keep=values>values[-1]*1e-7
    if keep.sum()<3: raise ValueError('Insufficient independent covariance modes')
    precision=(vectors[:,keep]/values[keep])@vectors[:,keep].T
    X=np.column_stack([np.ones(mask.sum()),tau[mask]])
    transform=np.linalg.solve(X.T@precision@X,X.T@precision)
    def parameters(total):
        beta=np.log(normalize(total,m)[...,mask])@transform.T
        E=m['mu']-beta[...,1]
        Z=np.exp(beta[...,0]-bin_log_factor(E-m['mu'],dt))
        return np.stack([E,Z],axis=-1)
    point=parameters(rows.sum(axis=0)); errors=[]
    for factor in [1,2,4,8]:
        if any(np.sum(owners==i)%factor for i in np.unique(owners)): continue
        group=reblock(rows,owners,factor); jk=parameters(group.sum(axis=0)[None,:]-group)
        errors.append(np.sqrt((len(group)-1)/len(group)*np.sum((jk-jk.mean(axis=0))**2,axis=0)))
    if len(np.unique(owners))>1:
        group=np.array([rows[owners==i].sum(axis=0) for i in np.unique(owners)])
        jk=parameters(group.sum(axis=0)[None,:]-group)
        errors.append(np.sqrt((len(group)-1)/len(group)*np.sum((jk-jk.mean(axis=0))**2,axis=0)))
    error=np.max(errors,axis=0)
    residual=np.log(G[mask])-X@(transform@np.log(G[mask]))
    result=dict(tau_min=lo,tau_max=hi,E=point[0],E_error=error[0],Z0=point[1],Z0_error=error[1],
        covariance_rank=int(keep.sum()),chi2_per_dof=float(residual@precision@residual/(keep.sum()-2)))
    for field,name in [(1,'direct_E'),(2,'Zb'),(3,'nph'),(4,'p0'),(5,'r2')]:
        value,err=conservative_ratio(rows,owners,
            lambda t,f=field:t[...,2+f*bins:2+(f+1)*bins][...,mask].sum(axis=-1),
            lambda t:t[...,2:2+bins][...,mask].sum(axis=-1))
        result[name],result[name+'_error']=value,err
    return public(result)


def summarize(paths, out, tail=True):
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    rows,profiles,owners,metas=read(paths);m=metas[0]
    G,C,jk=green_estimate(rows,owners,m)
    tau=(np.arange(m['bins'])+.5)*m['tau_max']/m['bins']
    np.savez_compressed(out/'green.npz',tau=tau,G=G,covariance=C,blocks=rows,owners=owners,profiles=profiles)
    np.savetxt(out/'green.csv',np.column_stack([tau,G,np.sqrt(np.diag(C))]),delimiter=',',header='tau,G_mu,error',comments='')
    result=dict(parameters={k:m[k] for k in ['t','omega','g','U','mu','origin','tau_max','bins','radius']},
        chains=len(metas),production_steps=sum(x['blocks']*x['steps_per_block'] for x in metas),
        reference_fraction=float(rows[:,0].sum()/rows[:,1].sum()),
        arc_cap_attempts=sum(x['arc_cap_attempts'] for x in metas),hop_cap_attempts=sum(x['hop_cap_attempts'] for x in metas),
        runs=[str(p) for p in paths],seeds=[x['seed'] for x in metas],
        maximum_arcs=max(x['maximum_arcs'] for x in metas),maximum_hops=max(x['maximum_hops'] for x in metas),
        maximum_extent=max(x['maximum_extent'] for x in metas),fit_windows=[],profiles=[])
    if tail:
        for w in range(3):
            result['fit_windows'].append(fit_window(rows,owners,m,(w+1)*m['tau_max']/4,(w+2)*m['tau_max']/4))
        result['primary']=result['fit_windows'][-1]
    sites=2*m['radius']+2
    for w in range(3):
        start=w*2*sites
        value,error=conservative_ratio(profiles,owners,lambda t:t[...,start:start+2*sites],
            lambda t:t[...,start:start+sites].sum(axis=-1)[...,None])
        result['profiles'].append(dict(tau_min=(w+1)*m['tau_max']/4,tau_max=(w+2)*m['tau_max']/4,
            x=list(range(-m['radius'],m['radius']+1)),p=value[:sites-1],p_error=error[:sites-1],
            Zi=value[sites:2*sites-1],Zi_error=error[sites:2*sites-1],
            p_overflow=value[sites-1],Zi_overflow=value[-1]))
    result['chain_estimates']=[fit_window(rows[owners==i],np.zeros(np.sum(owners==i),int),m,.75*m['tau_max'],m['tau_max']) for i in range(len(metas))] if tail else []
    (out/'summary.json').write_text(json.dumps(public(result),indent=2)+'\n')
    return public(result)


def exact_green(tau, m, kind):
    t,w,g,U,mu=(m[k] for k in ['t','omega','g','U','mu'])
    base=np.exp((mu+U)*tau)
    if kind=='atomic': return base*np.exp(g*g/w*tau-(g/w)**2*(-np.expm1(-w*tau)))
    if kind=='one_phonon': return base*(1+g*g*(tau/w+np.expm1(-w*tau)/w**2))
    if kind=='one_hop_pair':
        integral=tau*tau/2 if U==0 else tau/U+np.expm1(-U*tau)/U**2
        return base*(1+2*t*t*integral)
    if kind=='free': return np.exp((mu+2*t)*tau)*i0e(2*t*tau)
    if kind=='impurity':
        e=np.sqrt(4*t*t+U*U)
        continuum=quad(lambda q:(2*t*np.sin(q))**2/(np.pi*(U*U+(2*t*np.sin(q))**2))*np.exp((mu+2*t*np.cos(q))*tau),0,np.pi,epsabs=1e-12)[0]
        return U/e*np.exp((mu+e)*tau)+continuum
    raise ValueError(kind)


def validate(paths, out, kind):
    result=summarize(paths,out,tail=False)
    rows,_,owners,metas=read(paths);m=metas[0]
    G,C,_=green_estimate(rows,owners,m);dt=m['tau_max']/m['bins']
    exact=np.array([quad(lambda tau:exact_green(tau,m,kind),i*dt,(i+1)*dt,epsabs=1e-11)[0]/dt for i in range(m['bins'])])
    errors=[]
    for factor in [1,2,4,8]:
        _,cov,_=jackknife(reblock(rows,owners,factor),lambda t:normalize(t,m));errors.append(np.sqrt(np.diag(cov)))
    z=abs(G-exact)/np.maximum(np.max(errors,axis=0),1e-10*exact)
    check=dict(kind=kind,max_bin_standardized_error=float(z.max()),rms_standardized_error=float(np.sqrt(np.mean(z*z))),
               passed=bool(z.max()<5.5 and np.sqrt(np.mean(z*z))<2.5))
    if kind=='atomic':
        # The midpoint projector is not the vacuum-injection pole residue at
        # finite time. Test its exact coherent-state value independently.
        lo=.75*m['tau_max']; hi=m['tau_max']; mask=(np.arange(m['bins'])+.5)*dt>=lo
        denominator=quad(lambda t:exact_green(t,m,kind),lo,hi)[0]
        alpha=(m['g']/m['omega'])**2
        mean=lambda t:alpha*(-np.expm1(-m['omega']*t/2))**2
        check['projector_checks']={}
        for field,name,fn in [(2,'zero_ph',lambda t:np.exp(-mean(t))),(3,'nph',mean)]:
            expected=quad(lambda t:exact_green(t,m,kind)*fn(t),lo,hi)[0]/denominator
            value,error=conservative_ratio(rows,owners,
                lambda total,f=field:total[...,2+f*m['bins']:2+(f+1)*m['bins']][...,mask].sum(axis=-1),
                lambda total:total[...,2:2+m['bins']][...,mask].sum(axis=-1))
            passed=abs(value-expected)<max(5.5*error,1e-8)
            check['projector_checks'][name]=dict(value=float(value),error=float(error),expected=expected,passed=bool(passed))
            check['passed'] &= bool(passed)
    np.savez_compressed(Path(out)/'analytic_check.npz',G=G,exact=exact,error=np.max(errors,axis=0),z=z)
    (Path(out)/'validation.json').write_text(json.dumps(check,indent=2)+'\n')
    return check
