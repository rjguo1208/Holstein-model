"""Block-jackknife normalization and correlated long-time fits.

Histogram entries estimate BIN AVERAGES, not point values. The fit includes
the exact bin integral of the exponential, and propagates the shared n=0
normalization rather than treating bins as independent measurements.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import i0e


def read_runs(paths,filename='blocks.csv'):
    rows, owners, metadata = [], [], []
    keys = ('t','omega','g','k','mu','tau_max','max_order','bins')
    for i, path in enumerate(map(Path, paths)):
        meta = json.loads((path/'run.json').read_text())
        if not meta['complete']:
            raise ValueError(f'Incomplete run: {path}')
        data = np.genfromtxt(path/filename, delimiter=',', names=True)
        if len(data) != meta['blocks']:
            raise ValueError(f'Incomplete blocks: {path}')
        if metadata and any(meta[k] != metadata[0][k] for k in keys):
            raise ValueError('Combine only runs with identical physical/sampling parameters')
        rows.append(np.column_stack([data[n] for n in data.dtype.names[1:]]))
        owners.extend([i]*len(data)); metadata.append(meta)
    return np.concatenate(rows), np.asarray(owners), metadata


def reference_integral(meta):
    rate = -2*meta['t']*np.cos(meta['k'])-meta['mu']
    return meta['tau_max'] if abs(rate)<1e-12 else -np.expm1(-rate*meta['tau_max'])/rate


def combine_blocks(rows, owners, factor):
    groups=[]
    for i in np.unique(owners):
        a=rows[owners==i]
        if len(a)%factor:
            raise ValueError('Reblocking must divide each chain into equal blocks')
        groups.append(a.reshape(-1,factor,a.shape[1]).sum(axis=1))
    return np.concatenate(groups)


def normalize(totals, meta):
    n=meta['bins']; dt=meta['tau_max']/n
    return reference_integral(meta)*totals[...,3:3+n]/totals[...,0,None]/dt


def jackknife_error(values):
    b=len(values)
    return np.sqrt((b-1)/b*np.sum((values-values.mean(axis=0))**2,axis=0))


def estimate_green(rows, meta):
    total=rows.sum(axis=0)
    leave=normalize(total[None,:]-rows,meta)
    centered=leave-leave.mean(axis=0)
    cov=(len(rows)-1)/len(rows)*centered.T@centered
    return normalize(total,meta),cov,leave


def bin_log_factor(delta, dt):
    x=np.asarray(delta)*dt/2
    safe=np.where(np.abs(x)<1e-5,1.,x)
    return np.where(np.abs(x)<1e-5,x*x/6-x**4/180,np.log(np.sinh(safe)/safe))


def fit_tail(rows, owners, meta, lower=6., upper=20.):
    n=meta['bins']; dt=meta['tau_max']/n
    tau=(np.arange(n)+.5)*dt
    mask=(tau>=lower)&(tau<=upper)
    if mask.sum()<4: raise ValueError('Too few tail bins')
    green,cov,leave=estimate_green(rows,meta)
    if np.any(leave[:,mask]<=0): raise ValueError('Insufficient tail statistics')
    y=np.log(green[mask]); yj=np.log(leave[:,mask])
    centered=yj-yj.mean(axis=0)
    logcov=(len(rows)-1)/len(rows)*centered.T@centered
    vals,vec=np.linalg.eigh(logcov)
    floor=vals[-1]*1e-6
    precision=(vec/np.maximum(vals,floor))@vec.T
    X=np.column_stack([np.ones(mask.sum()),tau[mask]])
    estimator=np.linalg.solve(X.T@precision@X,X.T@precision)
    def parameters(values):
        beta=np.log(values[...,mask])@estimator.T
        E=meta['mu']-beta[...,1]
        Z=np.exp(beta[...,0]-bin_log_factor(E-meta['mu'],dt))
        return np.stack([E,Z],axis=-1)
    point=parameters(green)
    errors={}
    for factor in [1,2,4]:
        grouped=combine_blocks(rows,owners,factor)
        jack=normalize(grouped.sum(axis=0)[None,:]-grouped,meta)
        errors[str(factor)]=jackknife_error(parameters(jack)).tolist()
    error=np.max(np.asarray(list(errors.values())),axis=0)
    residual=y-X@(estimator@y)
    total=rows.sum(axis=0)
    energy=total[3+n:3+2*n][mask].sum()/total[3:3+n][mask].sum()
    grouped=combine_blocks(rows,owners,4)
    jack=grouped.sum(axis=0)[None,:]-grouped
    ej=jack[:,3+n:3+2*n][:,mask].sum(axis=1)/jack[:,3:3+n][:,mask].sum(axis=1)
    return dict(tau_min=lower,tau_max=upper,E=float(point[0]),Z=float(point[1]),
                E_error=float(error[0]),Z_error=float(error[1]),reblocking_errors=errors,
                chi2_per_dof=float(residual@precision@residual/(mask.sum()-2)),
                direct_energy=float(energy),direct_energy_error=float(jackknife_error(ej)),
                covariance_eigenvalue_floor=float(floor),bins=int(mask.sum()))


def exact_green(tau, meta, kind):
    t,omega,g,k,mu=(meta[n] for n in ('t','omega','g','k','mu'))
    eps=-2*t*np.cos(k)
    if kind=='free': return np.exp(-(eps-mu)*tau)
    if kind=='atomic':
        return np.exp((mu+g*g/omega)*tau-(g/omega)**2*(-np.expm1(-omega*tau)))
    if kind=='first_order':
        correction=quad(lambda l:(tau-l)*np.exp(-(omega-eps-2*t)*l)*i0e(2*t*l),0,tau,epsabs=1e-11)[0]
        return np.exp(-(eps-mu)*tau)*(1+g*g*correction)
    raise ValueError(kind)


def validate_run(path,kind,filename='blocks.csv'):
    rows,owners,metas=read_runs([path],filename); meta=metas[0]
    rows4=combine_blocks(rows,owners,4)
    green,cov,_=estimate_green(rows4,meta)
    dt=meta['tau_max']/meta['bins']; edges=np.arange(meta['bins']+1)*dt
    expected=np.asarray([quad(lambda tau:exact_green(tau,meta,kind),a,b,epsabs=1e-11)[0]/dt for a,b in zip(edges[:-1],edges[1:])])
    error=np.sqrt(np.diag(cov)); z=np.abs(green-expected)/np.maximum(error,1e-11*np.abs(expected))
    result=dict(kind=kind,max_bin_standardized_error=float(z.max()),
                rms_standardized_error=float(np.sqrt(np.mean(z*z))),
                crossing_measurements=int(rows[:,2].sum()),
                passed=bool(z.max()<5.5 and np.sqrt(np.mean(z*z))<2.0))
    if kind=='atomic' and filename=='blocks.csv':
        result['passed'] &= result['crossing_measurements']>0
    prefix='rb_' if filename=='rb_blocks.csv' else ''
    np.savetxt(Path(path)/(prefix+'validation_green.csv'),np.column_stack([(edges[:-1]+edges[1:])/2,green,error,expected]),delimiter=',',header='tau,G_mu,error,exact_G_mu',comments='')
    (Path(path)/(prefix+'validation.json')).write_text(json.dumps(result,indent=2)+'\n')
    return result


def summarize_conditional(paths,out):
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    rows,owners,metas=read_runs(paths,'rb_blocks.csv'); meta=metas[0]
    grouped=combine_blocks(rows,owners,4)
    green,cov,leave=estimate_green(grouped,meta)
    tau=(np.arange(meta['bins'])+.5)*meta['tau_max']/meta['bins']
    np.savez_compressed(out/'green.npz',tau=tau,G=green,covariance=cov,jackknife=leave,
                        blocks=grouped,owners=np.repeat(np.arange(len(metas)),len(grouped)//len(metas)))
    np.savetxt(out/'green.csv',np.column_stack([tau,green,np.sqrt(np.diag(cov))]),delimiter=',',header='tau,G_mu,error',comments='')
    result=dict(parameters={k:meta[k] for k in ('t','omega','g','k','mu','tau_max','max_order','bins')},
                method='conditional external-time integration',runs=[str(p) for p in paths],
                steps=sum(m['blocks']*m['steps_per_block'] for m in metas),
                cap_attempts=sum(m['cap_attempts'] for m in metas),
                maximum_order_observed=int(max(np.flatnonzero(np.loadtxt(Path(p)/'orders.csv',delimiter=',',skiprows=1)[:,1])[-1] for p in paths)))
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def summarize(paths,out):
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    rows,owners,metas=read_runs(paths); meta=metas[0]
    upper=min(20.,meta['tau_max']-1.)
    fits=[fit_tail(rows,owners,meta,start,upper) for start in [4.,6.,8.,10.] if start<upper-2]
    green,cov,leave=estimate_green(combine_blocks(rows,owners,4),meta)
    tau=(np.arange(meta['bins'])+.5)*meta['tau_max']/meta['bins']
    np.savez_compressed(out/'green.npz',tau=tau,G=green,covariance=cov,jackknife=leave)
    np.savetxt(out/'green.csv',np.column_stack([tau,green,np.sqrt(np.diag(cov))]),delimiter=',',header='tau,G_mu,error',comments='')
    histogram=sum(np.loadtxt(Path(p)/'orders.csv',delimiter=',',skiprows=1)[:,1] for p in paths)
    nz=np.flatnonzero(histogram)
    result=dict(parameters={k:meta[k] for k in ('t','omega','g','k','mu','tau_max','max_order','bins')},
                runs=[str(p) for p in paths],steps=sum(m['blocks']*m['steps_per_block'] for m in metas),
                seconds=sum(m['seconds'] for m in metas),primary=fits[1],fit_windows=fits,
                reference_fraction=float(rows[:,0].sum()/rows[:,1].sum()),
                maximum_order_observed=int(nz[-1]),cap_attempts=sum(m['cap_attempts'] for m in metas),
                crossing_measurements=int(rows[:,2].sum()),
                chain_estimates=[fit_tail(rows[owners==i],np.zeros(np.sum(owners==i),int),meta,6.,upper) for i in range(len(metas))])
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
