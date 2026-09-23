"""Audit frozen dense fits, compare independent VED, and draw vector maps."""
from __future__ import annotations
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from defect_analysis import public,jackknife,normalize
from defect_dense import prepare,estimate
from map_continuation import lorentz_map
from spectral_plot import vector_map,save_figure


def read_json(p):return json.loads(Path(p).read_text())
def digest(p):return sha256(Path(p).read_bytes()).hexdigest()

def relative_l1(a,b,e):return np.trapezoid(abs(a-b),e,axis=-1)/np.trapezoid(b,e,axis=-1)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for key,default in [('data','defect_dense_data_01'),('fit','defect_dense_fit_01'),('reference','defect_dense_reference_01'),('out','defect_dense_report_01')]:
        ap.add_argument('--'+key,type=Path,default=ROOT/'results'/default)
    args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=False)
    data,fit,ref=[read_json(p/'manifest.json') for p in [args.data,args.fit,args.reference]]
    assert all(m['complete'] for m in [data,fit,ref]);assert len(fit['points'])==34
    assert data['completed_chains']==272 and data['production_steps']==8704000000
    assert ref['frozen_manifest_sha256']==digest(args.fit/'manifest.json')
    for name,h in data['sources'].items():assert digest(args.data/'source'/name)==h
    assert digest(args.data/'source/defect_diagmc')==data['executable_sha256']
    for name,h in fit['source_hashes'].items():assert digest(args.fit/'source'/name)==h
    for r in ref['records']:assert digest(args.reference/r['file'])==r['sha256']
    axis=np.linspace(-4.5,5.5,1601);etas=np.array([.15,.25,.5,1.]);x=np.arange(-16,17)
    summary=dict(complete=False,production_steps=data['production_steps'],new_chains=272,
        measured_sites=list(range(17)),displayed_sites=x,energy=axis,eta=etas,raw_time_bins=192,
        bootstrap_per_point=64,bootstrap_total=2176,spectral_resolution_validated=False,
        optimization=read_json(ROOT/'results/defect_optimization_01/summary.json'),couplings={},plots=[],
        source_manifests={p.name:digest(p/'manifest.json') for p in [args.data,args.fit,args.reference]})
    old=read_json(ROOT/'results/defect_spectra_01/summary.json');allseeds=[];bundles=[]
    summary['baseline_input_sha256']={'results/defect_spectra_01/summary.json':digest(ROOT/'results/defect_spectra_01/summary.json')}
    summary['report_source_sha256']=digest(Path(__file__))
    summary['plot_helper_sha256']=digest(ROOT/'spectral_plot.py')
    for lam in [.25,.5]:
        central=[];low=[];high=[];candidates=[];reference=[[],[],[]];records=[];oldmaps=[]
        for site in range(17):
            label=f'lambda{lam:.2f}_i{site:02d}';pdir=args.fit/label;r=read_json(pdir/'summary.json')
            frozen=next(p for p in fit['points'] if p['coupling']==lam and p['site']==site)
            assert digest(pdir/'selection.json')==frozen['selected_sha256']
            assert r['complete'] and r['bootstrap_draws']==64 and not r['reference_used'] and r['order_caps']==0
            for name,h in r['input_sha256'].items():assert digest(args.data/'chains'/name)==h
            allseeds+=r['seeds']
            s=np.load(pdir/'spectra.npz');a=s['selected'];bs=s['bootstrap']
            assert a.shape==(4,1601) and bs.shape==(64,4,1601) and np.isfinite(bs).all() and np.all(bs>=0)
            np.testing.assert_allclose(s['bootstrap_quantiles'],np.quantile(bs,[.16,.5,.84],axis=0),rtol=0,atol=0)
            central.append(a);low.append(s['bootstrap_quantiles'][0]);high.append(s['bootstrap_quantiles'][2]);candidates.append(s['candidates'])
            notes=dict(site=site,configuration=r['config'],alpha=r['fit']['alpha'],cv=r['cv'][r['selected']]['score'],
                chi2_per_mode=r['fit']['chi2_per_mode'],covariance_rank=r['fit']['rank'],
                moment_residuals=r['fit']['moment_residuals'],upper_edge_weight=r['fit']['edge_weight'],
                configuration_names=[c['config']['name'] for c in r['cv']],
                configuration_scores=[c['score'] for c in r['cv']],
                configuration_relative_l1_vs_selected=relative_l1(s['candidates'],a,axis),
                standard_error_block16_over_block4=np.divide(s['blocking_errors'][4],s['blocking_errors'][2]),
                standard_error_chains_over_block4=np.divide(s['blocking_errors'][5],s['blocking_errors'][2]))
            if site<9:
                op=old[str(lam)]['records'][site];of=op['fit']
                oa=np.array([lorentz_map(axis,np.array(of['energies']),np.array(of['weights']),eta) for eta in etas]);oldmaps.append(oa)
                oldpath=ROOT/'results/defect_spectra_01/analysis'/label/'green.npz'
                summary['baseline_input_sha256'][str(oldpath.relative_to(ROOT))]=digest(oldpath)
                oldg=np.load(oldpath)
                newg=np.load(pdir/'green.npz');pp=r['parameters']
                _,b,_,bp=prepare(newg['blocks'],newg['owners'],pp,dict(factor=4,block=4));G,C=estimate(b,bp)
                ppold=dict(pp,bins=96)
                orows=oldg['blocks'][:,:98]
                _,ob,_,opb=prepare(orows,oldg['owners'],ppold,dict(factor=2,block=4));oG,oC=estimate(ob,opb)
                ratio=(np.sqrt(np.diag(C))/G)/(np.sqrt(np.diag(oC))/oG)
                notes['new_over_old_relative_green_error']=dict(median=float(np.median(ratio)),minimum=float(ratio.min()),maximum=float(ratio.max()))
            for j,(nh,radius) in enumerate([(10,40),(12,40),(12,64)]):
                v=np.load(args.reference/f'lambda{lam:.2f}_nh{nh}_R{radius}_i{site:02d}.npz')
                reference[j].append(v['A'])
            records.append(notes)
        A=np.stack(central,axis=1);V=np.stack(reference[2],axis=1)
        cloud=relative_l1(np.stack(reference[0],axis=1),np.stack(reference[1],axis=1),axis)
        space=relative_l1(np.stack(reference[1],axis=1),V,axis)
        errors=relative_l1(A,V,axis);old_error=relative_l1(np.stack(oldmaps,axis=1),V[:,:9],axis)
        lo=np.stack(low,axis=1);hi=np.stack(high,axis=1)
        np.savez_compressed(args.out/f'lambda{lam:.2f}_map.npz',x=x,measured_sites=np.arange(17),energy=axis,eta=etas,
            A=A[:,abs(x)],bootstrap_16=lo[:,abs(x)],bootstrap_84=hi[:,abs(x)],coupling=lam,U=1.,reflection_used=True,
            spectral_resolution_validated=False)
        np.savez_compressed(args.out/f'reference_lambda{lam:.2f}_map.npz',x=x,energy=axis,eta=etas,A=V[:,abs(x)],coupling=lam,U=1.)
        np.savez_compressed(args.out/f'lambda{lam:.2f}_diagnostics.npz',measured_sites=np.arange(17),energy=axis,eta=etas,
            A=A,A_ved=V,candidates=np.stack(candidates,axis=2),old_A=np.stack(oldmaps,axis=1),
            relative_l1=errors,old_relative_l1=old_error,reference_cloud_relative_l1=cloud,reference_spatial_relative_l1=space)
        summary['couplings'][str(lam)]=dict(records=records,relative_l1=errors,old_relative_l1=old_error,
            reference_cloud_relative_l1=cloud,reference_spatial_relative_l1=space,
            maximum_steps_relative_l1=np.max([r['steps_relative_l1'] for r in ref['records'] if r['coupling']==lam],axis=0),
            reference_maximum_moment_residual=max(max(abs(np.array(r['moment_residuals']))) for r in ref['records'] if r['coupling']==lam))
        bundles.append(dict(coupling=lam,A=A[:,abs(x)],A_ved=V[:,abs(x)]))
    assert len(set(allseeds))==272
    plt.rcParams.update({'font.size':10,'mathtext.fontset':'stix','svg.fonttype':'path','svg.hashsalt':'holstein-defect-dense-v1','savefig.facecolor':'white'})
    for eta in [.15,.25,1.]:
        ie=np.flatnonzero(etas==eta).item();vmax=max(b[key][ie].max() for b in bundles for key in ['A','A_ved']);norm=Normalize(0,vmax)
        fig,axes=plt.subplots(2,2,figsize=(10.4,8.9),layout='constrained',sharex=True,sharey=True);meshes=[]
        for row,key in enumerate(['A','A_ved']):
            for col,b in enumerate(bundles):
                ax=axes[row,col];mesh=vector_map(ax,x,axis,b[key][ie],'magma',norm);meshes.append(mesh)
                for edge in [-2,2]:ax.axhline(edge,color='white',ls='--',lw=.9,path_effects=[pe.Stroke(linewidth=1.5,foreground='black'),pe.Normal()])
                ax.axhline(-np.sqrt(5),color='#62d7ff',ls=(0,(2,2)),lw=.9,path_effects=[pe.Stroke(linewidth=1.5,foreground='black'),pe.Normal()])
                method='DiagMC reconstruction' if row==0 else 'Independent VED'
                ax.set(title=method+rf', $\lambda={b["coupling"]:g}$',xlim=(-16.5,16.5),ylim=(-4.5,5.5),xticks=[-16,-8,0,8,16])
                if row==1:ax.set_xlabel(r'Site $i$ (defect at $0$)')
                if col==0:ax.set_ylabel(r'Energy $\omega/t$ (vacuum = 0)')
        bar=fig.colorbar(mesh,ax=axes,label=r'$t A_{i,\eta}(\omega)$',shrink=.88,pad=.02);bar.solids.set_rasterized(False)
        fig.suptitle(rf'Single on-site defect: $t=\omega_0=U=1$, $\eta={eta:g}t$'+'\n'+
            r'White dashed: bare continuum edges $\pm2t$; blue dashed: $g=0$ bound state $-\sqrt{5}t$',fontsize=11)
        name=f'defect-dense-ldos-eta{eta:g}';save_figure(fig,args.out/name,meshes+[bar.solids]);plt.close(fig)
        summary['plots'].append(dict(name=name,eta=eta,vmax=vmax,vector_runs=[m.vector_audit for m in meshes],bare_edges=[-2.,2.],bare_bound=-np.sqrt(5)))
    fig,axes=plt.subplots(1,2,figsize=(10.4,4.1),layout='constrained',sharey=True)
    for ax,b in zip(axes,bundles):
        r=summary['couplings'][str(b['coupling'])]
        for ie,color in [(1,'#245b87'),(3,'#9b3d20')]:
            ax.plot(range(17),np.array(r['relative_l1'])[ie]*100,'o-',ms=3,color=color,label=rf'New, $\eta={etas[ie]:g}t$')
            ax.plot(range(9),np.array(r['old_relative_l1'])[ie]*100,'x--',color=color,label=rf'Previous, $\eta={etas[ie]:g}t$')
        ax.set(title=rf'$\lambda={b["coupling"]:g}$',xlabel=r'Measured site $i$',xticks=[0,4,8,12,16]);ax.legend(fontsize=8)
    axes[0].set_ylabel('Full-spectrum relative L1 difference (%)');save_figure(fig,args.out/'defect-dense-error');plt.close(fig)
    summary['complete']=True;(args.out/'summary.json').write_text(json.dumps(public(summary),indent=2)+'\n')
    shutil.copyfile(Path(__file__),args.out/'report_source.py')
    print(json.dumps({lam:dict(error_range_percent=(100*np.c_[np.min(s['relative_l1'],axis=1),np.max(s['relative_l1'],axis=1)]).tolist(),
        reference_cloud_max=np.max(s['reference_cloud_relative_l1'],axis=1).tolist()) for lam,s in summary['couplings'].items()},indent=2))


if __name__=='__main__':main()
