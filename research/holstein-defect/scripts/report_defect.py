"""Audit and plot the completed single-defect calculations without changing data."""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
from hashlib import sha256
import json
from pathlib import Path
import re
import shutil
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from defect_analysis import public,read,reblock,conservative_ratio


def audit(directory):
    manifest=json.loads((directory/'manifest.json').read_text())
    assert manifest['complete']
    for name,digest in manifest['sources'].items():
        assert sha256((directory/'source'/name).read_bytes()).hexdigest()==digest
    assert sha256((directory/'source/defect_diagmc').read_bytes()).hexdigest()==manifest['executable_sha256']
    for name,digest in manifest.get('inputs',{}).items():
        assert sha256((directory/name).read_bytes()).hexdigest()==digest
    seeds=[];steps=0;caps=0
    for task in manifest['tasks']:
        path=directory/'chains'/task['name'];r=json.loads((path/'run.json').read_text())
        assert r['complete']
        for name,value in task['parameters'].items():
            assert r[name.replace('-','_')]==value
        for name in ['blocks.csv','profiles.csv']:
            assert (path/name).read_bytes().count(b'\n')==r['blocks']+1
        seeds.append(r['seed']);steps+=r['blocks']*r['steps_per_block']
        if not task['name'].startswith('check_'):caps+=r['arc_cap_attempts']+r['hop_cap_attempts']
    assert len(set(seeds))==len(seeds) and caps==0
    return dict(chains=len(seeds),production_steps=steps,seeds=seeds,
        raw_blocks_verified=True,parameters_verified=True,source_hashes_verified=True,
        executable_hash_verified=True,input_hashes_verified=True,production_cap_attempts=caps)


def localization(case,directory):
    # Resolve archived chain names against the selected scan, not the original
    # cluster's absolute paths recorded in the immutable sampling summary.
    paths=[directory/'chains'/Path(p).name for p in case['diagmc']['runs']]
    raw,profiles,owners,metas=read(paths)
    m=metas[0];sites=2*m['radius']+2;start=4*sites
    data=profiles[:,start:start+sites]
    x=np.arange(-m['radius'],m['radius']+1)
    low,high=(3,9) if m['U']<1 else (2,6) if m['U']<2 else (1,4)
    result=[]
    # Several fixed tail windows expose dependence on the fitted spatial range.
    for lo,hi in [(low,high),(low+1,high+1),(low+2,high+2)]:
        radii=np.arange(lo,hi+1)
        indices_plus=radii+m['radius'];indices_minus=-radii+m['radius']
        def probability(total):
            return (total[...,indices_plus]+total[...,indices_minus])/(2*total.sum(axis=-1)[...,None])
        value,error=conservative_ratio(data,owners,lambda t:(t[...,indices_plus]+t[...,indices_minus])/2,
                                      lambda t:t.sum(axis=-1)[...,None])
        keep=(value>3*error)&(value>0)
        if keep.sum()<3:continue
        radii=radii[keep];X=np.column_stack([np.ones(len(radii)),radii])
        weights=(value[keep]/error[keep])**2
        transform=np.linalg.solve(X.T@(weights[:,None]*X),X.T*weights)
        def xi(total):
            p=probability(total)[...,keep]
            if np.any(p<=0):raise ValueError('Too few tail counts for logarithmic fitting')
            slope=(np.log(p)@transform.T)[...,1]
            return -2/slope
        point=xi(data.sum(axis=0));errors=[]
        for factor in [1,2,4,8]:
            grouped=reblock(data,owners,factor);j=xi(grouped.sum(axis=0)[None,:]-grouped)
            errors.append(float(np.sqrt((len(grouped)-1)/len(grouped)*np.sum((j-j.mean())**2))))
        grouped=np.array([data[owners==i].sum(axis=0) for i in np.unique(owners)])
        j=xi(grouped.sum(axis=0)[None,:]-grouped)
        errors.append(float(np.sqrt((len(grouped)-1)/len(grouped)*np.sum((j-j.mean())**2))))
        ref=case['reference'];refx=np.array(ref['x']);refp=np.array(ref['p'])
        p=np.array([(refp[refx==r][0]+refp[refx==-r][0])/2 for r in radii])
        refxi=-2/(np.log(p)@transform.T)[1]
        result.append(dict(radii=radii.tolist(),xi=float(point),error=max(errors),reference=float(refxi)))
    return dict(fits=result,meaning='finite-window exponential fit p_i proportional to exp(-2|i|/xi); not an assumed exact envelope')


def savefig(fig,path,meshes=()):
    # The stem may contain a decimal broadening (eta0.25).
    svg=Path(str(path)+'.svg')
    fig.savefig(svg)
    svg.write_text(re.sub(r'(<g id="QuadMesh_[^"]+")>',r'\1 shape-rendering="crispEdges">',svg.read_text()))
    for mesh in meshes:mesh.set_edgecolor('face');mesh.set_linewidth(.25)
    fig.savefig(Path(str(path)+'.pdf'))
    for mesh in meshes:mesh.set_edgecolor('none');mesh.set_linewidth(0)
    fig.savefig(Path(str(path)+'.png'),dpi=180)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scan',type=Path,default=ROOT/'results/defect_scan_01')
    parser.add_argument('--spectra',type=Path,default=ROOT/'results/defect_spectra_01')
    parser.add_argument('--out',type=Path,default=ROOT/'results/defect_report_01')
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=False)
    ground=json.loads((args.scan/'summary.json').read_text())
    spectral=json.loads((args.spectra/'summary.json').read_text())
    assert ground['validation_passed'] and ground['first_point_passed']
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),audit={},cases={},spectra={})
    allseeds=[]
    for directory in [ROOT/'results/defect_pilot_01',args.scan,args.spectra]:
        a=audit(directory);allseeds+=a['seeds'];report['audit'][directory.name]=a
    assert len(allseeds)==len(set(allseeds))
    report['total_chains']=len(allseeds)
    report['production_steps']=sum(a['production_steps'] for a in report['audit'].values())
    clean={}
    for lam in [.25,.5]:
        p=ROOT/f'results/mass_direct_02/ground/lambda{lam:.2f}/summary.json'
        # Locate the saved merged clean-system ground fit without inserting a
        # VED energy into the reported DiagMC binding energy.
        if not p.exists():
            candidates=list((ROOT/'results/mass_direct_02/ground').glob(f'*{lam:.2f}*/summary.json'))
            if len(candidates)!=1:raise ValueError('Need the existing clean DiagMC ground reference')
            p=candidates[0]
        clean[str(lam)]=dict(json.loads(p.read_text())['primary'],source=str(p),sha256=sha256(p.read_bytes()).hexdigest())
    report['clean_diagmc']=clean
    for name,case in ground['cases'].items():
        p=case['diagmc']['primary'];m=case['diagmc']['parameters'];lam=m['g']**2/(2*m['t']*m['omega'])
        lam=min([.25,.5],key=lambda v:abs(v-lam));bulk=clean[str(lam)]
        profile=case['diagmc']['profiles'][-1];mid=profile['x'].index(0)
        assert abs(sum(profile['p'])+profile['p_overflow']-1)<1e-10
        assert abs(sum(profile['Zi'])+profile['Zi_overflow']-p['Zb'])<1e-10
        report['cases'][name]=dict(coupling=lam,U=m['U'],E=p['E'],E_error=p['E_error'],
            binding=bulk['E']-p['E'],binding_error=float(np.hypot(bulk['E_error'],p['E_error'])),
            p0=p['p0'],p0_error=p['p0_error'],Zb=p['Zb'],Zb_error=p['Zb_error'],
            Z0_projected=profile['Zi'][mid],Z0_projected_error=profile['Zi_error'][mid],
            rms_radius=float(np.sqrt(p['r2'])),rms_radius_error=float(p['r2_error']/(2*np.sqrt(p['r2']))),
            nph=p['nph'],nph_error=p['nph_error'],localization=localization(case,args.scan),
            reference={k:case['reference'][k] for k in ['E','E_clean','binding','Zb','p0','Z0','r2','nph','convergence']},
            projection_window_changes=case['projection_window_changes'],profile_overflow=profile['p_overflow'],
            chain_estimates=case['diagmc']['chain_estimates'],fit_windows=case['diagmc']['fit_windows'])
    plt.rcParams.update({'font.size':11,'mathtext.fontset':'stix','svg.fonttype':'path','savefig.facecolor':'white'})
    fig,axes=plt.subplots(1,2,figsize=(10.2,4.5),layout='constrained')
    for lam,color in [(.25,'#245b87'),(.5,'#9b3d20')]:
        cases=sorted([r for r in report['cases'].values() if r['coupling']==lam],key=lambda r:r['U'])
        U=np.array([r['U'] for r in cases])
        for ax,key,err in [(axes[0],'binding','binding_error'),(axes[1],'Zb','Zb_error')]:
            ax.plot(U,[r['reference'][key] for r in cases],color=color,label=rf'VED $\lambda={lam:g}$')
            ax.errorbar(U,[r[key] for r in cases],yerr=[r[err] for r in cases],fmt='o',color=color,capsize=3,label=rf'DiagMC $\lambda={lam:g}$')
            ax.set_xlabel(r'$U/t$')
    axes[0].set_ylabel(r'Binding energy $\Delta_b/t$');axes[1].set_ylabel(r'Zero-phonon weight $Z_b$')
    axes[0].legend(fontsize=9);axes[1].legend(fontsize=9)
    savefig(fig,args.out/'defect-binding')
    fig,axes=plt.subplots(1,2,figsize=(10.2,4.5),layout='constrained',sharey=True)
    for ax,lam in zip(axes,[.25,.5]):
        case=ground['cases'][f'lambda{lam:.2f}_U1'];p=case['diagmc']['profiles'][-1];ref=case['reference']
        for key,color,label in [('p','#245b87',r'$p_i$'),('Zi','#9b3d20',r'$Z_i$')]:
            ax.plot(ref['x'],ref[key],color=color,label='VED '+label)
            ax.errorbar(p['x'],p[key],yerr=p[key+'_error'],fmt='o',ms=3,color=color,capsize=2,label='DiagMC '+label)
        ax.set(xlim=(-8,8),xlabel=r'Site $i$',title=rf'$U/t=1,\;\lambda={lam:g}$');ax.legend(fontsize=9)
    axes[0].set_ylabel('Probability / pole weight')
    savefig(fig,args.out/'defect-density')
    for eta in [.25,1.]:
        fig,axes=plt.subplots(2,2,figsize=(10.2,8.4),layout='constrained',sharex=True,sharey=True)
        entries=[]
        for row,prefix in enumerate(['','reference_']):
            for col,lam in enumerate([.25,.5]):
                with np.load(args.spectra/f'{prefix}lambda{lam:.2f}_map.npz') as d:
                    ie=int(np.flatnonzero(np.isclose(d['eta'],eta))[0]);A=d['A'][ie];x=d['x'];energy=d['energy']
                    assert A.shape==(17,801) and np.isfinite(A).all() and np.all(A>=0)
                    entries.append((row,col,lam,A,x,energy))
        vmax=np.ceil(max(e[3].max() for e in entries)*10)/10;meshes=[]
        for row,col,lam,A,x,energy in entries:
            ax=axes[row,col];mesh=ax.pcolormesh(x,energy,A.T,shading='nearest',cmap='magma',vmin=0,vmax=vmax,rasterized=False,antialiased=False,linewidth=0)
            meshes.append(mesh);method='DiagMC (preliminary)' if row==0 else 'VED reference'
            ax.set(title=method+rf', $\lambda={lam:g}$',xlim=(-8.5,8.5),ylim=(-4.5,5.5),xticks=[-8,-4,0,4,8])
            if row==1:ax.set_xlabel(r'Site $i$ (defect at $0$)')
            if col==0:ax.set_ylabel(r'Energy $\omega/t$ (vacuum = 0)')
        cb=fig.colorbar(mesh,ax=axes.ravel().tolist(),fraction=.035,pad=.02);cb.set_label(r'$t A_{i,\eta}(\omega)$');cb.solids.set_rasterized(False);meshes.append(cb.solids)
        fig.suptitle(rf'Single on-site defect: $t=\omega_0=U=1,\;\eta={eta:g}t$')
        savefig(fig,args.out/f'defect-ldos-eta{eta:g}',meshes)
    for lam in ['0.25','0.5']:
        s=spectral[lam]
        report['spectra'][lam]=dict(comparison=s['comparison'],
            reference_checks_passed=s['reference_checks_passed'],
            reference_max_cloud_change=float(np.max(s['reference_cloud_relative_l1'])),
            reference_max_spatial_change=float(np.max(s['reference_spatial_relative_l1'])),
            reference_max_step_change=max(c['step_error'] for c in s['reference_steps']),
            maximum_moment_residual=max(max(abs(np.array(r['fit']['moment_residuals']))) for r in s['records']),
            maximum_cv_score=max(float(np.mean(r['cv']['folds'],axis=0)[r['cv']['selected']]) for r in s['records']),
            spectral_resolution_validated=False)
    (args.out/'summary.json').write_text(json.dumps(public(report),indent=2)+'\n')
    shutil.copyfile(Path(__file__),args.out/'report_source.py')
    print(json.dumps({k:v for k,v in report.items() if k not in ['audit','cases','spectra','clean_diagmc']},indent=2))


if __name__=='__main__':main()
