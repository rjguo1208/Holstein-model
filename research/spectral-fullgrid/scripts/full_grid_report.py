"""Audit all 41 frozen DiagMC spectra, then compare with independent VED.

The 10 completed pilot momenta are reused without changes. The same frozen
selection protocol is applied to the remaining 31; VED never selects fits.
The representative new momenta (indices 8, 20, 32) were chosen before these
new VED comparisons. Existing six-point synthetic tests are reused, not
extrapolated into a full-grid resolution certificate.
"""
from collections import Counter
from datetime import datetime,timezone
from hashlib import sha256
import json
from pathlib import Path
import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from continuation import bin_kernel
from map_continuation import lorentz_map
from spectral_cv import whitening
from pilot_diagnostics import public
from pilot_resolution import peak_test

def save(fig,out,name):
    for extension in ['svg','pdf','png']:
        fig.savefig(out/f'{name}.{extension}',dpi=170,bbox_inches='tight')
    plt.close(fig)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=ROOT/'results/report')
    args=parser.parse_args()
    plan=json.loads((ROOT/'plans/pilot.json').read_text())
    assert plan['indices']==list(range(41))
    pipeline=json.loads((ROOT/'results/pipeline_status.json').read_text())
    assert pipeline['complete'] and pipeline['points']==82
    reused=json.loads((ROOT/'results/reuse_manifest.json').read_text())
    assert reused['complete'] and reused['chains']==160
    for name,expected in reused['sha256'].items():
        assert sha256((ROOT/name).read_bytes()).hexdigest()==expected, name
    frozen_source=json.loads((ROOT/'results/inference_source.json').read_text())
    for name,expected in frozen_source['sha256'].items():
        for path in [ROOT/name,ROOT/'results/inference_source'/name]:
            assert sha256(path.read_bytes()).hexdigest()==expected, str(path)
    checkpoints=[]
    for path in sorted((ROOT/'results/checkpoints').glob('job*/retained_outputs.json')):
        checkpoint=json.loads(path.read_text())
        for name,expected in checkpoint['sha256'].items():
            assert sha256((ROOT/name).read_bytes()).hexdigest()==expected,name
        for name,expected in checkpoint.get('saved_failed_batches_sha256',{}).items():
            assert sha256((ROOT/name).read_bytes()).hexdigest()==expected,name
            recovered_initial=ROOT/'results/bootstrap_initial'/Path(name).name
            assert sha256(recovered_initial.read_bytes()).hexdigest()==expected,str(recovered_initial)
        checkpoints.append(dict(stopped_job_id=checkpoint['stopped_job_id'],
            retained_files_verified=len(checkpoint['sha256']),
            retained_bootstrap_points=len(checkpoint['completed_bootstrap_points']),
            saved_failed_batch_files_verified=len(checkpoint.get('saved_failed_batches_sha256',{}))))
    freeze=json.loads((ROOT/'results/evaluation_freeze.json').read_text())
    for name,digest in freeze['comparison_hashes'].items():
        assert sha256((ROOT/'results/comparison'/name).read_bytes()).hexdigest()==digest
    assert json.loads((ROOT/'results/comparison/provenance.json').read_text())['complete']
    boot=json.loads((ROOT/'results/bootstrap/summary.json').read_text());assert boot['complete']
    recovery=json.loads((ROOT/'results/bootstrap_recovery.json').read_text())
    assert recovery['complete'] and recovery['unchanged_central_selections']
    assert recovery['unchanged_successful_bootstrap_draws'] and not recovery['reference_used']
    for name,expected in recovery['source_hashes'].items():
        assert sha256((ROOT/name).read_bytes()).hexdigest()==expected,name
    numerical_source=json.loads((ROOT/'results/numerical_recovery-source.json').read_text())
    assert numerical_source['sha256']==recovery['source_hashes']
    prior_solver_snapshots=[]
    for path in sorted((ROOT/'results/numerical_recovery').glob('source-job*/source-manifest.json')):
        earlier=json.loads(path.read_text())
        for name,expected in earlier['sha256'].items():
            assert sha256((path.parent/Path(name).name).read_bytes()).hexdigest()==expected,name
        prior_solver_snapshots.append(dict(path=str(path.relative_to(ROOT)),files_verified=len(earlier['sha256'])))
    assert len(recovery['analytic_validation'])==3 and recovery['independent_controls']
    assert all(p['maximum_weight_error']<2e-7 for p in recovery['analytic_validation'])
    assert all(p['maximum_broadened_curve_difference']<1e-5 for p in recovery['independent_controls'])
    assert all(p['gradient']<1e-8 and p['moment_error']<1e-6 for p in recovery['stable_solver_recoveries'])
    resolution=json.loads((ROOT/'results/resolution/summary.json').read_text());assert resolution['complete']
    raw=json.loads((ROOT/'results/raw_audit.json').read_text())
    diag=json.loads((ROOT/'results/diagnostics/summary.json').read_text())
    assert raw['chains']==raw['complete_chains']==656
    assert raw['new_chains']==496 and raw['reused_chains']==160
    assert len(raw['raw_sha256'])==656*4
    for name,expected in raw['raw_sha256'].items():
        assert sha256((ROOT/name).read_bytes()).hexdigest()==expected, name
    assert len(boot['points'])==82
    assert all(p['requested']==p['successful']==128 and not p['failures'] for p in boot['points'])
    assert len(diag['points'])==164 and len(diag['comparisons'])==82
    out=args.out;out.mkdir(parents=True,exist_ok=False)
    points=[];bundles={};refhash={};resolution_cases=[]
    for lam in plan['couplings']:
        baseline_dir=Path(plan['baseline_spectra'])
        if (ROOT/'inputs/baseline_maps').exists():baseline_dir=ROOT/'inputs/baseline_maps'
        baseline=np.load(baseline_dir/f'lambda{lam:.2f}_map.npz')
        rp=baseline_dir.parent/'reference_map_01'/f'lambda{lam:.2f}_map.npz'
        reference=np.load(rp);refhash[rp.name]=sha256(rp.read_bytes()).hexdigest()
        axis=baseline['energy'];etas=baseline['eta']
        np.testing.assert_array_equal(axis,reference['energy'])
        bundle={k:[] for k in ['baseline','sampling_only','selected','reference','quantiles','sensitivity_low','sensitivity_high']}
        for ik in plan['indices']:
            stem=f'lambda{lam:.2f}_ik{ik:03d}'
            path=ROOT/'results/comparison'/f'{stem}.json'
            record=json.loads(path.read_text())
            for source,digest in record['input_hashes'].items():
                source_path=ROOT/'results/pilot_data'/source.split('/pilot_data/',1)[1]
                if not source_path.exists():source_path=Path(source)
                assert sha256(source_path.read_bytes()).hexdigest()==digest
            d=np.load(ROOT/'results/comparison'/f'{stem}.npz')
            b=np.load(ROOT/f'results/bootstrap/lambda{lam:.2f}_ik{ik:03d}.npz')
            bm=json.loads((ROOT/f'results/bootstrap/{stem}.json').read_text())
            assert bm['comparison_sha256']==sha256(path.read_bytes()).hexdigest()
            assert b['samples'].shape==(bm['successful'],3,701)
            assert bm['successful']==128 and not bm['failures']
            for key in ['primary','sampling_only','candidates','sensitivity_low','sensitivity_high']:
                assert np.isfinite(d[key]).all() and (d[key]>=0).all()
            assert np.isfinite(b['samples']).all() and (b['samples']>=0).all()
            assert np.all(np.trapezoid(b['samples'],axis,axis=2)<1+1e-5)
            mr=float(np.max(abs(np.asarray(record['primary_fit']['moment_residuals']))))
            assert mr<1e-5
            refs=np.array([reference['A'][int(np.flatnonzero(np.isclose(reference['eta'],eta))[0]),ik] for eta in etas])
            curves=np.array([baseline['A'][:,ik],d['sampling_only'],d['primary']])
            errors=np.trapezoid(abs(curves-refs),axis,axis=2)/np.trapezoid(refs,axis,axis=1)
            diagnostic=next(p for p in diag['points'] if (p['dataset'],p['coupling'],p['index'])==('baseline',lam,ik))
            tt=np.asarray(diagnostic['tau']);G=np.asarray(diagnostic['G']);C=np.asarray(diagnostic['covariance']);use=tt<=8
            W=whitening(G[use],C[np.ix_(use,use)],1e-8)
            external=[]
            base=next(r for r in record['candidates'] if r['name']=='baseline')
            for key,fit in [('sampling_only',base['cv']['minimum_fit']),('primary',record['primary_fit'])]:
                regenerated=np.array([lorentz_map(axis,np.array(fit['energies']),np.array(fit['weights']),eta) for eta in etas])
                np.testing.assert_allclose(d[key],regenerated,atol=1e-12,rtol=1e-12)
                pred=bin_kernel(tt[use],.25,fit['energies'],-2.8)@np.array(fit['weights'])
                residual=W@(pred-G[use]);external.append(float(residual@residual/len(W)))
            point=dict(coupling=lam,index=ik,k=baseline['k'][ik],configuration=record['primary'],
                alpha=record['primary_fit']['alpha'],moment_residual=mr,
                errors=dict(baseline=errors[0],sampling_only=errors[1],selected=errors[2]),
                selected_cv_score=record['primary_cv']['score'],
                independent_old_data_scores=dict(sampling_only=external[0],selected=external[1]),
                bootstrap_successful=len(b['samples']),eligible_candidates=record['eligible_candidates'])
            points.append(point)
            for key,value in [('baseline',curves[0]),('sampling_only',curves[1]),('selected',curves[2]),('reference',refs),
                ('quantiles',b['quantiles']),('sensitivity_low',d['sensitivity_low']),('sensitivity_high',d['sensitivity_high'])]:bundle[key].append(value)
        arrays={key:np.stack(value,axis=1 if key!='quantiles' else 2) for key,value in bundle.items()}
        arrays.update(k=baseline['k'][plan['indices']],indices=plan['indices'],energy=axis,eta=etas,
            bootstrap_probabilities=[.025,.16,.84,.975],coupling=lam,
            spectral_resolution_validated=False,
            method='all 41 measured momenta: independent 4x DiagMC; frozen chain-CV selection; 128 block bootstrap draws')
        np.savez_compressed(out/f'lambda{lam:.2f}_fullgrid.npz',**arrays)
        bundles[lam]=arrays
    for path in sorted((ROOT/'results/resolution').glob('lambda*.json')):
        r=json.loads(path.read_text());cases=[]
        truth_arrays=np.load(path.with_suffix('.npz'))
        for ie,eta in enumerate(r['eta']):
            errors=[p['relative_l1'][ie] for p in r['records']]
            recovered=sum(p['peaks'][ie].get('doublet_recovered',False) for p in r['records'])
            false_split=sum(p['peaks'][ie].get('false_split',False) for p in r['records'])
            cases.append(dict(eta=eta,median_relative_l1=float(np.median(errors)) if errors else None,
                max_relative_l1=max(errors) if errors else None,doublet_recovered=recovered,
                false_splits=false_split,denominator=r['replicates'],
                truth_doublet_visible=peak_test(truth_arrays['energy'],truth_arrays['truth'][ie],r['truth']).get('doublet_recovered',False)))
        resolution_cases.append(dict(coupling=r['coupling'],index=r['index'],case=r['case'],
            successful=r['successful'],requested=r['replicates'],results=cases))
    summary={}
    for lam in plan['couplings']:
        selected=[p for p in points if p['coupling']==lam];stats={}
        for key in ['baseline','sampling_only','selected']:
            values=np.array([p['errors'][key] for p in selected])
            stats[key]=dict(median=np.median(values,axis=0),maximum=np.max(values,axis=0),minimum=np.min(values,axis=0))
        ratios=[np.median(np.array(p['sigma_ratio_pilot_over_baseline'])[:32]) for p in diag['comparisons'] if p['coupling']==lam]
        reblock=[p['reblocking'][-1]['median_sigma_over_factor4'] for p in diag['points'] if p['dataset']=='pilot' and p['coupling']==lam]
        stats.update(improved_points_at_eta025=sum(p['errors']['selected'][0]<p['errors']['baseline'][0] for p in selected),
            median_sigma_ratio=float(np.median(ratios)),point_sigma_ratio_range=[min(ratios),max(ratios)],
            block16_over4_median_sigma_range=[min(reblock),max(reblock)],
            maximum_selected_cv=max(p['selected_cv_score'] for p in selected),
            independent_old_data_scores={key:dict(median=float(np.median([p['independent_old_data_scores'][key] for p in selected])),
               maximum=max(p['independent_old_data_scores'][key] for p in selected)) for key in ['sampling_only','selected']})
        summary[str(lam)]=stats
    report=dict(complete=True,generated_utc=datetime.now(timezone.utc).isoformat(),indices=plan['indices'],
        couplings=plan['couplings'],eta=[.25,.5,1.],energy_window=[-3.25,5.5],
        metric='whole-spectrum relative L1 on identical energy/k/eta arrays; all 41 momenta per coupling',
        new_chains=raw['new_chains'],reused_chains=raw['reused_chains'],total_chains=raw['chains'],
        new_production_steps=raw['new_production_steps'],
        reused_production_steps=reused['production_steps'],
        total_production_steps=raw['production_steps'],sampling_ratio=4,
        new_indices=plan['new_indices'],reused_indices=plan['reused_indices'],
        reused_files_verified=len(reused['sha256']),raw_files_verified=len(raw['raw_sha256']),
        frozen_source_files_verified=len(frozen_source['sha256']),
        resumed_job_checkpoints=checkpoints,
        representative_new_indices=[8,20,32],
        bootstrap_replicates_per_point=128,bootstrap_points=82,bootstrap_draws=82*128,
        bootstrap_numerical_recoveries=sum(len(p.get('initial_failures',[])) for p in boot['points']),bootstrap_failures=sum(len(p['failures']) for p in boot['points']),
        bootstrap_stable_solver_recoveries=len(recovery['stable_solver_recoveries']),
        numerical_recovery_validation=recovery,
        previous_numerical_solver_source_snapshots=prior_solver_snapshots,
        synthetic_cases=len(resolution_cases),synthetic_draws=sum(p['requested'] for p in resolution_cases),
        synthetic_scope='Reused pilot tests: 2 couplings x indices 25,36,39; no new synthetic draws',
        synthetic_indices=[25,36,39],new_synthetic_draws=0,
        synthetic_failures=sum(p['requested']-p['successful'] for p in resolution_cases),
        selected_configurations=dict(Counter(p['configuration'] for p in points)),
        reference_used_in_selection=False,selection_freeze=freeze,reference_hashes=refhash,
        source_manifest_sha256=sha256((ROOT/'results/inference_source.json').read_bytes()).hexdigest(),
        reuse_manifest_sha256=sha256((ROOT/'results/reuse_manifest.json').read_bytes()).hexdigest(),
        report_script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        maximum_moment_residual=max(p['moment_residual'] for p in points),
        spectral_resolution_validated=False,summary=summary,points=points,resolution=resolution_cases,
        limitations=['All 41 measured momenta use the same refinement protocol; momentum-grid convergence and narrow-peak accuracy are not thereby certified.',
          'Bootstrap reselects alpha within the frozen representation; parameterization bias is shown separately.',
          'Reused synthetic recovery rates cover only six pilot coupling/momentum points, conditional on chosen spectra and measured covariance, with 24 draws per case.',
          'One-standard-error selection is a heuristic; shared training folds do not give independent confidence tests.',
          'Pointwise spectral error need not decrease after refinement; VED comparisons do not tune or replace selected fits.'])
    (out/'summary.json').write_text(json.dumps(public(report),indent=2,allow_nan=False)+'\n')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    fig,axes=plt.subplots(1,2,figsize=(11,4.3),constrained_layout=True)
    for ax,lam in zip(axes,plan['couplings']):
        ps=[p for p in points if p['coupling']==lam];ks=np.array([p['k'] for p in ps])/np.pi
        for key,label,color,marker in [('baseline','Original','#727b87','o'),('sampling_only','4x sampling, original selection','#2580a0','s'),('selected','4x sampling + refinement','#cb6238','D')]:
            ax.plot(ks,[100*p['errors'][key][0] for p in ps],marker=marker,color=color,label=label,lw=1,ms=4)
        ax.set(xlabel=r'$k/\pi$',ylabel='Relative spectral L1 error (%)',title=rf'$\lambda={lam}$, $\eta=0.25t$')
        ax.grid(alpha=.2);ax.set_ylim(bottom=0)
    axes[0].legend(fontsize=8);fig.suptitle('Independent VED comparison at all 41 measured momenta')
    save(fig,out,'full-grid-error-comparison')
    fig,axes=plt.subplots(2,3,figsize=(12,7.2),constrained_layout=True)
    for row,lam in enumerate(plan['couplings']):
        b=bundles[lam]
        for col,ik in enumerate([8,20,32]):
            j=plan['indices'].index(ik);ax=axes[row,col]
            ax.fill_between(b['energy'],b['quantiles'][1,0,j],b['quantiles'][2,0,j],color='#cb6238',alpha=.22,label='128-draw 16-84% range')
            for key,label,color,style in [('baseline','Original','#8a8f96','-'),('sampling_only','4x only','#2580a0',':'),('selected','Refined','#cb6238','-'),('reference','VED reference','#171717','--')]:
                ax.plot(b['energy'],b[key][0,j],color=color,ls=style,label=label,lw=1.25)
            ax.set(title=rf'$\lambda={lam},\ k/\pi={ik/40:g}$',xlabel=r'$\omega/t$',ylabel=r'$tA_{0.25t}$',xlim=(-3.25,5.5),ylim=(0,None))
            ax.grid(alpha=.15)
    axes[0,0].legend(fontsize=7);fig.suptitle('Newly refined momenta: absolute spectral weight and statistical uncertainty')
    save(fig,out,'full-grid-problem-spectra')
    fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    for ax,lam in zip(axes,plan['couplings']):
        for ik in [8,20,32]:
            r=next(p for p in diag['comparisons'] if p['coupling']==lam and p['index']==ik)
            tau=np.array(next(p for p in diag['points'] if p['dataset']=='pilot' and p['coupling']==lam and p['index']==ik)['tau'])
            ax.plot(tau,r['sigma_ratio_pilot_over_baseline'],label=rf'$k/\pi={ik/40:g}$')
        ax.axhline(.5,color='k',ls='--',lw=1,label='Expected for 4x effective samples')
        ax.set(xlabel=r'$\tau t$',ylabel='New / original standard error',title=rf'$\lambda={lam}$',ylim=(0,None));ax.grid(alpha=.2)
    axes[0].legend(fontsize=8);fig.suptitle('Measured imaginary-time uncertainty')
    save(fig,out,'full-grid-green-errors')
    fig,axes=plt.subplots(1,2,figsize=(10,5),constrained_layout=True)
    rowkeys=[(lam,ik) for lam in plan['couplings'] for ik in [25,36,39]]
    for ax,eta in zip(axes,[.1,.25]):
        data=np.zeros((6,4))
        for j,(lam,ik) in enumerate(rowkeys):
            for i,name in enumerate(['double_015','double_030','double_060','double_090']):
                rec=next(p for p in resolution_cases if (p['coupling'],p['index'],p['case'])==(lam,ik,name))
                v=next(p for p in rec['results'] if p['eta']==eta)
                data[j,i]=v['doublet_recovered']/v['denominator'] if v['truth_doublet_visible'] else np.nan
                label=f"{v['doublet_recovered']}/24" if v['truth_doublet_visible'] else 'n/a'
                ax.text(i,j,label,ha='center',va='center',color='white' if data[j,i]<.5 else 'black',fontsize=9)
        cmap=plt.get_cmap('viridis').copy();cmap.set_bad('#eeeeee')
        image=ax.pcolormesh(np.arange(5)-.5,np.arange(7)-.5,data,vmin=0,vmax=1,
                            cmap=cmap,shading='flat',rasterized=False)
        ax.set_xlim(-.5,3.5);ax.set_ylim(5.5,-.5)
        ax.set(xticks=range(4),xticklabels=['0.15','0.30','0.60','0.90'],yticks=range(6),
            yticklabels=[f'lambda={lam}, k/pi={ik/40:g}' for lam,ik in rowkeys],xlabel='True doublet separation / t',title=rf'Evaluated at $\eta={eta}t$')
    bar=fig.colorbar(image,ax=axes,label='Empirical doublet recovery fraction',shrink=.8)
    bar.solids.set_rasterized(False)
    fig.suptitle('Reused pilot synthetic test; failures count as non-recovery\nn/a: truth itself does not resolve the doublet at this broadening')
    save(fig,out,'full-grid-resolution')
    fig,axes=plt.subplots(2,2,figsize=(10,6.5),constrained_layout=True)
    titles={'single':'Single narrow peak','double_030':r'Doublet, $\Delta\omega=0.30t$',
            'double_060':r'Doublet, $\Delta\omega=0.60t$','continuum':'Broad continuum'}
    for ax,name in zip(axes.flat,['single','double_030','double_060','continuum']):
        d=np.load(ROOT/f'results/resolution/lambda0.50_ik036_{name}.npz')
        if len(d['reconstructed']):
            lo,med,hi=np.quantile(d['reconstructed'][:,0],[.16,.5,.84],axis=0)
            ax.fill_between(d['energy'],lo,hi,color='#cb6238',alpha=.25,label='16-84% of recovered spectra')
            ax.plot(d['energy'],med,color='#cb6238',label='Reconstruction median')
        ax.plot(d['energy'],d['truth'][0],'--',color='#171717',label='Known synthetic truth')
        ax.set(title=titles[name],xlabel=r'$\omega/t$',ylabel=r'$tA_{0.1t}$',xlim=(-3.25,5.5),ylim=(0,None));ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8);fig.suptitle(r'Reused pilot synthetic examples: noise from $\lambda=0.5,\ k/\pi=0.9$')
    save(fig,out,'full-grid-synthetic-examples')
    print(json.dumps(public({k:report[k] for k in ['new_chains','total_chains','new_production_steps','total_production_steps','bootstrap_failures','synthetic_failures','selected_configurations','summary']}),indent=2))

if __name__=='__main__':main()
