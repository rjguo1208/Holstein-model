"""Recover only failed bootstrap draws, preserving their RNG draws and CV rule."""
import concurrent.futures as cf
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from analysis import read_runs
from refinement import cross_validate
from pilot_compare import curves
from pilot_diagnostics import public

def repair(task):
    lam,ik,failed=task
    selected=json.loads((ROOT/f'results/comparison/lambda{lam:.2f}_ik{ik:03d}.json').read_text())
    config=selected['primary_cv']['config']
    paths=sorted((ROOT/'results/pilot_data').glob(f'ik{ik:03d}/lambda{lam:.2f}_*/run.json'))
    rows,owners,metas=read_runs([p.parent for p in paths],'rb_blocks.csv')
    rng=np.random.default_rng(912731+int(lam*100)*1000+ik);recovered={}
    for rep in range(max(failed)+1):
        blocks=[]
        for owner in np.unique(owners):
            group=rows[owners==owner];chunks=group.reshape(-1,16,group.shape[1])
            blocks.append(chunks[rng.integers(len(chunks),size=len(chunks))].reshape(group.shape))
        if rep not in failed:continue
        sample=np.concatenate(blocks)
        fit,cv=cross_validate(sample,owners,metas[0],config,recover=True)
        assert fit['converged'] and max(abs(fit['moment_residuals']))<1e-6
        recovered[rep]=dict(curves=curves(fit),alpha=fit['alpha'],fit=fit,
            resampled_rows_sha256=sha256(sample.tobytes()).hexdigest(),selected_cv=cv['selected'])
        print(f'Recovered lambda={lam} ik={ik:03d} replicate={rep}: '
              f'alpha={fit["alpha"]}, gradient={fit["dual_gradient"]:.4g}, '
              f'moment error={max(abs(fit["moment_residuals"])):.4g}',flush=True)
    return lam,ik,recovered

if __name__=='__main__':
    summary=ROOT/'results/bootstrap/summary.json';initial=json.loads(summary.read_text())
    tasks=[(p['coupling'],p['index'],[f['replicate'] for f in p['failures']]) for p in initial['points'] if p['failures']]
    if not tasks:
        print('No failed bootstrap draws require repair',flush=True);sys.exit(0)
    backup=ROOT/'results/bootstrap_initial';backup.mkdir(exist_ok=False)
    shutil.copyfile(summary,backup/'summary.json')
    for lam,ik,_ in tasks:
        for ext in ['json','npz']:
            name=f'lambda{lam:.2f}_ik{ik:03d}.{ext}';shutil.copyfile(ROOT/'results/bootstrap'/name,backup/name)
    with cf.ProcessPoolExecutor(max_workers=min(len(tasks),int(os.environ.get('SLURM_CPUS_PER_TASK','1')))) as pool:
        results=list(pool.map(repair,tasks))
    recoveries=[]
    for lam,ik,recovered in results:
        stem=f'lambda{lam:.2f}_ik{ik:03d}'
        with np.load(backup/f'{stem}.npz') as saved:data={k:saved[k] for k in saved.files}
        failed=set(recovered);successful=[r for r in range(128) if r not in failed]
        spectra=np.empty((128,)+data['samples'].shape[1:]);alphas=np.empty(128)
        spectra[successful]=data['samples'];alphas[successful]=data['alphas']
        for rep,record in recovered.items():spectra[rep]=record['curves'];alphas[rep]=record['alpha']
        np.testing.assert_array_equal(spectra[successful],data['samples'])
        data.update(samples=spectra,alphas=alphas,quantiles=np.quantile(spectra,[.025,.16,.84,.975],axis=0))
        np.savez_compressed(ROOT/'results/bootstrap'/f'{stem}.npz',**data)
        meta=json.loads((backup/f'{stem}.json').read_text())
        meta.update(initial_failures=meta['failures'],failures=[],successful=128,
            numerical_recoveries=public({str(rep):{k:v for k,v in r.items() if k!='curves'} for rep,r in recovered.items()}))
        (ROOT/'results/bootstrap'/f'{stem}.json').write_text(json.dumps(meta,indent=2)+'\n')
        entry=next(p for p in initial['points'] if p['coupling']==lam and p['index']==ik);entry.update(meta)
        recoveries.append(dict(coupling=lam,index=ik,replicates=sorted(recovered)))
    initial['recovery_job_id']=os.environ.get('SLURM_JOB_ID');initial['numerical_recoveries']=recoveries
    summary.write_text(json.dumps(initial,indent=2)+'\n')
    (ROOT/'results/bootstrap_recovery.json').write_text(json.dumps(dict(complete=True,
        job_id=os.environ.get('SLURM_JOB_ID'),recoveries=recoveries,
        source_hashes={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in [ROOT/'refinement.py',Path(__file__)]},
        unchanged_central_selections=True,unchanged_successful_bootstrap_draws=True,
        reference_used=False),indent=2)+'\n')
