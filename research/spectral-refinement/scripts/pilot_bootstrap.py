"""128 block bootstrap draws, reselecting alpha within the frozen representation."""
import concurrent.futures as cf
from hashlib import sha256
import json
import os
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from analysis import read_runs
from refinement import cross_validate
from pilot_compare import curves,AXIS,ETAS
from pilot_diagnostics import public

def run(task):
    lam,ik=task
    path=ROOT/f'results/comparison/lambda{lam:.2f}_ik{ik:03d}.json'
    record=json.loads(path.read_text())
    config=record['primary_cv']['config']
    paths=sorted((ROOT/'results/pilot_data').glob(f'ik{ik:03d}/lambda{lam:.2f}_*/run.json'))
    rows,owners,metas=read_runs([p.parent for p in paths],'rb_blocks.csv')
    # Moving whole groups of 16 raw blocks protects within-group correlations.
    # This changes the resampling unit, not the error bars used by model selection.
    # The final report includes the independent reblocking/chain diagnostics.
    chunk=16
    rng=np.random.default_rng(912731+int(lam*100)*1000+ik)
    samples=[];alphas=[];failures=[]
    for rep in range(128):
        blocks=[]
        for owner in np.unique(owners):
            group=rows[owners==owner]
            assert len(group)%chunk==0
            chunks=group.reshape(-1,chunk,group.shape[1])
            blocks.append(chunks[rng.integers(len(chunks),size=len(chunks))].reshape(group.shape))
        sample=np.concatenate(blocks)
        try:
            fit,cv=cross_validate(sample,owners,metas[0],config)
            samples.append(curves(fit));alphas.append(fit['alpha'])
        except (ValueError,RuntimeError,np.linalg.LinAlgError) as error:
            failures.append(dict(replicate=rep,error=str(error)))
        if (rep+1)%32==0:print(f'bootstrap lambda={lam} ik={ik:03d}: {rep+1}/128',flush=True)
    out=ROOT/'results/bootstrap'
    if not samples:raise ValueError('All bootstrap replicates failed')
    samples=np.asarray(samples)
    quantiles=np.quantile(samples,[.025,.16,.84,.975],axis=0)
    np.savez_compressed(out/f'lambda{lam:.2f}_ik{ik:03d}.npz',energy=AXIS,eta=ETAS,
        samples=samples,quantiles=quantiles,probabilities=[.025,.16,.84,.975],alphas=alphas)
    result=dict(coupling=lam,index=ik,requested=128,successful=len(samples),failures=failures,
        resampling='within each independent chain, nonoverlapping groups of 16 raw blocks',
        resampled_raw_steps=1600000,configuration=config,alpha_reselected_each_draw=True,
        configuration_reselected=False,comparison_sha256=sha256(path.read_bytes()).hexdigest())
    (out/f'lambda{lam:.2f}_ik{ik:03d}.json').write_text(json.dumps(public(result),indent=2)+'\n')
    return result

if __name__=='__main__':
    plan=json.loads((ROOT/'plans/pilot.json').read_text())
    manifest=json.loads((ROOT/'results/comparison/provenance.json').read_text())
    assert manifest['complete']
    out=ROOT/'results/bootstrap';out.mkdir(exist_ok=False)
    tasks=[(lam,ik) for lam in plan['couplings'] for ik in plan['indices']]
    with cf.ProcessPoolExecutor(max_workers=min(20,int(os.environ.get('SLURM_CPUS_PER_TASK','1')))) as pool:
        result=list(pool.map(run,tasks))
    (out/'summary.json').write_text(json.dumps(public(dict(complete=True,points=result,
        job_id=os.environ.get('SLURM_JOB_ID'))),indent=2)+'\n')
