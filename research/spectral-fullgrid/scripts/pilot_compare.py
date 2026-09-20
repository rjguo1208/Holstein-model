"""Compare continuation representations on the new chains; never load VED."""
import concurrent.futures as cf
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import sys
import traceback
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from analysis import read_runs
from map_continuation import lorentz_map
from refinement import candidates,cross_validate
from pilot_diagnostics import public

AXIS=np.linspace(-3.25,5.5,701)
ETAS=np.array([.25,.5,1.])

def curves(fit):
    return np.array([lorentz_map(AXIS,fit['energies'],fit['weights'],eta) for eta in ETAS])

def run(task):
    lam,ik=task
    paths=sorted((ROOT/'results/pilot_data').glob(f'ik{ik:03d}/lambda{lam:.2f}_*/run.json'))
    rows,owners,metas=read_runs([p.parent for p in paths],'rb_blocks.csv')
    p=metas[0];records=[];errors=[]
    for config in candidates():
        try:
            fit,cv=cross_validate(rows,owners,p,config)
            records.append(dict(name=config['name'],fit=fit,cv=cv))
            print(f'lambda={lam} ik={ik:03d} {config["name"]}: CV={cv["score"]:.5g}',flush=True)
        except (ValueError,RuntimeError,np.linalg.LinAlgError) as error:
            errors.append(dict(config=config,error=str(error),traceback=traceback.format_exc()))
            print(f'FAILED candidate lambda={lam} ik={ik:03d} {config["name"]}: {error}',flush=True)
    baseline=next(r for r in records if r['name']=='baseline')
    best=min(records,key=lambda r:r['cv']['score'])
    base_scores=baseline['cv']['folds'][:,baseline['cv']['selected']]
    best_scores=best['cv']['folds'][:,best['cv']['selected']]
    difference=base_scores-best_scores
    benefit=float(difference.mean());se=float(difference.std(ddof=1)/np.sqrt(len(difference)))
    primary=best if benefit>se else baseline
    # The spread of near-optimal representations is a sensitivity diagnostic,
    # not a posterior distribution or a calibrated confidence interval.
    threshold=best['cv']['score']+best['cv']['se']
    eligible=[r for r in records if r['cv']['score']<=threshold]
    if primary['name'] not in [r['name'] for r in eligible]:eligible.append(primary)
    result=dict(coupling=lam,index=ik,k=p['k'],parameters=p,primary=primary['name'],
        primary_fit=primary['fit'],primary_cv=primary['cv'],candidates=records,failures=errors,
        selection=dict(rule='one-SE alpha in each representation; switch from baseline only if paired CV gain exceeds one fold-SE',
                       gain=benefit,paired_se=se,best_candidate=best['name']),
        eligible_candidates=[r['name'] for r in eligible],
        input_hashes={str(path.parent/name):sha256((path.parent/name).read_bytes()).hexdigest()
                      for path in paths for name in ['run.json','rb_blocks.csv']},
        spectral_resolution_validated=False)
    out=ROOT/'results/comparison'
    (out/f'lambda{lam:.2f}_ik{ik:03d}.json').write_text(json.dumps(public(result),indent=2)+'\n')
    np.savez_compressed(out/f'lambda{lam:.2f}_ik{ik:03d}.npz',energy=AXIS,eta=ETAS,
        primary=curves(primary['fit']),sampling_only=curves(baseline['cv']['minimum_fit']),
        candidates=np.array([curves(r['fit']) for r in records]),names=[r['name'] for r in records],
        sensitivity_low=np.min([curves(r['fit']) for r in eligible],axis=0),
        sensitivity_high=np.max([curves(r['fit']) for r in eligible],axis=0))
    return dict(coupling=lam,index=ik,primary=primary['name'],score=primary['cv']['score'],
                baseline_min_score=float(baseline['cv']['mean'][baseline['cv']['best']]),failures=len(errors))

if __name__=='__main__':
    plan=json.loads((ROOT/'plans/pilot.json').read_text())
    out=ROOT/'results/comparison';out.mkdir(exist_ok=False)
    sources=[p for p in ROOT.glob('*.py')]+[p for p in (ROOT/'scripts').glob('*') if p.is_file()]
    hashes={}
    for path in sources:
        rel=path.relative_to(ROOT);dest=out/'source'/rel;dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(path,dest);hashes[str(rel)]=sha256(path.read_bytes()).hexdigest()
    manifest=dict(complete=False,job_id=os.environ.get('SLURM_JOB_ID'),source_hashes=hashes,
                  reference_used_in_inference=False,candidates=candidates())
    (out/'provenance.json').write_text(json.dumps(public(manifest),indent=2)+'\n')
    tasks=[(lam,ik) for lam in plan['couplings'] for ik in plan['indices']]
    with cf.ProcessPoolExecutor(max_workers=min(20,int(os.environ.get('SLURM_CPUS_PER_TASK','1')))) as pool:
        records=list(pool.map(run,tasks))
    manifest.update(complete=True,points=records)
    (out/'provenance.json').write_text(json.dumps(public(manifest),indent=2)+'\n')
