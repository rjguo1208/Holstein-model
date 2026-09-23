"""Run analytic checks and independent real-space VED/DiagMC on allocated CPUs."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import numpy as np
from scipy.sparse.linalg import eigsh

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from defect_analysis import public, summarize, validate
from defect_reference import Cloud, DefectBasis, Parameters, lanczos


def reference(coupling,U,out):
    out.mkdir(parents=True,exist_ok=False)
    p=Parameters(g=np.sqrt(2*coupling),U=U)
    radius=24 if U<1 else 16
    runs=[]
    for nh,r in [(8,radius),(10,radius),(12,radius),(12,radius+8)]:
        start=time.monotonic();basis=DefectBasis(nh,r);result=basis.ground(p)
        vector=result.pop('vector')
        record=dict(generations=nh,radius=r,seconds=time.monotonic()-start,**result)
        runs.append(public(record))
        (out/'convergence.json').write_text(json.dumps(runs,indent=2)+'\n')
        print(f'REFERENCE lambda={coupling} U={U} Nh={nh} R={r} dim={basis.dim} E={result["E"]:.11f}',flush=True)
    bulk=float(eigsh(basis.cloud.bulk(p),k=1,which='SA',tol=1e-11,return_eigenvectors=False)[0])
    result.update(E_clean=bulk,binding=bulk-result['E'])
    energy,weight,alpha,beta=lanczos(basis.hamiltonian(p),basis.source(),240)
    moments=np.array([weight.sum(),energy@weight,energy**2@weight])
    expected=np.array([1,-U,U*U+2*p.t**2+p.g**2])
    np.testing.assert_allclose(moments,expected,atol=2e-10,rtol=0)
    np.savez_compressed(out/'reference.npz',**result,energies=energy,weights=weight,alpha=alpha,beta=beta)
    checks=dict(energy_cloud_change=abs(runs[2]['E']-runs[1]['E']),
        energy_radius_change=abs(runs[3]['E']-runs[2]['E']),
        Zb_cloud_change=abs(runs[2]['Zb']-runs[1]['Zb']),
        Zb_radius_change=abs(runs[3]['Zb']-runs[2]['Zb']),
        maximum_moment_residual=float(max(abs(moments-expected))))
    result.update(parameters=asdict(p),convergence=checks,generations=12,radius=radius+8)
    result['converged']=bool(checks['energy_cloud_change']<1e-5 and checks['energy_radius_change']<1e-5 and checks['Zb_cloud_change']<2e-4 and checks['Zb_radius_change']<2e-4 and result['residual']<1e-7)
    (out/'summary.json').write_text(json.dumps(public(result),indent=2)+'\n')
    return public(result)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--stage',choices=['pilot','scan'],default='pilot')
    parser.add_argument('--jobs',type=int,default=4)
    parser.add_argument('--chains',type=int,default=4)
    parser.add_argument('--blocks',type=int,default=80)
    parser.add_argument('--steps-per-block',type=int,default=100000)
    parser.add_argument('--seed-offset',type=int,default=0)
    parser.add_argument('--pilot',type=Path)
    args=parser.parse_args()
    if args.stage=='scan':
        if args.pilot is None: parser.error('Scan requires --pilot with passed checks')
        pilot=json.loads((args.pilot/'summary.json').read_text())
        if not pilot['validation_passed'] or not pilot['first_point_passed']: raise ValueError('Pilot validation has not passed')
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    sources=[p for name in ['src','scripts','tests'] for p in (ROOT/name).glob('*') if p.is_file()]
    sources+=list(ROOT.glob('*.py'))+[ROOT/'Makefile']
    manifest=dict(started_utc=datetime.now(timezone.utc).isoformat(),job_id=os.environ.get('SLURM_JOB_ID'),
        arguments={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},complete=False,
        executable_sha256=sha256((ROOT/'build/defect_diagmc').read_bytes()).hexdigest(),sources={},tasks=[])
    for p in sources:
        dest=out/'source'/p.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
        manifest['sources'][str(p.relative_to(ROOT))]=sha256(p.read_bytes()).hexdigest()
    # Execute the recorded binary even if development continues in the workspace.
    executable=out/'source/defect_diagmc';shutil.copyfile(ROOT/'build/defect_diagmc',executable);executable.chmod(0o755)
    def save(): (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    save()
    def run(task):
        name,parameters=task
        command=[str(executable)]
        for k,v in parameters.items(): command.extend(['--'+k,str(v)])
        dest=out/'chains'/name;dest.parent.mkdir(exist_ok=True)
        command.extend(['--out',str(dest)])
        with (dest.parent/(dest.name+'.log')).open('w') as stream:
            subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,check=True)
        print('Finished',name,flush=True)
        return dest
    def acquire(label,parameters,blocks,offset):
        tasks=[]
        for i in range(args.chains):
            p=dict(parameters,blocks=blocks,**{'steps-per-block':args.steps_per_block,'thin':100,'warmup':300000,'seed':730_000_001+offset+1009*i+args.seed_offset})
            tasks.append((f'{label}_s{i:02d}',p))
        manifest['tasks'] += [dict(name=n,parameters=p) for n,p in tasks];save()
        with ThreadPoolExecutor(max_workers=args.jobs) as pool: return list(pool.map(run,tasks))
    validations=[]
    if args.stage=='pilot':
        cases=[('free',dict(t=1,omega=1,g=0,U=0,mu=-2.3)),
            ('impurity',dict(t=1,omega=1,g=0,U=1,mu=-2.45)),
            ('atomic',dict(t=0,omega=1,g=.7,U=1,mu=-1.8)),
            ('one_phonon',dict(t=0,omega=1,g=.7,U=1,mu=-1.8,**{'max-arcs':1})),
            ('one_hop_pair',dict(t=1,omega=1,g=0,U=1,mu=-1.8,**{'max-hop-pairs':1}))]
        for i,(kind,p) in enumerate(cases):
            paths=acquire('check_'+kind,p|{'tau-max':12,'bins':48},32,100_000*i)
            check=validate(paths,out/'validation'/kind,kind);validations.append(check)
            print('VALIDATION',json.dumps(check),flush=True)
            (out/'validation.json').write_text(json.dumps(validations,indent=2)+'\n')
            if not check['passed']: raise RuntimeError('Analytic validation failed: '+kind)
    else:
        validations=pilot['validation']
    cases=[(.25,1.),(.5,1.)] if args.stage=='pilot' else [(lam,U) for lam in [.25,.5] for U in [.5,1.,2.]]
    report=dict(stage=args.stage,validation=validations,validation_passed=all(r['passed'] for r in validations),cases={},first_point_passed=False)
    for i,(lam,U) in enumerate(cases):
        label=f'lambda{lam:.2f}_U{U:g}'
        ref=reference(lam,U,out/'reference'/label)
        if not ref['converged']: raise RuntimeError('VED cutoffs not converged: '+label)
        T=(128 if U<1 else 64 if U<2 else 32) if args.stage=='scan' else 48
        # The reference energy sets an importance-sampling shift only; no
        # reference energy/residue enters the MC fit or projector estimator.
        gap=2/T if args.stage=='scan' else .10
        p=dict(t=1,omega=1,g=np.sqrt(2*lam),U=U,mu=ref['E']-gap,
               **{'tau-max':T,'bins':48,'max-arcs':256,'max-hop-pairs':512})
        paths=acquire(label,p,args.blocks,2_000_000+100_000*i)
        result=summarize(paths,out/'analysis'/label)
        point=result['primary'];prev=result['fit_windows'][-2]
        comparisons={}
        for key in ['E','Z0','Zb','p0','r2','nph']:
            difference=point[key]-ref[key]
            # A conservative acceptance screen; precision is reported separately.
            tolerance=max(5*point[key+'_error'],.003 if key=='E' else .02*max(abs(ref[key]),.1))
            comparisons[key]=dict(mc=point[key],error=point[key+'_error'],reference=ref[key],difference=difference,tolerance=tolerance,passed=abs(difference)<tolerance)
        passed=all(c['passed'] for c in comparisons.values()) and result['arc_cap_attempts']==result['hop_cap_attempts']==0
        report['cases'][label]=dict(reference=ref,diagmc=result,comparison=comparisons,passed=passed,
            projection_window_changes={key:point[key]-prev[key] for key in ['E','Zb','p0','r2']})
        report['first_point_passed']=all(c['passed'] for c in report['cases'].values())
        (out/'summary.json').write_text(json.dumps(public(report),indent=2)+'\n')
        print('COMPARISON',label,json.dumps(public(comparisons)),flush=True)
    manifest['complete']=True;manifest['finished_utc']=datetime.now(timezone.utc).isoformat();save()
    if not report['first_point_passed']: raise RuntimeError('First-point comparison needs review; scan is gated')


if __name__=='__main__': main()
