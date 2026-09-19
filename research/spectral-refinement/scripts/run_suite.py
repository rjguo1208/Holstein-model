"""Run bounded independent chains on the allocated compute node."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from analysis import summarize, summarize_conditional, validate_run

parser=argparse.ArgumentParser()
parser.add_argument('suite',choices=['validation','k0','dispersion','controls','spectral_data','spectral_map'])
parser.add_argument('--out',required=True)
parser.add_argument('--jobs',type=int,default=4)
parser.add_argument('--momenta',type=float,nargs='+',default=[.15,.30,.45])
parser.add_argument('--couplings',type=float,nargs='+',default=[.25,.5])
parser.add_argument('--omega',type=float,default=1.)
parser.add_argument('--chains',type=int)
parser.add_argument('--blocks',type=int,default=240)
parser.add_argument('--steps-per-block',type=int,default=200000)
parser.add_argument('--seed-offset',type=int,default=0)
parser.add_argument('--rb-thin',type=int,default=100)
parser.add_argument('--map-points',type=int,default=17)
parser.add_argument('--map-index',type=int)
args=parser.parse_args()
if args.suite=='spectral_map':
    if args.map_points < 2 or (args.map_index is not None and not 0 <= args.map_index < args.map_points):
        parser.error('Require map-points >= 2 and 0 <= map-index < map-points')
output=Path(args.out).resolve()
if output.exists(): raise SystemExit('Refusing to overwrite an existing suite')
output.mkdir(parents=True)
sources=[p for folder in ['src','scripts'] for p in (ROOT/folder).glob('*') if p.is_file()]+[ROOT/'analysis.py',ROOT/'Makefile']
manifest=dict(suite=args.suite,python=sys.version,job_id=os.environ.get('SLURM_JOB_ID'),
    source_hashes={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sources},
    executable_sha256=sha256((ROOT/'build/diagmc').read_bytes()).hexdigest())
for p in sources:
    dest=output/'source'/p.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)

tasks=[]
if args.suite=='validation':
    for kind,params in [
        ('free',dict(t=1,omega=1,g=0,k=.37,mu=-2.1,**{'tau-max':12,'max-order':32,'bins':48})),
        ('atomic',dict(t=0,omega=1,g=.7,k=0,mu=-.65,**{'tau-max':16,'max-order':64,'bins':64})),
        ('first_order',dict(t=1,omega=1,g=.7,k=.45,mu=-2.2,**{'tau-max':12,'max-order':1,'bins':48}))]:
        tasks.append((kind,params|{'blocks':160,'steps-per-block':50000,'warmup':100000,'thin':5,'rb-thin':100,'seed':11939}))
else:
    for coupling in args.couplings:
        mu={.25:-2.32,.5:-2.58}.get(coupling,-2-1.2*coupling-.1)
        momenta=[0.] if args.suite!='dispersion' else args.momenta
        if args.suite=='spectral_map':
            from math import pi
            indices=range(args.map_points) if args.map_index is None else [args.map_index]
            momenta=[pi*i/(args.map_points-1) for i in indices]
        for k in momenta:
            groups=[('base',{})] if args.suite!='controls' else [
                ('mu',{'mu':mu-.15}),('tau',{'tau-max':32,'bins':128}),('order',{'max-order':48})]
            for label,extra in groups:
                seeds=([101,211,307,401] if args.suite in ['k0','dispersion','spectral_data','spectral_map'] else [1009,2003])
                if args.chains is not None:seeds=[17011+1009*i for i in range(args.chains)]
                for seed in seeds:
                    name=f'lambda{coupling:.2f}_k{k:.2f}_{label}_s{seed}'
                    if args.suite=='spectral_map':
                        ik=round(k*(args.map_points-1)/pi)
                        name=f'lambda{coupling:.2f}_ik{ik:03d}_{label}_s{seed}'
                    params={'t':1,'omega':args.omega,'g':(2*coupling*args.omega)**.5,'k':k,'mu':mu,
                        'tau-max':24,'max-order':max(96,int(128*coupling)),'bins':96,'blocks':args.blocks,
                        'steps-per-block':args.steps_per_block,'warmup':500000,'thin':5,'seed':seed}|extra
                    if args.suite=='dispersion':
                        params['seed']=seed+round(k*1000)*10000+round(coupling*100)*1000000
                    if args.suite in ['spectral_data','spectral_map']:
                        params.update({'mu':min(-2.8,mu-.1),'tau-max':12,'bins':96,'rb-thin':args.rb_thin,
                                       'seed':seed+60000000+round(coupling*100)*100000})
                        if args.suite=='spectral_map':
                            params['seed']+=100000000+ik*1000000
                    params['seed']+=args.seed_offset
                    tasks.append((name,params))
manifest['tasks']=[dict(name=name,parameters=parameters) for name,parameters in tasks]
(output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

def run(task):
    name,params=task
    command=[str(ROOT/'build/diagmc')]
    for key,value in params.items(): command.extend(['--'+key,str(value)])
    command.extend(['--out',str(output/name)])
    with (output/(name+'.log')).open('w') as log:
        subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    print('Finished',name,flush=True)
    return output/name

with ThreadPoolExecutor(max_workers=args.jobs) as pool:
    paths=list(pool.map(run,tasks))
if args.suite=='validation':
    results=[validate_run(path,path.name) for path in paths]
    results += [dict(validate_run(path,path.name,'rb_blocks.csv'),measurement='conditional') for path in paths]
    (output/'validation.json').write_text(json.dumps(results,indent=2)+'\n')
    if not all(r['passed'] for r in results): raise SystemExit('Analytic validation failed')
    print(json.dumps(results,indent=2),flush=True)
else:
    groups={}
    for path in paths: groups.setdefault(path.name.rsplit('_s',1)[0],[]).append(path)
    # Large-k states need not have a well-resolved isolated ground pole.
    # A spectral-map acquisition therefore saves G and its covariance without
    # forcing an exponential E/Z tail fit that can fail or bias continuation.
    if args.suite=='spectral_map':
        results={name:summarize_conditional(group,output/'conditional'/name) for name,group in groups.items()}
    else:
        results={name:summarize(group,output/'analysis'/name) for name,group in groups.items()}
    if args.suite=='spectral_data':
        for name,group in groups.items(): summarize_conditional(group,output/'conditional'/name)
    (output/'summary.json').write_text(json.dumps(results,indent=2)+'\n')
    for name,result in results.items(): print(name,result.get('primary',result),flush=True)
