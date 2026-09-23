"""Parallel, immutable acquisition of denser local imaginary-time spectra."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--jobs',type=int,default=16)
    args=parser.parse_args()
    optimization=json.loads((ROOT/'results/defect_optimization_01/summary.json').read_text())
    if not optimization['complete'] or not all(r['trajectory_identical'] and r['green_blocks_bitwise_identical'] for r in optimization['cases']):raise ValueError('Optimization verification failed')
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    groundfile=ROOT/'results/defect_scan_01/summary.json';ground=json.loads(groundfile.read_text())
    if not ground['validation_passed'] or not ground['first_point_passed']:raise ValueError('Ground checks failed')
    shutil.copyfile(groundfile,out/'ground_input.json')
    source=out/'source';source.mkdir()
    paths=list((ROOT/'src').glob('*'))+[Path(__file__),ROOT/'scripts/python.sh',ROOT/'scripts/defect_dense_sample.slurm',ROOT/'Makefile']
    hashes={}
    for p in paths:
        if not p.is_file():continue
        name=p.relative_to(ROOT);dest=source/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
        hashes[str(name)]=sha256(dest.read_bytes()).hexdigest()
    executable=source/'defect_diagmc';shutil.copy2(ROOT/'build/defect_diagmc',executable)
    assert sha256(executable.read_bytes()).hexdigest()==sha256((ROOT/'results/defect_optimization_01/optimized/defect_diagmc').read_bytes()).hexdigest()
    tasks=[]
    site_order=[0,4,8]+[i for i in range(17) if i not in [0,4,8]]
    for site in site_order:
        for il,lam in enumerate([.25,.5]):
            for chain in range(8):
                p=dict(t=1,omega=1,g=float(np.sqrt(2*lam)),U=1,
                    mu=ground['cases'][f'lambda{lam:.2f}_U1']['reference']['E']-.10,origin=site,
                    **{'tau-max':12,'bins':192,'blocks':320,'steps-per-block':100000,'warmup':500000,
                       'thin':100,'radius':24,'max-arcs':256,'max-hop-pairs':512,'observables':'green',
                       'seed':960000001+10000000*il+100000*site+1009*chain})
                tasks.append(dict(name=f'lambda{lam:.2f}_i{site:02d}_s{chain:02d}',parameters=p))
    assert len({t['parameters']['seed'] for t in tasks})==272
    manifest=dict(complete=False,started_utc=datetime.now(timezone.utc).isoformat(),job_id=os.environ.get('SLURM_JOB_ID'),
        threads=args.jobs,measured_sites=list(range(17)),reflection_used=True,sources=hashes,tasks=tasks,
        ground_sha256=sha256(groundfile.read_bytes()).hexdigest(),executable_sha256=sha256(executable.read_bytes()).hexdigest(),
        production_steps=272*320*100000,completed_chains=0,
        note='New independent seeds; 4x production steps per position, 2x measured time bins. Extra sites extend spatial range, not lattice resolution.')
    def save():
        p=out/'manifest.tmp';p.write_text(json.dumps(manifest,indent=2)+'\n');p.replace(out/'manifest.json')
    save();(out/'chains').mkdir()
    def run(task):
        path=out/'chains'/task['name'];command=[str(executable)]
        for k,v in task['parameters'].items():command+=['--'+k,str(v)]
        command+=['--out',str(path)]
        with (path.parent/(path.name+'.log')).open('w') as stream:subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,check=True)
        m=json.loads((path/'run.json').read_text());assert m['complete'] and m['arc_cap_attempts']==m['hop_cap_attempts']==0
        for k,v in task['parameters'].items():assert m[k.replace('-','_')]==v
        assert (path/'blocks.csv').read_bytes().count(b'\n')==m['blocks']+1
        return task['name']
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for future in as_completed([pool.submit(run,task) for task in tasks]):
            name=future.result();manifest['completed_chains']+=1;save()
            print('COMPLETE',manifest['completed_chains'],'/272',name,flush=True)
    manifest.update(complete=True,finished_utc=datetime.now(timezone.utc).isoformat());save()


if __name__=='__main__':main()
