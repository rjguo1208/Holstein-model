"""Same-seed correctness and speed checks on an allocated compute node."""
from pathlib import Path
import json
import os
import subprocess
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/defect_optimization_01'
records=[]
for case,(lam,site) in enumerate([(.25,0),(.25,8),(.5,0),(.5,8)]):
    runs=[]
    for mode,binary in [('baseline',OUT/'baseline/defect_diagmc'),
                        ('all',OUT/'optimized/defect_diagmc'),('green',OUT/'optimized/defect_diagmc')]:
        path=OUT/f'case{case}_{mode}'
        command=[str(binary),'--out',str(path),'--g',str(np.sqrt(2*lam)),
            '--U','1','--mu',str(-2.606 if lam==.25 else -2.903),'--origin',str(site),
            '--tau-max','12','--bins','192','--blocks','30','--steps-per-block','100000',
            '--warmup','300000','--seed',str(930000001+case*1009)]
        if mode!='baseline':command+=['--observables',mode]
        subprocess.run(command,check=True)
        runs.append((json.loads((path/'run.json').read_text()),np.loadtxt(path/'blocks.csv',delimiter=',',skiprows=1)))
    old,full,green=runs
    for key in ['attempted','accepted','invalid_proposals','maximum_arcs','maximum_hops','maximum_extent']:
        assert old[0][key]==full[0][key]==green[0][key],key
    np.testing.assert_array_equal(old[1][:,:195],full[1][:,:195])
    np.testing.assert_array_equal(old[1][:,:195],green[1])
    np.testing.assert_allclose(old[1],full[1],rtol=2e-12,atol=2e-10)
    np.testing.assert_array_equal(np.loadtxt(OUT/f'case{case}_baseline/profiles.csv',delimiter=',',skiprows=1),
                                  np.loadtxt(OUT/f'case{case}_all/profiles.csv',delimiter=',',skiprows=1))
    r=dict(coupling=lam,site=site,baseline_seconds=old[0]['seconds'],all_seconds=full[0]['seconds'],
           green_seconds=green[0]['seconds'],trajectory_identical=True,green_blocks_bitwise_identical=True,
           maximum_absolute_difference_all_fields=float(np.max(abs(old[1]-full[1]))),
           green_speedup=old[0]['seconds']/green[0]['seconds'])
    records.append(r);print(json.dumps(r),flush=True)
summary=dict(complete=True,job_id=os.environ.get('SLURM_JOB_ID'),cases=records,
             combined_green_speedup=sum(r['baseline_seconds'] for r in records)/sum(r['green_seconds'] for r in records))
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
