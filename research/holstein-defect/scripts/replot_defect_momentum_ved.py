"""Refresh plot labels from audited VED arrays without recalculating spectra."""
from hashlib import sha256
import json
from pathlib import Path
import shutil
import sys
import numpy as np
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from report_defect_momentum import draw
from defect_analysis import public

old = ROOT/'results/defect_momentum_ved_report_01'
out = ROOT/'results/defect_momentum_ved_report_02'
s = json.loads((old/'summary.json').read_text())
assert s['complete']
out.mkdir(exist_ok=False)
for f in old.glob('*.npz'):
    shutil.copyfile(f, out/f.name)
bundles = []
for lam in [.25, .5]:
    b = {}
    for U in [0, 1]:
        with np.load(out/f'lambda{lam:.2f}_U{U}_Nh12_R96_W32.npz') as f:
            b[U] = dict(f)
    bundles.append(b)
plt.rcParams.update({'font.size':10, 'mathtext.fontset':'stix', 'svg.fonttype':'path',
                     'svg.hashsalt':'holstein-defect-momentum-v1', 'savefig.facecolor':'white'})
s['plots'] = [draw(out, bundles, eta) for eta in [.1, .15, .25]]
s['plots'].append(draw(out, bundles, .15, difference=True))
s['plot_source_sha256'] = sha256((ROOT/'scripts/report_defect_momentum.py').read_bytes()).hexdigest()
s['original_report_sha256'] = sha256((old/'summary.json').read_bytes()).hexdigest()
s['plot_update'] = 'Close the math delimiter in clean-system panel titles; all numerical arrays unchanged.'
(out/'summary.json').write_text(json.dumps(public(s), indent=2)+'\n')
print('Replotted VED titles; numerical arrays unchanged.', flush=True)
