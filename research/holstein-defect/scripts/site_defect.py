"""Publish the audited single-defect data, simple HTML and reproducible archives."""
from __future__ import annotations
import argparse
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
import shutil
from zipfile import ZipFile, ZIP_DEFLATED
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def files(directory):
    return [p for p in directory.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc']


def archive(destination, paths):
    # Published archives are immutable. To release changed content choose a new name.
    with ZipFile(destination,'x',compression=ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(set(paths)):
            z.write(p,Path('holstein-diagmc')/p.relative_to(ROOT))
    with ZipFile(destination) as z:
        assert z.testzip() is None
        count=len(z.infolist())
    assert destination.stat().st_size < 95_000_000
    return dict(file=destination.name,bytes=destination.stat().st_size,
                sha256=sha256(destination.read_bytes()).hexdigest(),members=count)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--site',type=Path,default=ROOT.parent/'Holstein-model')
    parser.add_argument('--report',type=Path,default=ROOT/'results/defect_report_02')
    parser.add_argument('--package',action='store_true')
    parser.add_argument('--archives-dir',type=Path,default=ROOT/'results/defect_release_01')
    args=parser.parse_args()
    s=json.loads((args.report/'summary.json').read_text())
    scan=json.loads((ROOT/'results/defect_scan_01/summary.json').read_text())
    assert s['total_chains']==148 and s['production_steps']==2816000000
    assert all(c['passed'] for c in scan['cases'].values())
    plots=args.site/'site/results';data=args.site/'site/data'
    for p in files(args.report):
        if p.suffix in ['.svg','.pdf','.png']:shutil.copyfile(p,plots/p.name)
    shutil.copyfile(args.report/'summary.json',data/'defect-summary.json')
    shutil.copyfile(ROOT/'results/defect_jobs_01.json',data/'defect-jobs.json')
    shutil.copyfile(ROOT/'DEFECT.md',data/'defect-method.md')
    downloads=[]
    for lam in [.25,.5]:
        for prefix,label in [('', 'DiagMC'),('reference_', 'VED')]:
            source=ROOT/f'results/defect_spectra_01/{prefix}lambda{lam:.2f}_map.npz'
            name=f'defect-{label.lower()}-lambda{lam:.2f}_map.npz'
            shutil.copyfile(source,data/name)
            with np.load(source) as d:
                assert d['A'].shape==(3,17,801)
                xx,ww=np.meshgrid(d['x'],d['energy'],indexing='ij')
                csvname=name.replace('.npz','_eta0.25.csv.gz')
                with gzip.open(data/csvname,'wt') as f:
                    np.savetxt(f,np.column_stack([xx.ravel(),ww.ravel(),d['A'][0].ravel()]),
                        delimiter=',',header='site,energy_over_t,t_A_eta0.25',comments='',fmt='%.12g')
            downloads.append(f'<li>{label}，λ={lam:g}：<a href="data/{name}">NPZ</a> · <a href="data/{csvname}">η=0.25t 的 CSV（gzip）</a></li>')
    with (data/'defect-ground.csv').open('w') as f:
        keys=['coupling','U','E','E_error','binding','binding_error','p0','p0_error','Zb','Zb_error','Z0_projected','Z0_projected_error','rms_radius','rms_radius_error','nph','nph_error']
        writer=csv.DictWriter(f,fieldnames=keys,extrasaction='ignore');writer.writeheader();writer.writerows(s['cases'].values())
    def pm(c,key):return f'{c[key]:.4f} ± {c[key+"_error"]:.4f}'
    ground=[];comparison=[];local=[]
    for c in s['cases'].values():
        ground.append(f'<tr><td>{c["coupling"]:g}</td><td>{c["U"]:g}</td>'+''.join(f'<td>{pm(c,k)}</td>' for k in ['E','binding','p0','Zb'])+'</tr>')
        ref=c['reference'];diff=c['E']-ref['E']
        comparison.append(f'<tr><td>{c["coupling"]:g}</td><td>{c["U"]:g}</td><td>{ref["E"]:.7f}</td><td>{diff:+.5f}</td><td>{abs(diff)/c["E_error"]:.2f}</td></tr>')
        if c['U']==1:
            x=c['localization']['fits'][0]
            local.append(f'<tr><td>{c["coupling"]:g}</td><td>{pm(c,"Z0_projected")}</td><td>{pm(c,"rms_radius")}</td><td>{x["xi"]:.3f} ± {x["error"]:.3f}</td><td>{x["reference"]:.3f}</td></tr>')
    spectral=[]
    for lam,c in s['spectra'].items():
        spectral.append(f'<tr><td>{lam}</td>'+''.join(f'<td>{min(e["relative_l1"]):.1%}–{max(e["relative_l1"]):.1%}</td>' for e in c['comparison'])+f'<td>{c["maximum_cv_score"]:.2f}</td></tr>')
    def figure(stem,alt,caption,height=450):
        return f'<figure><div class="figure-scroll data"><img src="results/{stem}.svg" width="1020" height="{height}" alt="{alt}"></div><figcaption>{caption}</figcaption><p class="figure-links"><a href="results/{stem}.svg">SVG</a><a href="results/{stem}.pdf">PDF</a><a href="results/{stem}.png">PNG</a></p></figure>'
    html=(ROOT/'scripts/defect_page.html').read_text()
    values={
        'GROUND':''.join(ground),'COMPARISON':''.join(comparison),'LOCAL':''.join(local),'SPECTRAL':''.join(spectral),'DOWNLOADS':''.join(downloads),
        'BINDING':figure('defect-binding','两种耦合下束缚能与总零声子权重随吸引势变化，DiagMC误差条与VED参考比较。','点及误差条为 DiagMC；线连接相同参数处的独立 VED 结果。'),
        'DENSITY':figure('defect-density','缺陷位于零点；电子密度与局域谱极点权重是不同的曲线，分别比较DiagMC与VED。',r'\(U=t\)：蓝色为电子密度 \(p_i\)，棕色为局域束缚态极点权重 \(Z_i\)。两者都直接测量；本图保留正、负位置独立统计的差异。'),
        'FINE':figure('defect-ldos-eta0.25','单缺陷局域谱图：横轴位置，纵轴能量，颜色为绝对谱强度。上排DiagMC初步重建，下排VED参考。',r'\(U=t,\ \eta=0.25t\)。四图使用同一色标和真空能量零点。上排出现的部分细峰没有得到 VED 支持，不能据此确认新的激发。',840),
        'COARSE':figure('defect-ldos-eta1','同一局域谱的较粗展宽图，上排DiagMC与下排VED使用共同色标。',r'\(\eta=t\)：同一组重建谱增加显示展宽，宽谱包络与 VED 更接近。没有重新拟合或逐位置归一化。',840),
    }
    for key,value in values.items():html=html.replace('@'+key+'@',value)
    release='https://github.com/rjguo1208/Holstein-model/releases/download/single-defect-20260923/'
    for name in ['completed','ground-raw','spectra-raw','ground-green','spectra-green']:
        filename=f'defect-{name}-01.zip'
        html=html.replace('data/'+filename,release+filename)
    assert not any('@'+key+'@' in html for key in values)
    (args.site/'src/defect.html').write_text(html)
    if args.package:
        args.archives_dir.mkdir(parents=True,exist_ok=True)
        source=list(ROOT.glob('*.py'))+list(ROOT.glob('*.md'))+[ROOT/'Makefile',ROOT/'pytest.ini']
        for name in ['src','scripts','tests']:source+=files(ROOT/name)
        mainpaths=source+[ROOT/'results/defect_jobs_01.json']
        raw_ground=[];raw_spectra=[];green_ground=[];green_spectra=[]
        for name in ['defect_pilot_01','defect_scan_01','defect_spectra_01','defect_projector_check_01',args.report.name]:
            for p in files(ROOT/'results'/name):
                if p.name=='green.npz':
                    (green_spectra if name=='defect_spectra_01' else green_ground).append(p)
                elif 'chains' in p.relative_to(ROOT/'results'/name).parts:
                    (raw_spectra if name=='defect_spectra_01' else raw_ground).append(p)
                else:mainpaths.append(p)
        for c in s['clean_diagmc'].values():mainpaths.append(Path(c['source']))
        for jid in ['20881037','20881090','20881164']:mainpaths+=list((ROOT/'logs').glob('*'+jid+'.out'))
        records=[]
        for name,paths in [('defect-completed-01.zip',mainpaths),('defect-ground-raw-01.zip',raw_ground),
                           ('defect-spectra-raw-01.zip',raw_spectra),('defect-ground-green-01.zip',green_ground),
                           ('defect-spectra-green-01.zip',green_spectra)]:
            record=archive(args.archives_dir/name,paths)
            record.update(url=release+name,zip_integrity_verified=True)
            records.append(record);print(record,flush=True)
        (data/'defect-archives.json').write_text(json.dumps(dict(root='holstein-diagmc',archives=records),indent=2)+'\n')
    print('Generated defect.html, four figures, ground CSV and local spectra.')


if __name__=='__main__':main()
