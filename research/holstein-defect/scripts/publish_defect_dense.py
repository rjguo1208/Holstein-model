"""Publish completed local-spectrum refinement in both website languages."""
from __future__ import annotations
import argparse
from hashlib import sha256
import gzip
import json
from pathlib import Path
import shutil
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--website',type=Path,required=True)
    ap.add_argument('--report',type=Path,default=ROOT/'results/defect_dense_report_01')
    ap.add_argument('--resolution',type=Path,default=ROOT/'results/defect_dense_resolution_01');a=ap.parse_args()
    report=json.loads((a.report/'summary.json').read_text());assert report['complete']
    resolution=json.loads((a.resolution/'summary.json').read_text());assert resolution['complete']
    repo=a.website;data=repo/'site/data';plots=repo/'site/results'
    for p in a.report.glob('*'):
        if p.suffix in ['.svg','.pdf','.png']:shutil.copyfile(p,plots/p.name)
    shutil.copyfile(a.report/'summary.json',data/'defect-dense-summary.json')
    shutil.copyfile(a.resolution/'summary.json',data/'defect-dense-resolution.json')
    shutil.copyfile(ROOT/'DEFECT_DENSE.md',data/'defect-dense-method.md')
    shutil.copyfile(ROOT/'results/defect_dense_job_records_01.json',data/'defect-dense-jobs.json')
    for lam in [.25,.5]:
        for prefix,label in [('', 'diagmc'),('reference_','ved')]:
            source=a.report/f'{prefix}lambda{lam:.2f}_map.npz';name=f'defect-dense-{label}-lambda{lam:.2f}'
            shutil.copyfile(source,data/(name+'.npz'));d=np.load(source)
            for eta in [.15,.25]:
                ie=np.flatnonzero(d['eta']==eta).item();xx,ee=np.meshgrid(d['x'],d['energy'],indexing='ij')
                with gzip.open(data/f'{name}-eta{eta:g}.csv.gz','wt') as stream:
                    np.savetxt(stream,np.column_stack([xx.ravel(),ee.ravel(),d['A'][ie].ravel()]),delimiter=',',header='site,energy_over_t,tA',comments='')
        shutil.copyfile(a.report/f'lambda{lam:.2f}_diagnostics.npz',data/f'defect-dense-diagnostics-lambda{lam:.2f}.npz')
    catalog=json.loads((repo/'src/locales/en.json').read_text())
    def tr(zh,en):catalog[' '.join(zh.split())]=en;return zh
    def p(zh,en):return '<p>'+tr(zh,en)+'</p>\n'
    def figure(stem,zh,en,alt,alt_en,width=1040,height=890):
        return ('<figure><div class="figure-scroll data"><img src="results/'+stem+'.svg" width="'+str(width)+'" height="'+str(height)+'" alt="'+tr(alt,alt_en)+'"></div><figcaption>'+tr(zh,en)+'</figcaption><p class="figure-links">'+''.join('<a href="results/'+stem+'.'+ext+'">'+ext.upper()+'</a>' for ext in ['svg','pdf','png'])+'</p></figure>\n')
    s='<section id="spectra"><h2>'+tr('加密局域谱图：位置 × 能量 × 谱强度','Denser local spectra: position × energy × spectral intensity')+'</h2>\n'
    s+=r'''\[A_i(\omega)=\sum_\nu|\langle\nu|c_i^\dagger|\mathrm{vac}\rangle|^2\delta(\omega-E_\nu),\qquad
\mathcal G_{ii,\mu}(\tau)=\int d\omega\,A_i(\omega)e^{-(\omega-\mu)\tau}.\]
'''
    s+=p(r'本轮优先加密单缺陷局域谱：每种耦合独立测量 \(i=0,\ldots,16\) 共17个位置，反射后显示33列；位置间不插值。横轴为格点，纵轴为相对真空的能量，颜色为绝对强度 \(tA_{i,\eta}\)。增加位置扩展了空间范围，格距仍为 \(a=1\)。',r'This refinement focuses on the single-defect local spectrum: 17 sites \(i=0,\ldots,16\) are measured independently for each coupling and reflected to display 33 columns, without interpolation between sites. The horizontal axis is lattice position, the vertical axis is energy relative to vacuum, and color gives the absolute intensity \(tA_{i,\eta}\). More sites extend the spatial range; the spacing remains \(a=1\).')
    s+=p(r'每个位置由4条链增至8条链，每条生产步从1600万增至3200万；每个位置的生产统计增至四倍。虚时间箱从96增至192，\(\tau_{\max}=12/t\)；绘图能量点从801增至1601。能量网格加密和减小展宽本身不代表真实谱峰分辨率已经提高。',r'Each site increases from four to eight chains and from 16 million to 32 million production steps per chain, giving four times the production statistics per site. Imaginary-time bins increase from 96 to 192 with \(\tau_{\max}=12/t\); plotted energy points increase from 801 to 1601. A denser energy grid and smaller broadening do not by themselves establish improved peak resolution.')
    s+=r'''\[A_{i,\eta}(\omega)=\int d\omega'\,A_i(\omega')\frac{\eta/\pi}{(\omega-\omega')^2+\eta^2}.\]
'''
    s+=p(r'白色虚线为裸电子连续带边界 \(\omega=\pm2t\)，蓝色虚线为 \(g=0\) 时单缺陷的束缚能 \(E_b=-\sqrt{4t^2+U^2}\)。局域谱横轴是位置；余弦裸电子色散 \(\varepsilon_k=-2t\cos k\) 已加到干净体系的动量谱主图中。',r'White dashed lines mark the bare-electron continuum edges \(\omega=\pm2t\); the blue dashed line marks the single-defect bound energy \(E_b=-\sqrt{4t^2+U^2}\) at \(g=0\). The local spectrum has position on its horizontal axis. The cosine bare-electron dispersion \(\varepsilon_k=-2t\cos k\) is overlaid on the main clean-system momentum maps.')
    for eta in [.25,.15,1.]:
        zh=rf'\(U=t,\ \eta={eta:g}t\)。上排为新 DiagMC 重建，下排为独立 VED；四图共用色标、展宽和能量零点。'
        en=rf'\(U=t,\ \eta={eta:g}t\). The upper row shows the new DiagMC reconstruction and the lower row independent VED. All four panels share the color scale, broadening and energy origin.'
        if eta==.15:zh+='较小展宽用于探索，窄峰尚未通过分辨率验证。';en+=' This smaller broadening is exploratory; narrow peaks have not passed resolution validation.'
        s+=figure(f'defect-dense-ldos-eta{eta:g}',zh,en,
            f'加密单缺陷局域谱，展宽{eta:g}t；横轴为33个格点，虚线标出裸电子带边界与无声子缺陷束缚能。',
            f'Denser local spectrum at broadening {eta:g}t with 33 lattice sites; dashed lines mark bare continuum edges and the phonon-free defect bound energy.')
    s+='<div class="table-scroll"><table><caption>'+tr('全部17个独立测量位置：整谱相对 L1 差异范围，DiagMC 与 VED','All 17 independently measured sites: range of full-spectrum relative L1 differences between DiagMC and VED')+'</caption><thead><tr><th scope="col">λ</th>'+''.join(f'<th scope="col">η={e:g}t</th>' for e in [.15,.25,.5,1.])+'</tr></thead><tbody>'
    for lam,c in report['couplings'].items():
        v=np.array(c['relative_l1'])*100
        s+='<tr><td>'+lam+'</td>'+''.join(f'<td>{row.min():.1f}%–{row.max():.1f}%</td>' for row in v)+'</tr>'
    s+='</tbody></table></div>\n'
    s+=r'''\[\delta_i(\eta)=\frac{\int_{-4.5t}^{5.5t}|A^{\mathrm{MC}}_{i,\eta}-A^{\mathrm{VED}}_{i,\eta}|\,d\omega}
{\int_{-4.5t}^{5.5t}A^{\mathrm{VED}}_{i,\eta}\,d\omega}.\]
'''
    s+=figure('defect-dense-error','新旧结果与同一加密 VED 参考比较。旧数据仅有0至8号测量位置；连线只引导阅读离散误差值。','Both versions are compared with the same dense VED reference. The previous data cover measured sites 0 through 8 only; lines guide the eye between discrete error values.','两种耦合下，新旧局域谱在0.25t和t展宽时与VED的积分误差。','Integrated differences from VED for the old and new local spectra at broadenings 0.25t and t for both couplings.',1040,410)
    ratios=[r['new_over_old_relative_green_error']['median'] for c in report['couplings'].values() for r in c['records'] if 'new_over_old_relative_green_error' in r]
    s+=p(f'在相同的时间箱和块长下，原有九个位置的 Green 函数相对统计误差，新/旧比值的逐位置中位数为{min(ratios):.2f}–{max(ratios):.2f}。增加统计确实降低了测量噪声，但逆拉普拉斯变换仍会放大误差。图中差异同时包含新增数据和延拓方案变化的影响。',f'With common time bins and block length, the per-site median ratios of new to old relative Green-function statistical errors are {min(ratios):.2f}–{max(ratios):.2f} across the original nine sites. More statistics reduce measurement noise, while inverse Laplace transformation still amplifies errors. The plotted spectral differences reflect both new data and changes to continuation.')
    s+=p('细谱尚未通过定量分辨率验证。新的窄峰图仍应结合候选方案敏感性、bootstrap 范围和 VED 差异阅读；不能仅凭更细的图像确认新激发。','Fine spectral structure has not passed quantitative resolution validation. Read the narrower maps together with representation sensitivity, bootstrap ranges and differences from VED; a finer image alone does not establish a new excitation.')
    s+='<h3>'+tr('延拓、交叉验证与误差','Continuation, cross-validation and uncertainty')+'</h3>\n'
    s+=r'''\[M_0=1,\quad M_1=V_i,\quad M_2=V_i^2+2t^2+g^2,\qquad V_i=-U\delta_{i0}.\]
\[H=H_{e,\mathrm{def}}+\omega_0\sum_i\left(b_i+\frac{g}{\omega_0}n_i\right)^\dagger\left(b_i+\frac{g}{\omega_0}n_i\right)-\frac{g^2}{\omega_0},\qquad
E\ge-\sqrt{4t^2+U^2}-\frac{g^2}{\omega_0}.\]
'''
    s+=p('最大熵延拓使用上述局域谱矩和解析下界；后者由声子项配方及无声子缺陷哈密顿量的最低能量得到，无需 VED 极点先验。主设置为641个谱支撑点，上界10t；另比较321点、上界14t、时间箱、时间窗、块长和协方差截断。旧的321点非负平滑重建也作为候选，保留原来的较松下界。','Maximum-entropy continuation uses these local moments and the analytic lower bound, obtained by completing the phonon square and using the lowest energy of the phonon-free defect Hamiltonian. No VED pole prior is required. The main representation uses 641 support energies with upper support 10t; alternatives test 321 energies, upper support 14t, time bins, time windows, block lengths and covariance cutoffs. The previous 321-energy nonnegative smooth reconstruction is also a candidate, retaining its looser lower bound.')
    s+=p('每折留出两条独立链，四折使用相同验证时间箱和协方差规则。每个方案选择位于最低平均分数一个折间标准误差以内的最大正则化参数，再按所选分数选择方案；所有选择在 VED 对照前冻结。这是模型选择启发式，不是分辨率或置信区间保证。','Each of four folds holds out two independent chains and uses common validation bins and covariance rules. Within each representation, the largest regularization parameter within one fold-standard-error of the minimum mean score is selected; the representation is then chosen by its selected score. All choices are frozen before VED comparison. This is a model-selection heuristic, not a resolution or confidence-interval guarantee.')
    s+=p('每个位置完成64次分块 bootstrap，共2176次；在固定方案内重新选择正则化参数，保存逐点16–84%经验范围。候选方案的谱与分数、多种块长及整链误差均可下载。这些经验范围不包括全部谱网格、支撑上界和方案选择偏差。','Each site has 64 block-bootstrap draws, totaling 2,176. The regularization parameter is reselected within the fixed representation, and pointwise empirical 16–84% ranges are saved. Candidate spectra and scores, multiple block lengths and whole-chain errors are downloadable. These empirical ranges do not include all biases from the energy grid, upper support or representation selection.')
    s+='<h3>'+tr('已知谱峰的恢复检查','Recovery checks on prescribed spectral peaks')+'</h3>\n'
    s+=p('使用满足精确局域谱矩的人工单峰与双峰谱，分别加入由新旧数据测得的相关噪声。两种噪声水平采用相同的641点最大熵表示及留出链选参；双峰间距测试0.3t、0.6t和t，每例16次。该检查不使用 VED 的峰形；固定协方差和有限重复次数限制了它的适用范围。','Prescribed single-peak and doublet spectra obeying the exact local moments are perturbed with correlated noise measured from the old and new data. Both noise levels use the same 641-energy MaxEnt representation and held-out-chain selection. Doublet separations are 0.3t, 0.6t and t, with 16 trials per case. VED line shapes are not used. Fixed covariance and the finite number of trials limit the diagnostic.')
    s+='<div class="table-scroll"><table><caption>'+tr('人工双峰的成功恢复次数 / 64；合并两种耦合与 i=0、8，不能解释为统一物理分辨率','Successful prescribed-doublet recoveries / 64; pooled over both couplings and i=0, 8, not a universal physical resolution')+'</caption><thead><tr><th scope="col">Δω/t</th><th scope="col">'+tr('旧噪声','Old noise')+'</th><th scope="col">'+tr('新噪声','New noise')+'</th></tr></thead><tbody>'
    for separation in [.3,.6,1.]:
        s+=f'<tr><td>{separation:g}</td>'
        for version in ['old','new']:
            rr=[r for r in resolution['records'] if r['noise']==version and r['separation']==separation]
            assert len(rr)==64;s+=f'<td>{sum(r["success"] for r in rr)}/64</td>'
        s+='</tr>'
    s+='</tbody></table></div>\n'
    false=[]
    for version in ['old','new']:
        rr=[r for r in resolution['records'] if r['noise']==version and r['separation']==0];assert len(rr)==64;false.append(sum(r['false_split'] for r in rr))
    s+=p(f'单峰对照被误分裂的次数：旧噪声{false[0]}/64，新噪声{false[1]}/64。成功判据要求窗口内恰有两个足够突出的峰，且峰位分别距设定值不超过0.10t；完整设置与逐次结果见下载文件。',f'The single-peak control is falsely split in {false[0]}/64 old-noise trials and {false[1]}/64 new-noise trials. Doublet success requires exactly two sufficiently prominent peaks in the window, each within 0.10t of its prescribed position. Full settings and per-trial results are available in the download.')
    s+='</section>\n'
    page=repo/'src/defect.html';html=page.read_text();start=html.index('<section id="spectra">');end=html.index('<section id="checks">');html=html[:start]+s+html[end:]
    html=html.replace('局域谱图也已生成，但细峰仍为初步重建；下文给出与 VED 的实际差异。',tr('局域谱已完成四倍统计加密，并扩展至33个显示位置；细峰仍需结合下文的 VED 对照和恢复检查判断。','The local spectra now use four times the production statistics and extend to 33 displayed positions. Fine peaks still require the VED comparison and recovery checks below.'))
    html=html.replace('本轮共148条链','初版共148条链').replace('Slurm 完成记录；热化步不计入生产步','初版 Slurm 完成记录；热化步不计入生产步')
    # Replace outdated spectral convergence only, preserving ground-state text.
    old='局域谱在全部九个位置比较10→12代、半径24→40及 Lanczos 200→400步。\n测试展宽下最大谱差异分别约为0.21%、0.0078%和0.000073%；这是收敛差异，不是无限基底误差的严格上界。'
    cloud=max(np.max(c['reference_cloud_relative_l1']) for c in report['couplings'].values())*100
    space=max(np.max(c['reference_spatial_relative_l1']) for c in report['couplings'].values())*100
    step=max(np.max(c['maximum_steps_relative_l1']) for c in report['couplings'].values())*100
    new=f'加密局域谱在全部17个位置比较10→12代、半径40→64及 Lanczos 200→400步；覆盖四种展宽的最大谱差异依次为{cloud:.3g}%、{space:.3g}%和{step:.3g}%。这是收敛差异，不是无限基底误差的严格上界。'
    html=html.replace(old,tr(new,f'The dense local spectra compare generations 10→12, radius 40→64 and Lanczos steps 200→400 at all 17 sites. The largest spectral differences across all four broadenings are {cloud:.3g}%, {space:.3g}% and {step:.3g}%, respectively. These are convergence differences, not rigorous infinite-basis error bounds.'))
    html=html.replace('C++ 更新核和17项 Python 测试通过。','C++ 更新核和20项 Python 测试通过。')
    extra=p('新增采样为272条链、87.04亿生产步，highmem 作业20882157耗时2分35秒。采样器通过同种子轨迹与 Green 函数逐位一致性检查，所测用例提速约1.64–1.68倍。内存复用、特殊函数复用和精简观测输出减少单链开销；独立链和位置并行运行，每个数值进程使用一个 BLAS 线程。','The new acquisition contains 272 chains and 8.704 billion production steps; highmem job 20882157 took 2 minutes 35 seconds. Fixed-seed checks preserve trajectories and Green-function blocks bit for bit; the tested cases run about 1.64–1.68 times faster. Reused memory and special-function calculations, plus compact measurements, reduce per-chain cost. Independent chains and positions run in parallel, with one BLAS thread per numerical process.')
    html=html.replace('</section>\n<section id="download">',extra+'</section>\n<section id="download">')
    download='<section id="download"><h2>'+tr('数据、代码与复现','Data, code and reproduction')+'</h2>\n'
    download+=p('加密结果、完整选择诊断和恢复检查均已保存；原基态结果与旧版发布包仍可追溯。','Dense results, full selection diagnostics and recovery checks are saved; the original ground-state results and previous release remain available.')
    links=[('defect-dense-summary.json','加密谱误差与收敛 JSON','Dense-spectrum errors and convergence JSON'),('defect-dense-resolution.json','人工谱恢复检查 JSON','Prescribed-spectrum recovery JSON'),('defect-dense-method.md','优化、方法和复现说明','Optimization, methods and reproduction'),('defect-dense-jobs.json','新作业记录','New job records'),('defect-ground.csv','基态表 CSV','Ground-state CSV'),('defect-summary.json','原始基态与初步谱报告','Original ground-state and preliminary-spectrum report')]
    download+='<ul>'+''.join('<li><a href="data/'+name+'">'+tr(zh,en)+'</a></li>' for name,zh,en in links)+'</ul>\n'
    download+='<ul>'
    for lam in [.25,.5]:
        for kind in ['diagmc','ved']:
            name=f'defect-dense-{kind}-lambda{lam:.2f}'
            download+=f'<li>{"DiagMC" if kind=="diagmc" else "VED"}, λ={lam:g}: <a href="data/{name}.npz">NPZ</a> · <a href="data/{name}-eta0.25.csv.gz">η=0.25t CSV (gzip)</a> · <a href="data/{name}-eta0.15.csv.gz">η=0.15t CSV (gzip)</a></li>'
        download+=f'<li>λ={lam:g}: <a href="data/defect-dense-diagnostics-lambda{lam:.2f}.npz">'+tr('候选谱、旧谱与参考截断差异 NPZ','Candidate spectra, previous spectra and reference-cutoff differences NPZ')+'</a></li>'
    download+='</ul>\n'
    download+=p('局域谱 NPZ 的 A 轴为 (eta, site, energy)，形状 (4,33,1601)，eta=[0.15,0.25,0.5,1]。DiagMC 文件含 bootstrap_16、bootstrap_84 和17个独立测量位置 measured_sites。范围为 −4.5t 至5.5t，每列保持绝对谱重；没有逐列归一化或能量平移。','The local-spectrum NPZ array A has axes (eta, site, energy), shape (4,33,1601), and eta=[0.15,0.25,0.5,1]. DiagMC files include bootstrap_16, bootstrap_84 and the 17 independently measured sites in measured_sites. The energy range is −4.5t to 5.5t; every column retains absolute spectral weight, without per-column normalization or energy shifts.')
    download+='<p><a href="https://github.com/rjguo1208/Holstein-model/releases/tag/defect-spectra-dense-20260923">'+tr('本轮源码、原始链、传播子和全部 bootstrap 下载','Download this run’s source, raw chains, Green functions and all bootstrap draws')+'</a> · <a href="data/defect-dense-archives.json">'+tr('发布包大小与 SHA-256','Release sizes and SHA-256')+'</a> · <a href="https://github.com/rjguo1208/Holstein-model/releases/tag/single-defect-20260923">'+tr('保留的初版发布包','Preserved initial release')+'</a>。</p>\n'
    download+=p('所有新 ZIP 解压至同一父目录，共用 holstein-diagmc/ 根目录。脚本拒绝覆盖既有结果；复算时使用新的输出目录。矢量图仅合并色标下颜色完全相同的相邻格块，已验证颜色逐格一致，未平滑谱数据。','Extract all new ZIP files into the same parent directory; they share the holstein-diagmc/ root. Scripts refuse to overwrite existing results, so reruns require fresh output directories. Vector plots merge only adjacent cells with exactly identical colormap colors; cell-by-cell color equality is checked and spectral data are not smoothed.')
    download+='</section>'
    start=html.index('<section id="download">');end=html.index('<footer>',start);html=html[:start]+download+html[end:]
    page.write_text(html)
    # Harvest complete text nodes to preserve catalog strictness. Some retained
    # paragraphs combine the edited sentence with existing text; translate those
    # whole nodes explicitly below rather than weakening localization checks.
    (repo/'src/locales/en.json').write_text(json.dumps(catalog,ensure_ascii=False,indent=2)+'\n')
    source=repo/'research/holstein-defect'
    files=list(ROOT.glob('*.py'))+[ROOT/'DEFECT_DENSE.md',ROOT/'Makefile',ROOT/'requirements.txt']
    files+=[p for d in ['src','tests'] for p in (ROOT/d).rglob('*') if p.is_file() and '__pycache__' not in str(p)]
    files+=[p for p in (ROOT/'scripts').glob('*') if p.is_file() and ('defect' in p.name or p.name=='python.sh')]
    for p in files:
        target=source/p.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
    print('Copied dense data and source; updated shared page and English catalog. Review combined translation nodes before building.')


if __name__=='__main__':main()
