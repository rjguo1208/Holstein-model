"""Publish completed momentum spectra, shared bilingual text and compact previews."""
from __future__ import annotations
import argparse
import importlib.util
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
TAG = 'defect-momentum-20260923'
RELEASE = 'https://github.com/rjguo1208/Holstein-model/releases/tag/'+TAG
DOWNLOAD = 'https://github.com/rjguo1208/Holstein-model/releases/download/'+TAG+'/'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=ROOT.parent/'Holstein-model')
    args = parser.parse_args()
    ved = ROOT/'results/defect_momentum_ved_report_02'
    mc = ROOT/'results/defect_momentum_mc_report_01'
    vs = json.loads((ved/'summary.json').read_text())
    ms = json.loads((mc/'summary.json').read_text())
    assert vs['complete'] and ms['complete']
    verification = json.loads((ROOT/'results/defect_momentum_release_02/verification.json').read_text())
    assert verification['passed']
    catalog_path = args.repo/'src/locales/en.json'
    catalog = json.loads(catalog_path.read_text())
    def tr(zh, en):
        catalog[' '.join(zh.split())] = en
        return zh
    parts = []
    def prose(zh, en):
        parts.append('<p>'+tr(zh, en)+'</p>')
    def heading(tag, zh, en, anchor=None):
        parts.append('<'+tag+(f' id="{anchor}"' if anchor else '')+'>'+tr(zh, en)+'</'+tag+'>')
    def figure(folder, name, height, zh, en, alt_zh, alt_en):
        shutil.copyfile(folder/(name+'.png'), args.repo/'site/results'/(name+'.png'))
        parts.append(f'<figure><div class="figure-scroll data"><img src="results/{name}.png" width="1040" height="{height}" alt="'+tr(alt_zh, alt_en)+'"></div><figcaption>'+tr(zh, en)+'</figcaption><p class="figure-links">'+
                     ''.join(f'<a href="{DOWNLOAD}{name}.{ext}">{ext.upper()}</a>' for ext in ['svg', 'pdf', 'png'])+'</p></figure>')
    parts.append('<section id="momentum">')
    heading('h2', '缺陷体系动量谱：201个 k 点', 'Defect momentum spectra: 201 k points')
    prose(r'新增结果使用相同参数 \(t=\omega_0=U=1\)、\(\lambda=0.25,0.5\)。先完成实空间 VED，再独立重新采样非局域 bare DiagMC。横轴为 \(k/\pi\)，纵轴为相对真空的能量，颜色为谱强度；白色虚线为裸宿主余弦色散 \(\varepsilon_k=-2t\cos k\)。',
          r'The new results use the same parameters \(t=\omega_0=U=1\) and \(\lambda=0.25,0.5\). Real-space VED was completed first, followed by independent fresh sampling of nonlocal bare DiagMC. The horizontal axis is \(k/\pi\), the vertical axis is energy relative to the vacuum, and color gives spectral intensity. The white dashed curve is the bare-host cosine dispersion \(\varepsilon_k=-2t\cos k\).')
    heading('h3', '相干观测窗口与动量分辨率', 'Coherent observation window and momentum resolution')
    parts.append(r'\[c^\dagger_{k,L}=\frac{1}{\sqrt L}\sum_{j=-16}^{16}e^{ikj}c_j^\dagger,\qquad L=33.\]')
    parts.append(r'\[A_{L,\eta}(k,\omega)=-\frac{1}{\pi L}\operatorname{Im}\sum_{i,j=-16}^{16}e^{ik(j-i)}\langle\mathrm{vac}|c_i(\omega+i\eta-H)^{-1}c_j^\dagger|\mathrm{vac}\rangle.\]')
    prose('单缺陷破坏平移对称性，k 在这里标记入射探测态。主图采用以缺陷为中心的33格点相干窗口，保留所有非局域传播贡献；电子在传播过程中可以离开窗口。它与33格点封闭链、周期缺陷阵列、对局域谱直接做傅里叶变换的结果不同。',
          'A single defect breaks translation symmetry, so k labels the incident probe state here. The main plots use a coherent 33-site window centered on the defect and retain all nonlocal propagation contributions. The electron can leave the window during propagation. This differs from an isolated 33-site chain, a periodic defect array, or a Fourier transform of local spectral intensities.')
    prose(r'\(k\in[0,\pi]\) 共201点，\(\Delta k=\pi/200\)；能量范围 \([-4.5t,5.5t]\) 共1601点，\(\Delta\omega=0.00625t\)。每个 k 都直接由传播矩阵或端点位移统计求和，没有动量插值。相邻点共享实空间信息与统计误差；33格点矩形窗口的动量宽度约为 \(2\pi/33\)，不能仅凭网格更密宣称物理分辨率提高。',
          r'There are 201 points over \(k\in[0,\pi]\), with \(\Delta k=\pi/200\), and 1601 energies over \([-4.5t,5.5t]\), with \(\Delta\omega=0.00625t\). Each k is evaluated directly from the propagation matrix or measured endpoint displacements without momentum interpolation. Neighboring points share real-space information and statistical errors. The 33-site rectangular window has a momentum width of order \(2\pi/33\); a denser grid alone does not establish higher physical resolution.')
    prose('VED 另计算65格点观测窗口。扩大观测窗口会改变物理观测量，并稀释单个缺陷的相对谱重；这与扩大数值边界以检查收敛是两件事。全部图保留绝对谱重，未逐列归一化或平移能量。',
          'VED also evaluates a 65-site observation window. Enlarging this window changes the observable and dilutes the relative spectral weight of one defect; it is distinct from moving the numerical boundary to check convergence. Every plot retains absolute spectral weights without normalizing individual columns or shifting energies.')
    heading('h3', 'VED：含缺陷与干净体系', 'VED: defect and clean systems', 'momentum-ved')
    figure(ved, 'defect-momentum-ved-eta0.25', 880,
           r'\(\eta=0.25t\)。上排为单缺陷，下排为相同观测窗口的干净体系；两列对应两种耦合，四图共用色标。VED 使用声子云12代、电子半径96、Lanczos 800步。',
           r'\(\eta=0.25t\). The top row shows one defect and the bottom row shows the clean system with the same observation window. Columns show the two couplings and all four panels share a color scale. VED uses cloud generation 12, electron radius 96 and 800 Lanczos steps.',
           '201个动量点的VED谱图；比较单缺陷与干净体系，白色虚线为弯曲的裸电子余弦色散。',
           'VED spectra at 201 momenta comparing a single defect with the clean system; the white dashed curve is the bare-electron cosine dispersion.')
    figure(ved, 'defect-momentum-ved-eta0.15-difference', 480,
           r'\(\eta=0.15t\) 的缺陷减干净体系差值：红色表示增重，蓝色表示减重，零点为白色。低能束缚态的能量固定，但探测谱重随 k 改变；水平束缚态信号与余弦裸带参照可以同时存在。',
           r'Defect-minus-clean difference at \(\eta=0.15t\): red indicates increased weight, blue indicates decreased weight, and white denotes zero. A localized bound state has a fixed energy but its probe weight varies with k; a horizontal bound-state signal can coexist with the cosine bare-band reference.',
           '缺陷与干净体系的动量谱差值；红蓝色分别显示谱重增加和减少，保留余弦裸带虚线。',
           'Momentum-spectrum difference between defect and clean systems; red and blue show increased and decreased weight, with a dashed cosine bare-band guide.')
    prose(r'独立检查覆盖全部动量：声子云10→12代、电子半径64→96、递推400→800步；另在每组参数的五个动量直接注入余弦/正弦态，检查非局域矩阵重建。声子截断造成的最大整谱相对 L1 差异，在 \(\eta/t=0.25,0.15,0.1\) 时分别为0.846%、3.021%、6.749%。更窄展宽图显示截断敏感性，不能当作全部细结构已收敛的证明。',
          r'Independent checks cover every momentum: cloud generations 10→12, electron radius 64→96 and recursion length 400→800. Direct cosine/sine injection at five momenta per parameter pair also checks the nonlocal reconstruction. The largest full-spectrum relative L1 differences from the cloud cutoff are 0.846%, 3.021% and 6.749% at \(\eta/t=0.25,0.15,0.1\), respectively. The narrower maps reveal cutoff sensitivity and do not establish convergence of all fine structure.')
    parts.append('<p>'+tr('较窄展宽的完整 VED 图：', 'Complete VED maps with narrower broadening: ')+
                 f'<a href="{DOWNLOAD}defect-momentum-ved-eta0.15.pdf">η=0.15t PDF</a> · '+
                 f'<a href="{DOWNLOAD}defect-momentum-ved-eta0.1.pdf">η=0.1t PDF</a>。</p>')
    heading('h3', '重新采样的非局域 bare DiagMC', 'Freshly sampled nonlocal bare DiagMC', 'momentum-diagmc')
    prose('新的采样器允许电子起点和终点不同，增加单次跳跃的插入/删除与起点平移热浴更新，直接记录端点位移。两种耦合各8条独立链，每条1.28亿生产步，共20.48亿步；虚时间箱192个，最大虚时12/t。采样顺序在 VED 完成之后，16个进程在 highmem 并行，所有阶数上限触发次数为零。',
          'The new sampler allows different initial and final electron sites, adds single-hop insertion/removal and an origin-translation heat bath, and directly records endpoint displacement. Each coupling has eight independent chains with 128 million production steps per chain, totaling 2.048 billion steps, with 192 time bins and maximum imaginary time 12/t. Sampling followed completion of VED, using 16 parallel processes on highmem. No order cap was reached.')
    prose('正式采样前，自由电子、无声子缺陷和原子极限通过24条链的解析检查，最大偏差为2.21个统计标准误差。编译测试检查了详细平衡、奇数跳跃、声子局域性、起点热浴和中途越出观测窗口。已有闭合路径局域谱数据没有被重新标记成动量谱。',
          'Before production, 24 chains checked analytic free-electron, phonon-free impurity and atomic limits; the largest discrepancy was 2.21 statistical standard errors. Compiled tests check detailed balance, odd hopping orders, phonon locality, the origin heat bath and excursions outside the observation window. Earlier closed-path local spectra were not relabeled as momentum spectra.')
    parts.append(r'\[M_0=1,\qquad M_1=-2t\left(1-\frac1L\right)\cos k-\frac UL,\]')
    parts.append(r'\[M_2=2t^2+2t^2\left(1-\frac2L\right)\cos2k+\frac{4tU}{L}\cos k+\frac{U^2}{L}+g^2.\]')
    prose('延拓使用未展宽谱的精确窗口谱矩与解析能量下界，不输入 VED 极点或谱形。带符号 Green 函数和完整协方差均保留，不裁剪噪声负值。每折留出两条链进行四折选择，比较六种时间、谱网格、块长或协方差设置；各动量完成32次分块 bootstrap，共12864次。16–84%经验范围以所选表示为条件，不包括全部系统偏差，也不是相邻动量的独立误差。',
          'Continuation uses the exact window moments of the unbroadened spectrum and an analytic energy lower bound, with no VED pole or spectral-shape input. Signed Green functions and full covariance are retained without clipping noisy negative values. Four folds each hold out two chains, comparing six time, spectral-grid, block-length or covariance settings. Each momentum has 32 block-bootstrap repeats, totaling 12864. The 16–84% empirical ranges are conditional on the selected representation, exclude some systematic biases and are not independent errors at neighboring momenta.')
    figure(mc, 'defect-momentum-diagmc-eta0.25', 880,
           r'\(\eta=0.25t\)。左列为全新非局域 DiagMC 重建，右列为相同窗口和展宽下的 VED；每行对应一种耦合。颜色保持绝对强度。网格加密未被当作细峰分辨率验证。',
           r'\(\eta=0.25t\). The left column shows fresh nonlocal DiagMC reconstruction and the right column shows VED with the same window and broadening. Rows show the two couplings. Colors retain absolute intensity; grid refinement is not treated as validation of fine-peak resolution.',
           '201个动量点的非局域DiagMC与VED对照谱图，展宽0.25t，虚线为裸电子余弦色散。',
           'Nonlocal DiagMC versus VED spectra at 201 momenta and broadening 0.25t, with a dashed bare-electron cosine dispersion.')
    prose('DiagMC 图中个别竖直接缝与逐动量的离散选参有关，不能解释为物理能带断裂。这里保留原始重建结果，没有沿动量方向平滑这些变化。',
          'Some vertical seams in the DiagMC maps are associated with discrete parameter selection at individual momenta and should not be interpreted as physical band discontinuities. The original reconstructions are retained without smoothing these changes along momentum.')
    figure(mc, 'defect-momentum-diagmc-eta1', 880,
           r'\(\eta=t\) 的较宽展宽对照；窗口、动量点、能量零点与上图相同。窄峰仍需结合候选方案、bootstrap 和独立 VED 差异判断。',
           r'Comparison at the broader \(\eta=t\), with the same window, momenta and energy zero as above. Narrow peaks still require assessment using candidate sensitivity, bootstrap ranges and differences from independent VED.',
           '展宽t时的非局域DiagMC与VED动量谱对照。',
           'Comparison of nonlocal DiagMC and VED momentum spectra at broadening t.')
    heading('h3', '动量谱数据与复现', 'Momentum-spectrum data and reproduction')
    parts.append('<div class="table-scroll"><table><caption>'+tr('全部201个动量：DiagMC 与 VED 的整谱相对 L1 差异范围', 'All 201 momenta: range of full-spectrum relative L1 differences between DiagMC and VED')+'</caption><thead><tr><th scope="col">λ</th><th scope="col">η=0.15t</th><th scope="col">η=0.25t</th><th scope="col">η=t</th></tr></thead><tbody>')
    for lam in [.25, .5]:
        r = ms['couplings'][str(lam)]
        cells = [f'{100*r["minimum_relative_l1_by_eta"][i]:.1f}%–{100*r["maximum_relative_l1_by_eta"][i]:.1f}%' for i in [1, 2, 4]]
        parts.append('<tr><td>'+str(lam)+'</td>'+''.join('<td>'+x+'</td>' for x in cells)+'</tr>')
    parts.append('</tbody></table></div>')
    prose(r'在 \(\tau\leq6/t\) 内，未延拓 Green 函数与独立 VED 的最大差异为2.47个保守统计标准误差；这里取多种块长与整链估计中较大的误差。\(\eta=t\) 时两种耦合的整谱差异均小于1.75%。这些检查支持采样与宽谱形的一致性，不能证明细小峰已经分辨。',
          r'For \(\tau\leq6/t\), the largest discrepancy between uncontinued Green functions and independent VED is 2.47 conservative statistical standard errors, using the largest error from multiple block lengths and whole-chain estimates. At \(\eta=t\), full-spectrum differences remain below 1.75% for both couplings. These checks support consistency of sampling and broad spectral shapes, but do not establish resolution of fine peaks.')
    prose('误差在展示能量窗内积分；截断、统计噪声和延拓偏差都可能贡献差异。新的 DiagMC 窄峰尚未通过定量分辨率验证。完整数据保留未延拓的 Green 函数、不同块长误差、候选方案与 bootstrap 权重。',
          'Differences are integrated over the displayed energy window and can include cutoff effects, sampling noise and continuation bias. Narrow peaks in the new DiagMC spectra have not passed quantitative resolution validation. Complete data retain uncontinued Green functions, errors at different block lengths, candidate representations and bootstrap weights.')
    parts.append('<ul>')
    for lam in [.25, .5]:
        for U in [0, 1]:
            name = f'lambda{lam:.2f}_U{U}_Nh12_R96_W32.npz'
            parts.append(f'<li>VED, λ={lam:g}, U/t={U}: <a href="{DOWNLOAD}{name}">NPZ</a></li>')
        name = f'defect-momentum-diagmc-lambda{lam:.2f}.npz'
        parts.append(f'<li>DiagMC, λ={lam:g}, U/t=1: <a href="{DOWNLOAD}{name}">NPZ</a></li>')
    for file, zh, en in [('defect-momentum-summary.json', '动量谱检查与误差 JSON', 'Momentum-spectrum checks and errors JSON'),
                          ('defect-momentum-method.md', '观测量、优化与复现说明', 'Observable, optimization and reproduction'),
                          ('defect-momentum-jobs.json', 'highmem 作业记录', 'Highmem job records'),
                          ('defect-momentum-archives.json', '发布文件大小与 SHA-256', 'Release file sizes and SHA-256')]:
        parts.append(f'<li><a href="data/{file}">'+tr(zh, en)+'</a></li>')
    parts.append(f'<li><a href="{RELEASE}">'+tr('完整源码、原始链、VED递推与延拓结果发布包', 'Complete release: source, raw chains, VED recursions and continuation results')+'</a></li></ul>')
    prose('VED NPZ 的 A_W16、A_W32 分别对应33和65格点窗口；DiagMC NPZ 的 A、A_ved、bootstrap_16、bootstrap_84 均按 (eta, k, energy) 排列，形状为 (5,201,1601)，eta=[0.1,0.15,0.25,0.5,1]。网页使用清晰的 PNG 预览；链接提供原始矢量 SVG/PDF 与数值数组。',
          'VED NPZ arrays A_W16 and A_W32 correspond to 33-site and 65-site windows. DiagMC NPZ arrays A, A_ved, bootstrap_16 and bootstrap_84 use axes (eta, k, energy), shape (5,201,1601), and eta=[0.1,0.15,0.25,0.5,1]. Pages use clear PNG previews and link to original SVG/PDF vectors and numerical arrays.')
    parts.append('</section>')
    page = args.repo/'src/defect.html'
    html = page.read_text()
    if '<section id="momentum">' in html:
        raise ValueError('Momentum section already exists; review before replacing')
    html = html.replace('<section id="spectra">', '\n'.join(parts)+'\n<section id="spectra">', 1)
    toc = tr('动量–能量谱图', 'Momentum–energy spectra')
    html = html.replace('<li><a href="#spectra">', f'<li><a href="#momentum">{toc}</a></li><li><a href="#spectra">', 1)
    old_title = '单缺陷 Holstein：DiagMC 结果与局域谱'
    title = tr('单缺陷 Holstein：局域谱与动量谱', 'Single-defect Holstein: local and momentum spectra')
    html = html.replace(old_title, title)
    old_desc = '1D零温单电子Holstein模型的单点吸引势：实空间bare DiagMC、独立VED验证、束缚态与初步局域谱函数图。'
    desc = tr('1D零温单电子Holstein模型：单缺陷局域谱、201点VED动量谱、全新非局域bare DiagMC及误差检查。',
              '1D zero-temperature one-electron Holstein model: single-defect local spectra, 201-point VED momentum spectra, fresh nonlocal bare DiagMC and error checks.')
    html = html.replace(old_desc, desc)
    page.write_text(html)
    spec = importlib.util.spec_from_file_location('site_localize', args.repo/'scripts/localize.py')
    localize = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(localize)
    _, strings = localize.build()
    keys = {key for values in strings.values() for key in values}
    catalog = {k:v for k,v in catalog.items() if k in keys}
    missing = keys-catalog.keys()
    if missing:
        raise ValueError('Missing translations: '+repr(missing))
    catalog_path.write_text(json.dumps(catalog, ensure_ascii=False, indent=2)+'\n')
    compact = dict(ved=vs, diagmc={k:v for k,v in ms.items() if k != 'couplings'})
    compact['diagmc']['couplings'] = {k:{kk:vv for kk,vv in v.items() if kk != 'records'} for k,v in ms['couplings'].items()}
    (args.repo/'site/data/defect-momentum-summary.json').write_text(json.dumps(compact, indent=2)+'\n')
    shutil.copyfile(ROOT/'DEFECT_MOMENTUM.md', args.repo/'site/data/defect-momentum-method.md')
    for name in ['defect-momentum-jobs.json', 'defect-momentum-archives.json']:
        shutil.copyfile(ROOT/'results/defect_momentum_release_02'/name, args.repo/'site/data'/name)
    destination = args.repo/'research/holstein-defect'
    sources = [ROOT/'defect_momentum.py', ROOT/'defect_momentum_mc.py', ROOT/'defect_momentum_fit.py',
               ROOT/'defect_entropy.py', ROOT/'DEFECT_MOMENTUM.md', ROOT/'Makefile']
    sources += list((ROOT/'src').glob('defect_momentum*'))+list((ROOT/'tests').glob('*defect_momentum*'))
    sources += list((ROOT/'scripts').glob('*defect_momentum*'))
    for file in sources:
        if not file.is_file():
            continue
        target = destination/file.relative_to(ROOT)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(file, target)
    readme = args.repo/'README.md'
    addition = '''

新增 [缺陷动量谱](https://rjguo1208.github.io/Holstein-model/defect.html#momentum)
（[English](https://rjguo1208.github.io/Holstein-model/en/defect.html#momentum)）：
201个动量点、1601个能量点，采用以缺陷为中心的33格点相干观测窗口。
先完成 VED，再以开放电子端点重新采样20.48亿步 bare DiagMC；两种方法均叠加裸宿主余弦色散，
并提供干净体系对照、缺陷差值图及12864次条件 bootstrap。窄峰仍未通过定量分辨率验证。
方法与复现见 [DEFECT_MOMENTUM.md](research/holstein-defect/DEFECT_MOMENTUM.md)，
完整数据见 [动量谱发布包](https://github.com/rjguo1208/Holstein-model/releases/tag/defect-momentum-20260923)。

The new defect momentum spectra use 201 momenta and 1601 energies for a coherent
33-site observation window centered on one defect. VED was completed first,
followed by 2.048 billion fresh bare DiagMC steps with open electron endpoints.
Both methods include the bare-host cosine guide, with clean-system and defect-difference
comparisons and 12864 conditional bootstrap repeats. Narrow peaks are not certified
as resolved. Full definitions, cutoff checks, source and data are linked above.
'''
    if 'defect-momentum-20260923' not in readme.read_text():
        readme.write_text(readme.read_text().rstrip()+addition)
    print('Published both language sources, four compact preview plots and reproducible research code.')


if __name__ == '__main__':
    main()
