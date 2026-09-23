# Diagrammatic Monte Carlo 与 Holstein 模型

网站：<https://rjguo1208.github.io/Holstein-model/>

中英文理论笔记，讲解 DiagMC 的图空间、Holstein 单极化子的裸图展开、详细平衡、绝对归一化、物理量提取、收敛检查，以及有限电子密度的区别。

[中文](https://rjguo1208.github.io/Holstein-model/) · [English](https://rjguo1208.github.io/Holstein-model/en/)。理论笔记、数值结果、谱函数图和单缺陷模型均有完整英文版；每页右上角的“中文 / English”切换到当前页面的对应语言，页面目录继续使用所选语言。两种语言共享公式、数值、科学图和下载数据。

The theory notes, numerical results, spectral maps and single-defect model are available in Chinese
and English. Use the language links in the top-right corner to switch the current
page. Both versions share equations, numerical tables, figures and downloads.

新增 [完整采样实例](https://rjguo1208.github.io/Holstein-model/#sampling-example)：从零阶图插入一条声子线，再插入得到交叉图，逐段列出电子动量、传播子、顶点因子和图权重，并手算正反提议概率、接受率及详细平衡。四个页面共用页面目录，宽屏固定在左侧，窄屏显示在页首，当前页有明确标记。

新增 [数值结果页](https://rjguo1208.github.io/Holstein-model/results.html)：1D、零温、单电子 bare DiagMC，默认参数 t=ω₀=1、λ=0.25/0.5。基态和有效质量已完成验证，谱函数明确标为初步重建，并列出独立参考误差。原始块数据、代码和报告作为 ZIP 提供下载；所有曲线都有 SVG/PDF 版本。结果页的科学计算与网站排版检查分别记录。

第二轮增加直接质量估计量、真正小 k 的 VED 曲率对照、精确非线性极点、连续矩形 SOM 变体、
重新选择正则化参数的块 bootstrap 与谱分辨率检查。新增采样和后处理现已全部完成，
并保留未通过 η=0.25t 精细谱目标的结果。更新包与首轮数据包分别下载，均不依赖网页脚本。

页面采用白底黑字的单栏排版。公式由 LaTeX 源表达式在构建时排成 HTML + MathML，数学字体随站点托管。四幅费曼图由独立的 LaTeX/TikZ 文件编译成 PDF，再转换为字体已轮廓化的 SVG。阅读时无需 JavaScript、外部字体或 CDN；长公式与图在窄屏中可横向滚动。

新增 [动量–能量谱函数图](https://rjguo1208.github.io/Holstein-model/spectral-map.html)：
默认两种耦合、Γ→X 路径、线性与对数色标，提供 SVG/PDF/PNG、NPZ 和 CSV。
41点 DiagMC 加密计算已在 Kestrel 完成：复用9点、补算32点，共328条链；页面展示实际重建及全部41点与独立 VED 的比较。
细分辨率谱仍未通过稳定性验证，较粗展宽下的整体谱形更稳定。
数值来源、参考收敛和浏览器检查见 [谱图验证记录](docs/spectral-map-verification.md)。
最新作业、精度检查和复现方法见 [Kestrel 41点完成记录](docs/dense-map-kestrel-completed.md)；[17点完成记录](docs/completed-runs-verification.md)保留为历史版本。

[谱函数主图](https://rjguo1208.github.io/Holstein-model/spectral-map.html#diagmc-map)的两种耦合、全部41点均已使用改进结果，与相同坐标、展宽和绝对色标的 VED 并排展示。
每点4倍采样和128次 bootstrap；本轮补齐31点、新增496条链，完整网格共656条链、314.88亿生产步和10496次 bootstrap。
η=0.25t 的全41点误差中位数由 16.3%／23.3% 变为 10.2%／10.1%；70/82组误差下降。精细窄峰仍未普遍通过验证。
见 [全部41点完成报告](docs/full-grid-refinement-completed.md)、[科学源码与复现方法](research/spectral-fullgrid/README.md) 和 [完整复现数据](https://github.com/rjguo1208/Holstein-model/releases/tag/full-grid-refinement-20260919)。

前轮10点试验、864次模拟谱恢复与原始谱图保留在页面的历史记录中；
[前轮报告](docs/spectral-refinement-completed.md) 与 [前轮源码](research/spectral-refinement/README.md) 保持可访问。

新增 [单缺陷模型](https://rjguo1208.github.io/Holstein-model/defect.html)（[English](https://rjguo1208.github.io/Holstein-model/en/defect.html)）：
只改变一个格点的电子势能，实现连续时间实空间 bare DiagMC，计算两种耦合和三种吸引强度。
六组基态能量与独立 VED 的最大差异约为1.11个统计标准误差。局域谱图并列展示 DiagMC 与 VED，
细峰明确标为初步重建；源码和148条链的数据见 [方法](research/holstein-defect/README.md)、
[验证记录](docs/defect-verification.md) 与 [复现数据](https://github.com/rjguo1208/Holstein-model/releases/tag/single-defect-20260923)。

The new single-defect calculation expands hopping and electron–phonon vertices
in continuous imaginary time. Six bound-state energies agree with independent
VED within about 1.11 estimated statistical standard errors. Position–energy
local spectral maps include explicit reference comparisons and retain their
preliminary status. All 148 chains and reproducible source are available.

## 本地预览

已生成的页面与全部资源保存在 `site/`：

```bash
cd site
python3 -m http.server 8000 --bind 127.0.0.1
```

打开 <http://localhost:8000>。

## 修改正文或公式

需要 Node.js 22 和 Python 3。修改 `src/*.html` 中的正文和 LaTeX，或修改 `src/style.css`，然后运行：

```bash
npm ci --ignore-scripts
npm run build
npm run check
```

使用 `\(...\)` 写行内公式，`\[...\]` 写独立公式。构建采用固定版本的 [KaTeX](https://katex.org/docs/api.html)，遇到不支持的公式语法即报错。数学表达式中的小于号使用 `\lt`，不使用 HTML 实体。`site/*.html`、`site/en/*.html` 和 `site/assets/` 均为已提交的发布产物。

## 中英文同步维护 / Bilingual updates

以后每次更新网站，必须同时更新中文和英文；这项用户要求也记录在 [AGENTS.md](AGENTS.md)。
共同结构与中文内容保存在 `src/*.html`，经过人工翻译的英文保存在 `src/locales/en.json`。
键是合并空白后的中文原文，值是对应英文。查看当前所有待翻译文本：

```bash
python3 scripts/localize.py --extract
```

修改中文后，在同一提交中更新对应英文；标题、段落、表格说明、图注、下载说明、元数据、替代文本以及公式中的文字说明均须同步。构建会拒绝缺失、空白、仍含中文或已过时的翻译，并检查两种语言的数学表达式一致。
`npm run check` 检查所有页面的双语配对、当前语言、页面目录、本地链接、公式、锚点、数值表和共用资源。
科学图使用英文或数学标注，两种语言直接复用同一数据与图。提交重建后的中文和英文产物，再一起部署。

For every website update, edit the shared Chinese sources and their English
translations together. The build fails on missing or obsolete translation keys
and changed mathematics. Run both build and checks, review desktop/mobile
language switching, and commit and deploy both generated versions together.

## 修改费曼图

`site/figures/` 中的四份 `.tex` 均可独立编译，使用标准 LaTeX、AMS Math 和 TikZ。每幅图在网页中同时提供 SVG、PDF 和 LaTeX 源文件。

安装 [Tectonic](https://tectonic-typesetting.github.io/en-US/install.html) 和 Poppler（提供 `pdftocairo`）后：

```bash
npm run figures
npm run check
```

如果 Tectonic 不在 PATH，可设置 `TECTONIC=/path/to/tectonic`。构建脚本把 LaTeX 源文件的 SHA-256 写入 SVG；站点检查会拒绝源文件已改变但图未重新编译的情况。

## 发布

GitHub Pages 的发布来源为 GitHub Actions。推送到 `main` 会运行检查；发布仍由 `Deploy GitHub Pages` 工作流手动触发：

```bash
gh workflow run pages.yml --ref main
```

检查和发布流程都会重新渲染公式、验证本地链接与图文件，并确认发布产物与源文件一致。TikZ 产物已提交，因此普通页面构建无需下载 TeX 发行版。工作流只上传 `site/`，所有资源路径兼容 `/Holstein-model/` 子路径。

## 文件结构

```text
src/index.html                正文及 LaTeX 公式
src/results.html              数值结果、误差与初步谱函数
src/spectral-map.html         谱函数图、VED 比较与完整历史记录
src/locales/en.json           全站英文翻译；键为中文原文
src/style.css                 简洁排版与打印样式
site/index.html               预先渲染的静态网页
site/en/                      完整英文版三个页面
site/assets/katex/            数学样式、字体与许可证
site/figures/                 费曼图的 LaTeX、PDF 与 SVG
site/results/                 科学图的 SVG 与 PDF
site/data/                    结果摘要、CSV、代码和原始数据包
scripts/build.mjs             公式渲染与静态资源复制
scripts/localize.py           英文翻译、共用链接与公式一致性检查
scripts/build_figures.py      TikZ → PDF → SVG
scripts/check_site.py         数学、链接、字体与图来源检查
.github/workflows/            自动检查及 Pages 发布
AGENTS.md                     后续更新必须同时维护中英文
```

KaTeX 的许可证随发布资源保留在 `site/assets/katex/LICENSE`。科学文献链接见网页末尾。

本版的浏览器与排版检查见 [验证记录](docs/verification.md)。

采样手算实例与全站页面目录的检查见 [实例及导航验证记录](docs/sampling-example-verification.md)。

数值结果与新增页面的检查见 [结果验证记录](docs/results-verification.md)。

第二轮的完成范围、排队作业、数值限制与显示检查见 [第二轮验证记录](docs/results-v2-verification.md)。

中英文翻译、同步构建及桌面／手机语言切换的检查见 [双语网站验证记录](docs/bilingual-website-verification.md)。
