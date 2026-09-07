import { broaden, subtractLowest, relativeL1, selectLR } from './spectra.mjs';

const $ = id => document.getElementById(id);
const state = { case: 'B', k: 0, eta: 0.05, samples: 262144, rcond: 0.001, rest: false };
const cache = new Map();
let generation = 0;
let downloadable = null;
const colors = { ved: '#1b655d', lr: '#c36a43', grid: '#e4e6de', label: '#737b73' };
const kLabel = k => k === 0 ? '0' : k === 0.5 ? 'π/2' : 'π';
const percentage = value => value === null ? '未定义' : `${(value * 100).toFixed(2)}%`;
const escapeHTML = text => String(text).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));

function ticks(lo, hi, count = 5) {
  return Array.from({ length: count }, (_, i) => lo + (hi - lo) * i / (count - 1));
}

function formatAxis(value) {
  if (Math.abs(value) < 1e-10) return '0';
  if (Math.abs(value) < .01) return value.toExponential(1);
  return Number(value.toFixed(2)).toString();
}

function drawPlot(element, series, options) {
  const width = options.small ? 420 : 850;
  const height = options.small ? 225 : 340;
  const margin = { top: 18, right: 18, bottom: 41, left: options.small ? 54 : 55 };
  const iw = width - margin.left - margin.right;
  const ih = height - margin.top - margin.bottom;
  const xMin = options.xMin;
  const xMax = options.xMax;
  const allValues = series.flatMap(s => Array.from(s.y));
  let yMin = Math.min(...allValues);
  let yMax = Math.max(...allValues);
  if (options.fromZero) yMin = Math.min(0, yMin);
  const pad = Math.max((yMax - yMin) * .12, options.small ? .0001 : .002);
  yMax += pad;
  if (!options.fromZero || yMin < 0) yMin -= pad;
  if (yMax === yMin) yMax += .01;
  const px = x => margin.left + (x - xMin) * iw / (xMax - xMin);
  const py = y => margin.top + (yMax - y) * ih / (yMax - yMin);
  const id = element.id;
  let svg = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" role="img" aria-labelledby="${id}-title ${id}-desc"><title id="${id}-title">${escapeHTML(options.title)}</title><desc id="${id}-desc">${escapeHTML(options.description)}</desc><defs><clipPath id="${id}-clip"><rect x="${margin.left}" y="${margin.top}" width="${iw}" height="${ih}"/></clipPath></defs>`;
  for (const tick of ticks(yMin, yMax)) {
    svg += `<line x1="${margin.left}" x2="${width - margin.right}" y1="${py(tick)}" y2="${py(tick)}" stroke="${colors.grid}" stroke-width="1"/><text x="${margin.left - 10}" y="${py(tick) + 4}" text-anchor="end" fill="${colors.label}" font-family="sans-serif" font-size="11">${formatAxis(tick)}</text>`;
  }
  const xticks = options.small ? [0, .5, 1] : ticks(xMin, xMax, 7);
  for (const tick of xticks) {
    svg += `<line x1="${px(tick)}" x2="${px(tick)}" y1="${height - margin.bottom}" y2="${height - margin.bottom + 4}" stroke="${colors.label}"/><text x="${px(tick)}" y="${height - margin.bottom + 20}" text-anchor="middle" fill="${colors.label}" font-family="sans-serif" font-size="11">${options.small ? kLabel(tick) : formatAxis(tick)}</text>`;
  }
  svg += `<text x="${margin.left}" y="10" fill="${colors.label}" font-family="sans-serif" font-size="10">${escapeHTML(options.yLabel)}</text><text x="${width - margin.right}" y="${height - 4}" text-anchor="end" fill="${colors.label}" font-family="sans-serif" font-size="11">${options.small ? 'k' : 'ω / t'}</text>`;
  for (const s of series) {
    const path = s.x.map((x, i) => `${i === 0 ? 'M' : 'L'}${px(x).toFixed(3)},${py(s.y[i]).toFixed(3)}`).join(' ');
    if (!options.small && !s.dashed) {
      svg += `<path d="${path} L${px(s.x.at(-1))},${py(0)} L${px(s.x[0])},${py(0)} Z" fill="${s.color}" opacity=".055" clip-path="url(#${id}-clip)"/>`;
    }
    svg += `<path d="${path}" fill="none" stroke="${s.color}" stroke-width="${options.small ? 1.8 : 2}" ${s.dashed ? 'stroke-dasharray="6 4"' : ''} stroke-linejoin="round" clip-path="url(#${id}-clip)"/>`;
    if (options.small) {
      for (let i = 0; i < s.x.length; i++) {
        svg += s.dashed
          ? `<rect x="${px(s.x[i]) - 3}" y="${py(s.y[i]) - 3}" width="6" height="6" fill="#fffefa" stroke="${s.color}" stroke-width="1.6"><title>${s.name}: k=${kLabel(s.x[i])}, ${s.y[i].toPrecision(9)}</title></rect>`
          : `<circle cx="${px(s.x[i])}" cy="${py(s.y[i])}" r="3.1" fill="${s.color}"><title>${s.name}: k=${kLabel(s.x[i])}, ${s.y[i].toPrecision(9)}</title></circle>`;
      }
    }
  }
  element.innerHTML = svg + '</svg>';
}

async function getCase(name) {
  if (!cache.has(name)) {
    const response = await fetch(`data/${name}.json`);
    if (!response.ok) throw new Error(`数据请求失败（HTTP ${response.status}）`);
    const data = await response.json();
    if (data.case !== name || data.momenta.length !== 3) throw new Error('数据格式与参数不匹配');
    cache.set(name, data);
  }
  return cache.get(name);
}

async function render() {
  const version = ++generation;
  $('results-panel').setAttribute('aria-busy', 'true');
  $('download-spectrum').disabled = true;
  $('raw-data-link').removeAttribute('href');
  $('raw-data-link').setAttribute('aria-disabled', 'true');
  $('case-label').textContent = `CASE ${state.case} / 数据加载中`;
  $('load-status').textContent = '正在读取并展宽保存的谱…';
  downloadable = null;
  try {
    const data = await getCase(state.case);
    if (version !== generation) return;
    const point = data.momenta.find(p => p.k_pi === state.k);
    const lr = selectLR(point, state.samples, state.rcond);
    const refValues = broaden(point.reference.poles, data.omega, state.eta);
    const lrValues = broaden(lr.poles, data.omega, state.eta);
    const refRest = subtractLowest(refValues, data.omega, point.reference.E, point.reference.Z, state.eta);
    const lrRest = subtractLowest(lrValues, data.omega, lr.E, lr.Z, state.eta);
    const plottedRef = state.rest ? refRest : refValues;
    const plottedLR = state.rest ? lrRest : lrValues;
    const error = relativeL1(plottedRef, plottedLR, data.omega);
    $('case-label').textContent = `CASE ${state.case} / ω₀/t = ${data.model.omega0} / λ = ${data.model.lambda_} / k = ${kLabel(state.k)}`;
    $('spectrum-title').textContent = state.rest ? '扣除最低极点的余谱' : '电子加谱 A(k, ω)';
    $('spectrum-caption').textContent = `η/t = ${state.eta}；窗口 ω/t ∈ [${data.omega[0].toFixed(3)}, ${data.omega.at(-1).toFixed(3)}]。${state.rest ? '余谱仍可包含其他束缚态。' : '窗口外的 Lorentz 尾部未重新归一化。'}`;
    const series = [
      { x: data.omega, y: plottedRef, color: colors.ved, name: 'VED' },
      { x: data.omega, y: plottedLR, color: colors.lr, dashed: true, name: 'LR-VMC' },
    ];
    drawPlot($('spectrum-plot'), series, { xMin: data.omega[0], xMax: data.omega.at(-1), fromZero: true,
      title: `${state.case} 组 k=${kLabel(state.k)}，η=${state.eta} 的${state.rest ? '余谱' : '谱函数'}`,
      description: '绿色实线为 VED，橙色虚线为随机 LR-VMC，保留所有原始谱权重。可下载数值 CSV。',
      yLabel: state.rest ? 'Arest × t' : 'A(k, ω) × t' });
    $('energy-ref').textContent = point.reference.E.toFixed(6);
    $('energy-lr').textContent = lr.E.toFixed(6);
    $('z-ref').textContent = point.reference.Z.toFixed(6);
    $('z-lr').textContent = lr.Z.toFixed(6);
    $('l1-label').textContent = state.rest ? '余谱 L¹ 差异' : '全谱 L¹ 差异';
    $('spectral-error').textContent = percentage(error);
    $('spectral-error').classList.toggle('error-value', error === null || error >= (state.rest ? .05 : .03));
    $('spectral-weight').textContent = lr.total_weight.toFixed(4);
    const conv = point.convergence.find(r => r.eta === state.eta);
    const passed = conv.reference_certified;
    $('reference-note').classList.toggle('reference-passed', passed);
    $('reference-note').textContent = `VED 参考${passed ? '通过本次已测试收敛判据' : '尚未通过本次收敛判据'}：Nₕ = 18 → 20 的全谱 L¹ 变化为 ${percentage(conv.spectral_l1)}。${state.eta >= .1 ? '当前使用补充展宽，不能替代原定高分辨率目标。' : ''}LR-VMC 的整体物理收敛尚未认证；谱差异未附统计误差条。`;
    const allLR = data.momenta.map(p => selectLR(p, state.samples, state.rcond));
    const momentum = data.momenta.map(p => p.k_pi);
    for (const [element, prop, label] of [['dispersion-plot', 'E', 'E(k) / t'], ['residue-plot', 'Z', 'Z(k)']]) {
      drawPlot($(element), [
        { x: momentum, y: data.momenta.map(p => p.reference[prop]), color: colors.ved, name: 'VED' },
        { x: momentum, y: allLR.map(p => p[prop]), color: colors.lr, dashed: true, name: 'LR-VMC' },
      ], { small: true, fromZero: prop === 'Z', xMin: -.04, xMax: 1.04,
        yLabel: label, title: `${state.case} 组 ${label}`,
        description: '仅有三个已计算动量点。圆点为 VED，方点为 LR-VMC。折线仅辅助阅读。' });
    }
    $('vmc-table').innerHTML = data.momenta.map(p => `<tr><th scope="row">${kLabel(p.k_pi)}</th><td>${p.vmc.energy.mean.toFixed(7)} ± ${p.vmc.energy.chain_standard_error.toFixed(7)}</td><td>${p.vmc.Z.mean.toFixed(7)} ± ${p.vmc.Z.chain_standard_error.toFixed(7)}</td><td>${p.vmc.nph.mean.toFixed(5)}</td></tr>`).join('');
    $('raw-data-link').href = `data/${state.case}.json`;
    $('raw-data-link').setAttribute('aria-disabled', 'false');
    downloadable = { ...state, omega: data.omega, reference: refValues, candidate: lrValues,
      referenceRest: refRest, candidateRest: lrRest };
    $('download-spectrum').disabled = false;
    $('load-status').textContent = '';
  } catch (error) {
    if (version !== generation) return;
    $('case-label').textContent = `CASE ${state.case} / 加载失败`;
    $('load-status').replaceChildren(document.createTextNode(`图表未能更新：${error.message}。`));
    const retry = document.createElement('button');
    retry.className = 'text-button';
    retry.textContent = '重试';
    retry.addEventListener('click', render);
    $('load-status').append(retry);
    // Clear stale output so it cannot be mistaken for the newly selected case.
    for (const id of ['spectrum-plot', 'dispersion-plot', 'residue-plot', 'vmc-table', 'reference-note']) $(id).replaceChildren();
    for (const id of ['energy-ref', 'energy-lr', 'z-ref', 'z-lr', 'spectral-error', 'spectral-weight']) $(id).textContent = '—';
  } finally {
    if (version === generation) $('results-panel').setAttribute('aria-busy', 'false');
  }
}

document.querySelectorAll('[data-case]').forEach(button => button.addEventListener('click', () => {
  state.case = button.dataset.case;
  document.querySelectorAll('[data-case]').forEach(b => b.setAttribute('aria-pressed', String(b === button)));
  render();
}));
document.querySelectorAll('[data-k]').forEach(button => button.addEventListener('click', () => {
  state.k = Number(button.dataset.k);
  document.querySelectorAll('[data-k]').forEach(b => b.setAttribute('aria-pressed', String(b === button)));
  render();
}));
for (const id of ['eta', 'samples', 'rcond']) $(id).addEventListener('change', event => {
  state[id] = Number(event.target.value);
  render();
});
$('rest').addEventListener('change', event => { state.rest = event.target.checked; render(); });
$('download-spectrum').addEventListener('click', () => {
  if (!downloadable) return;
  const d = downloadable;
  const lines = ['omega_over_t,A_VED_times_t,A_LR_times_t,A_rest_VED_times_t,A_rest_LR_times_t'];
  for (let i = 0; i < d.omega.length; i++) {
    lines.push([d.omega[i], d.reference[i], d.candidate[i], d.referenceRest[i], d.candidateRest[i]].join(','));
  }
  const url = URL.createObjectURL(new Blob([lines.join('\n') + '\n'], { type: 'text/csv;charset=utf-8' }));
  const link = document.createElement('a');
  link.href = url;
  link.download = `Holstein_${d.case}_k${d.k}_eta${d.eta}_N${d.samples}_rc${d.rcond}.csv`;
  document.body.append(link);
  link.click();
  link.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
});

render();
