/** Lorentz broadening of every saved pole, without shifting or normalizing. */
export function broaden(poles, omega, eta) {
  if (!(eta > 0)) throw new Error('展宽必须大于零');
  const values = new Float64Array(omega.length);
  for (let i = 0; i < poles.energies.length; i++) {
    const e = poles.energies[i];
    const w = poles.weights[i] * eta / Math.PI;
    for (let j = 0; j < omega.length; j++) {
      values[j] += w / ((omega[j] - e) ** 2 + eta ** 2);
    }
  }
  return values;
}

export function subtractLowest(values, omega, energy, residue, eta) {
  return Float64Array.from(values, (value, j) =>
    value - residue * eta / Math.PI / ((omega[j] - energy) ** 2 + eta ** 2));
}

export function integrate(values, omega) {
  let result = 0;
  for (let i = 1; i < omega.length; i++) {
    result += (values[i - 1] + values[i]) * (omega[i] - omega[i - 1]) / 2;
  }
  return result;
}

export function relativeL1(reference, candidate, omega) {
  const mass = integrate(reference, omega);
  if (mass <= 1e-10) return null;
  return integrate(Float64Array.from(reference, (v, i) => Math.abs(v - candidate[i])), omega) / mass;
}

export function selectLR(point, samples, rcond) {
  const record = point.lr.find(r => r.samples === samples && r.rcond === rcond);
  if (!record) throw new Error('这组采样与截断设置没有已完成的数据');
  return record;
}
