/** Cross-check the browser's numerical transformations against saved Python results. */
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { broaden, subtractLowest, relativeL1 } from '../site/assets/spectra.mjs';

let count = 0;
let maxDifference = 0;
for (const name of ['A', 'B', 'C', 'D']) {
  const data = JSON.parse(await readFile(new URL(`../site/data/${name}.json`, import.meta.url)));
  for (const point of data.momenta) {
    for (const eta of [.025, .05]) {
      const ref = broaden(point.reference.poles, data.omega, eta);
      const rest = subtractLowest(ref, data.omega, point.reference.E, point.reference.Z, eta);
      for (const lr of point.lr) {
        const candidate = broaden(lr.poles, data.omega, eta);
        const candidateRest = subtractLowest(candidate, data.omega, lr.E, lr.Z, eta);
        const expected = lr.metrics.find(r => r.eta === eta);
        assert.ok(expected);
        for (const [actual, recorded] of [
          [relativeL1(ref, candidate, data.omega), expected.spectral_l1],
          [relativeL1(rest, candidateRest, data.omega), expected.rest_spectral_l1],
        ]) {
          const delta = Math.abs(actual - recorded);
          assert.ok(delta < 1e-9, `${name}, k=${point.k_pi}, eta=${eta}, N=${lr.samples}, cutoff=${lr.rcond}: ${delta}`);
          maxDifference = Math.max(maxDifference, delta);
        }
        count++;
      }
    }
  }
}
console.log(`${count} whole/rest-spectrum comparisons match the Python records; max difference = ${maxDifference.toExponential(3)}`);
