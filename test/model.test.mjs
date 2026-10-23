// Check the JavaScript model against values computed by the Julia model (see make_reference.jl).
// Run with `npm test`, or `node --test`.

import { test } from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import {
  cd69Fraction,
  discriminationPower,
  localP15Slope,
  DEFAULT_PARAMS,
  PEPTIDES,
} from "../js/model.js";

const ref = JSON.parse(readFileSync(new URL("./reference.json", import.meta.url), "utf8"));

// Agreement to a relative error of `rtol`; null in the reference stands for NaN.
function close(actual, expected, rtol, what) {
  if (expected === null) return assert.ok(Number.isNaN(actual), `${what}: expected NaN, got ${actual}`);
  const err = Math.abs(actual - expected) / Math.max(Math.abs(expected), 1e-300);
  assert.ok(err <= rtol, `${what}: ${actual} vs ${expected} (relative error ${err})`);
}

test("default values match the Julia set", () => {
  assert.deepEqual({ ...DEFAULT_PARAMS }, ref.sets[0].params);
});

test("peptide KDs match the Julia set", () => {
  assert.deepEqual(PEPTIDES.map((q) => q.KD), ref.KDs.slice(1, -1));
});

ref.sets.forEach((set, s) => {
  test(`parameter set ${s}: fraction CD69+ over a grid of pMHC and KD`, () => {
    let i = 0;
    for (const KD of ref.KDs)
      for (const x of ref.xs) close(cd69Fraction(x, KD, set.params), set.fraction[i++], 1e-12, `x=${x}, KD=${KD}`);
  });

  test(`parameter set ${s}: P15, discrimination power and local slopes`, () => {
    const dp = discriminationPower(set.params);
    // Bisection stops at a width near 1e-17 in log10, so these agree to about 1e-12 at worst.
    dp.log10P15.forEach((v, i) => close(v, set.log10_p15[i], 1e-10, `log10 P15 of peptide ${i}`));
    close(dp.slope, set.slope, 1e-9, "slope");
    close(dp.intercept, set.intercept, 1e-9, "intercept");
    // A central difference of bisection results amplifies their rounding, so allow more.
    PEPTIDES.forEach((q, i) => close(localP15Slope(q.KD, set.params), set.local_slopes[i], 1e-6, `local slope at ${q.name}`));
  });
});
