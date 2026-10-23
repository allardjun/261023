// Base KP-X-O model of T cell activation, in closed form.
// It follows `get_CT`, `get_CN`, `get_X` and `get_O1` in the Julia model, and test/model.test.mjs checks the two agree.
//
// Units: pMHC concentration x in µM, KD in µM, surface densities in µm⁻², rates in s⁻¹.
// In the code, `g1` is g̃1 = g1·kon3d, so that the off-rate is koff = g1·KD.

// Default values: one parameter set from fitting the model to CD69 data.
// The fitting procedure is not yet published, so these are a starting point, not a result.
export const DEFAULT_PARAMS = Object.freeze({
  fA: 0.0010760218351601953,
  alpha: 149.34436781578844,
  TT: 1.094630372285705,
  g1: 0.03204443607404042,
  kon2d: 0.08674689927211769,
  kP: 1.6476071717517033,
  nKP: 2.0248129398409676,
  lambdaX: 0.0177085161359101,
  nC: 0.8070635083920701,
  X_O1: 0.001188267819200772,
  O1max: 0.8571992797935549,
});

// Peptide KDs (µM), from Aleksic et al. 2010.
export const PEPTIDES = Object.freeze([
  { name: "9V", KD: 7.97 },
  { name: "6V", KD: 22.8 },
  { name: "3Y", KD: 47.3 },
  { name: "4D", KD: 140.0 },
  { name: "4A", KD: 299.0 },
]);

// Bound TCR–pMHC density C_T (µm⁻²), the smaller root of the Lauffenburger quadratic.
// This form equals (b − sqrt(b² − 4·PT·TT))/2, but keeps its precision when PT ≪ b, at low pMHC.
export function getCT(x, KD, p) {
  const PT = p.alpha * x;
  const b = PT + p.TT + (p.g1 * KD) / p.kon2d;
  return (p.fA * 2 * PT * p.TT) / (b + Math.sqrt(b * b - 4 * PT * p.TT));
}

// Fully proofread complex density C_N (µm⁻²) after nKP proofreading steps.
export function getCN(CT, KD, p) {
  return Math.pow(1 / (1 + (KD * p.g1) / p.kP), p.nKP) * CT;
}

// Output of the cytoplasmic switch X, in [0, 1].
export function getX(CN, p) {
  const a = Math.pow(CN, p.nC);
  return a / (a + Math.pow(p.lambdaX, p.nC));
}

// Fraction of CD69+ cells O1, in [0, O1max].
export function getO1(X, p) {
  return (p.O1max * X) / (p.X_O1 + X);
}

// Fraction of CD69+ cells for pMHC concentration x (µM) and dissociation constant KD (µM).
export function cd69Fraction(x, KD, p) {
  return getO1(getX(getCN(getCT(x, KD, p), KD, p), p), p);
}

// The four combinations of parameters that the output depends on.
export function compositeParams(p) {
  return {
    xstar: p.TT / p.alpha, // µM
    ystar: p.lambdaX / (p.fA * p.TT), // dimensionless
    KD_Lauf: (p.TT * p.kon2d) / p.g1, // µM
    KD_McK: p.kP / p.g1, // µM
  };
}

// Discrimination power, as in Pettmann et al. 2021 (eLife 10:e67092).
// Potency is P15, the pMHC concentration that makes 15% of cells CD69+ (an absolute level, not 15% of each peptide's own top).
// The discrimination power is the slope of a straight-line fit of log10 P15 against log10 KD, over the peptide panel.

// log10 of the pMHC concentration (µM) at which the response reaches `level`, found by bisection; NaN if it never does.
export function log10Potency(KD, p, { level = 0.15, lo = -12, hi = 4 } = {}) {
  if (!(cd69Fraction(10 ** hi, KD, p) >= level)) return NaN;
  for (let i = 0; i < 60; i++) {
    const mid = (lo + hi) / 2;
    if (cd69Fraction(10 ** mid, KD, p) < level) lo = mid;
    else hi = mid;
  }
  return (lo + hi) / 2;
}

// Fit log10 P15 = intercept + slope·log10 KD over the peptides that reach 15%.
// Peptides that never reach 15% are left out; with fewer than two left, the slope is NaN.
export function discriminationPower(p, KDs = PEPTIDES.map((q) => q.KD)) {
  const log10P15 = KDs.map((KD) => log10Potency(KD, p));
  const pts = KDs.map((KD, i) => [Math.log10(KD), log10P15[i]]).filter(([, y]) => Number.isFinite(y));
  if (pts.length < 2) return { slope: NaN, intercept: NaN, log10P15 };
  const n = pts.length;
  const mx = pts.reduce((s, [x]) => s + x, 0) / n;
  const my = pts.reduce((s, [, y]) => s + y, 0) / n;
  const sxy = pts.reduce((s, [x, y]) => s + (x - mx) * (y - my), 0);
  const sxx = pts.reduce((s, [x]) => s + (x - mx) ** 2, 0);
  const slope = sxy / sxx;
  return { slope, intercept: my - slope * mx, log10P15 };
}

// Local slope d log10 P15 / d log10 KD at one KD, by a central difference.
export function localP15Slope(KD, p, h = 1e-3) {
  return (log10Potency(KD * 10 ** h, p) - log10Potency(KD / 10 ** h, p)) / (2 * h);
}
