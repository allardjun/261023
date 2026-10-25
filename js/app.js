// Page logic: build the sliders, recompute the model on every change, and redraw the 3d plot with Plotly.
// Plotly is loaded as a classic script in index.html, so it is `window.Plotly` here.

import {
  DEFAULT_PARAMS,
  PEPTIDES,
  cd69Fraction,
  discriminationPower,
  localP15Slope,
  log10Potency,
} from "./model.js";

const Plotly = window.Plotly;

// One slider per parameter.
// A "log" slider moves in log10 steps; the hill exponents and the two fractions move in plain steps.
// The bounds cover the default values and a broad range around them.
const GROUPS = [
  { id: "binding", title: "Binding" },
  { id: "proofreading", title: "Proofreading" },
  { id: "switch", title: "Switch and output" },
];

const SLIDERS = [
  { name: "fA", label: "f<sub>A</sub>", unit: "", scale: "lin", bounds: [0, 1], group: "binding",
    text: "Fraction of receptors that are free to bind." },
  { name: "alpha", label: "α", unit: "µm⁻² per µM", scale: "log", bounds: [1e-1, 1e3], group: "binding",
    text: "Converts the pMHC concentration used to load the cells into a surface density." },
  { name: "TT", label: "T<sub>T</sub>", unit: "µm⁻²", scale: "log", bounds: [1e-1, 1e3], group: "binding",
    text: "Total density of T cell receptors on the cell surface." },
  { name: "g1", label: "g̃<sub>1</sub>", unit: "µM⁻¹ s⁻¹", scale: "log", bounds: [1e-3, 1e2], group: "binding",
    text: "Sets how fast a complex lets go: the off-rate is g̃₁ times KD." },
  { name: "kon2d", label: "k<sub>on</sub><sup>2d</sup>", unit: "µm² s⁻¹", scale: "log", bounds: [1e-3, 1e1], group: "binding",
    text: "Rate at which a receptor and a peptide meet and bind in the contact between cells." },
  { name: "kP", label: "k<sub>P</sub>", unit: "s⁻¹", scale: "log", bounds: [1e-2, 1e2], group: "proofreading",
    text: "Rate of each proofreading step." },
  { name: "nKP", label: "n<sub>KP</sub>", unit: "", scale: "lin", bounds: [0.1, 6], group: "proofreading",
    text: "Number of proofreading steps a complex must pass before it signals." },
  { name: "lambdaX", label: "λ<sub>X</sub>", unit: "µm⁻²", scale: "log", bounds: [1e-6, 1e1], group: "switch",
    text: "Density of signalling complexes at which the switch is half on." },
  { name: "nC", label: "n<sub>C</sub>", unit: "", scale: "lin", bounds: [0.1, 5], group: "switch",
    text: "Steepness of the switch." },
  { name: "X_O1", label: "X<sub>O1</sub>", unit: "", scale: "log", bounds: [1e-5, 1e1], group: "switch",
    text: "Switch level at which the output is half its top." },
  { name: "O1max", label: "O<sub>1</sub><sup>max</sup>", unit: "", scale: "lin", bounds: [0, 1], group: "switch",
    text: "Largest possible fraction of CD69+ cells." },
];

const SLIDER_STEPS = 500;

// The pMHC axis runs past the tested range (1e-8 to 0.1 µM), so that every curve reaches its top.
const LOG10X_RANGE = [-8, 2];
const LOG10KD_RANGE = [0, 3];
const LOG10X = linspace(...LOG10X_RANGE, 400);
const XS = LOG10X.map((v) => 10 ** v);
const KD_FINE = linspace(...LOG10KD_RANGE, 200).map((v) => 10 ** v);

// Viridis, reversed, so the weakest peptide is darkest; the dark theme uses a lighter stretch of it.
const PEPTIDE_COLORS = {
  light: ["#9ad93d", "#2bb17e", "#277e8e", "#3f4788", "#440154"],
  dark: ["#fde725", "#7dd24f", "#24aa83", "#287c8e", "#5b5fb0"],
};

// Each curve's projection on the floor is drawn in its own hue, darkened by this fraction, so it reads as a shadow rather than as more of the curve.
const SHADOW = { light: 0.3, dark: 0.45 };

let params = { ...DEFAULT_PARAMS, ...paramsFromHash() };
const plotEl = document.getElementById("plot");

function linspace(a, b, n) {
  return Array.from({ length: n }, (_, i) => a + ((b - a) * i) / (n - 1));
}

const toSlider = (s, v) => (s.scale === "log" ? Math.log10(v) : v);
const fromSlider = (s, v) => (s.scale === "log" ? 10 ** v : v);

// Three significant figures, in a form that also reads back as a number.
function fmt(v) {
  if (!Number.isFinite(v)) return "n/a";
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 1e-3 && a < 1e4) return String(Number(v.toPrecision(3)));
  return v.toExponential(2).replace("e+", "e");
}

function inBounds(s, v) {
  return Number.isFinite(v) && v >= s.bounds[0] && v <= s.bounds[1] && !(s.scale === "log" && v <= 0);
}

// Theme

function theme() {
  const css = getComputedStyle(document.documentElement);
  const get = (name) => css.getPropertyValue(name).trim();
  const dark = window.matchMedia("(prefers-color-scheme: dark)").matches
    ? document.documentElement.dataset.theme !== "light"
    : document.documentElement.dataset.theme === "dark";
  return {
    text: get("--text"),
    muted: get("--muted"),
    grid: get("--grid"),
    border: get("--border"),
    fit: get("--fit"),
    font: get("--font"),
    peptides: dark ? PEPTIDE_COLORS.dark : PEPTIDE_COLORS.light,
    shadow: dark ? SHADOW.dark : SHADOW.light,
  };
}

// The same hue, mixed toward black by the fraction k.
function darken(hex, k) {
  const n = parseInt(hex.slice(1), 16);
  const c = [(n >> 16) & 255, (n >> 8) & 255, n & 255].map((v) => Math.round(v * (1 - k)));
  return `rgb(${c.join(", ")})`;
}

// Controls

function buildControls() {
  const form = document.getElementById("controls");
  for (const g of GROUPS) {
    const fs = document.createElement("fieldset");
    fs.innerHTML = `<legend>${g.title}</legend>`;
    for (const s of SLIDERS.filter((s) => s.group === g.id)) {
      const [lo, hi] = s.bounds.map((b) => (s.scale === "log" ? Math.log10(b) : b));
      const row = document.createElement("div");
      row.className = "control";
      row.title = s.text;
      row.innerHTML = `
        <label for="range-${s.name}">${s.label}${s.unit ? `<span class="unit">${s.unit}</span>` : ""}</label>
        <input type="range" id="range-${s.name}" min="${lo}" max="${hi}" step="${(hi - lo) / SLIDER_STEPS}"
               aria-describedby="text-${s.name}">
        <input type="text" id="text-${s.name}" inputmode="decimal" spellcheck="false"
               aria-label="${s.name} value">`;
      fs.append(row);
      const range = row.querySelector('input[type="range"]');
      const text = row.querySelector('input[type="text"]');
      s.range = range;
      s.textInput = text;
      range.addEventListener("input", () => {
        params[s.name] = fromSlider(s, Number(range.value));
        text.value = fmt(params[s.name]);
        text.removeAttribute("aria-invalid");
        setFill(s);
        schedule();
      });
      const commit = () => {
        const v = Number(text.value);
        if (!inBounds(s, v)) {
          text.setAttribute("aria-invalid", "true");
          return;
        }
        text.removeAttribute("aria-invalid");
        params[s.name] = v;
        syncControl(s);
        schedule();
      };
      text.addEventListener("change", commit);
      text.addEventListener("keydown", (e) => {
        if (e.key === "Enter") {
          e.preventDefault();
          commit();
        }
      });
    }
    form.append(fs);
  }
  form.addEventListener("submit", (e) => e.preventDefault());
}

function setFill(s) {
  const r = s.range;
  const f = ((Number(r.value) - Number(r.min)) / (Number(r.max) - Number(r.min))) * 100;
  r.style.setProperty("--fill", `${f}%`);
}

function syncControl(s) {
  s.range.value = toSlider(s, params[s.name]);
  s.textInput.value = fmt(params[s.name]);
  setFill(s);
}

function buildGlossary() {
  const dl = document.getElementById("glossary");
  dl.innerHTML = SLIDERS.map((s) => `<dt>${s.label}</dt><dd>${s.text}${s.unit ? ` (${s.unit})` : ""}</dd>`).join("");
}

// Sharing: the values live in the address, so a copied link reopens them.

function paramsFromHash() {
  const out = {};
  const q = new URLSearchParams(location.hash.slice(1));
  for (const s of SLIDERS) {
    if (!q.has(s.name)) continue;
    const v = Number(q.get(s.name));
    if (inBounds(s, v)) out[s.name] = v;
  }
  return out;
}

// Only values that differ from the defaults go in, so links stay short.
function writeHash() {
  const changed = SLIDERS.filter((s) => params[s.name] !== DEFAULT_PARAMS[s.name]);
  const hash = changed.length ? "#" + new URLSearchParams(changed.map((s) => [s.name, String(params[s.name])])) : "";
  history.replaceState(null, "", location.pathname + location.search + hash);
}

// Plot

function traces(p, t) {
  const zeros = XS.map(() => 0);
  const dp = discriminationPower(p);
  const out = [];

  PEPTIDES.forEach((q, i) => {
    const color = t.peptides[i];
    const ys = XS.map(() => q.KD);
    const zs = XS.map((x) => cd69Fraction(x, q.KD, p));
    out.push({
      type: "scatter3d", mode: "lines", x: XS, y: ys, z: zeros,
      line: { color: darken(color, t.shadow), width: 2 }, hoverinfo: "skip",
    });
    out.push({
      type: "scatter3d", mode: "lines", x: XS, y: ys, z: zs, name: q.name,
      line: { color, width: 6 },
      hovertemplate: `<b>${q.name}</b>, KD ${q.KD} µM<br>pMHC %{x:.3~e} µM<br>CD69+ %{z:.3f}<extra></extra>`,
    });
    out.push({
      type: "scatter3d", mode: "text", x: [XS[0]], y: [q.KD], z: [zs[0]], text: [q.name],
      textposition: "middle left", textfont: { color, size: 13, family: t.font }, hoverinfo: "skip",
    });
  });

  // Potency P15 on the floor, with gaps where the response never reaches 15%.
  out.push({
    type: "scatter3d", mode: "lines",
    x: KD_FINE.map((KD) => {
      const v = log10Potency(KD, p);
      return Number.isFinite(v) ? 10 ** v : null;
    }),
    y: KD_FINE, z: KD_FINE.map(() => 0),
    line: { color: t.muted, width: 2, dash: "dash" },
    hovertemplate: "P<sub>15</sub> %{x:.3~e} µM<br>KD %{y:.3~f} µM<extra></extra>",
  });

  // Each peptide's P15, and the straight-line fit through them.
  out.push({
    type: "scatter3d", mode: "markers",
    x: dp.log10P15.map((v) => (Number.isFinite(v) ? 10 ** v : null)),
    y: PEPTIDES.map((q) => q.KD), z: PEPTIDES.map(() => 0),
    marker: { color: t.peptides, size: 4, line: { width: 0 } },
    text: PEPTIDES.map((q) => q.name),
    hovertemplate: "<b>%{text}</b> P<sub>15</sub> %{x:.3~e} µM<extra></extra>",
  });
  if (Number.isFinite(dp.slope)) {
    const ends = [PEPTIDES[0].KD, PEPTIDES[PEPTIDES.length - 1].KD];
    const fx = ends.map((KD) => 10 ** (dp.intercept + dp.slope * Math.log10(KD)));
    out.push({
      type: "scatter3d", mode: "lines", x: fx, y: ends, z: [0, 0],
      line: { color: t.fit, width: 4 }, hoverinfo: "skip",
    });
    out.push({
      type: "scatter3d", mode: "text", x: [fx[1]], y: [ends[1]], z: [0],
      text: [`slope ${dp.slope.toFixed(2)}`], textposition: "middle right",
      textfont: { color: t.fit, size: 13, family: t.font }, hoverinfo: "skip",
    });
  }
  return { data: out, dp };
}

// Tick labels written out as powers of ten.
// Plotly's own labels for log axes in 3d fail to draw the label "1", so every label is given explicitly.
const SUPERSCRIPT = { "-": "⁻", 0: "⁰", 1: "¹", 2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹" };
function powerTicks([lo, hi], step) {
  const exps = [];
  for (let e = lo; e <= hi; e += step) exps.push(e);
  return {
    tickvals: exps.map((e) => 10 ** e),
    ticktext: exps.map((e) => "10" + [...String(e)].map((c) => SUPERSCRIPT[c]).join("")),
  };
}

function layout(t) {
  const axis = (extra) => ({
    color: t.text,
    gridcolor: t.grid,
    linecolor: t.border,
    zerolinecolor: t.border,
    showbackground: false,
    showspikes: false,
    tickfont: { size: 11, color: t.muted },
    title: { font: { size: 13, color: t.text } },
    ...extra,
  });
  return {
    paper_bgcolor: "rgba(0,0,0,0)",
    font: { family: t.font, color: t.text },
    margin: { l: 0, r: 0, t: 0, b: 0 },
    showlegend: false,
    uirevision: "keep",
    hoverlabel: { font: { family: t.font } },
    scene: {
      uirevision: "keep",
      xaxis: axis({ type: "log", range: LOG10X_RANGE, ...powerTicks(LOG10X_RANGE, 2),
                    title: { text: "pMHC concentration (µM)", font: { size: 13, color: t.text } } }),
      yaxis: axis({ type: "log", range: LOG10KD_RANGE, ...powerTicks(LOG10KD_RANGE, 1),
                    title: { text: "Peptide KD (µM)", font: { size: 13, color: t.text } } }),
      // The zero tick is left off, since it collides with the corner ticks of the other axes; the floor marks zero.
      zaxis: axis({ range: [0, 1], tickvals: [0.25, 0.5, 0.75, 1], tickformat: ".2f",
                    title: { text: "Fraction CD69+", font: { size: 13, color: t.text } } }),
      aspectmode: "manual",
      aspectratio: { x: 1.6, y: 1.1, z: 0.8 },
      camera: defaultCamera(),
    },
  };
}

// Front right, looking slightly down, so pMHC runs along the front edge and KD recedes to the back right.
// On a narrow screen the camera starts further back, so the whole box fits.
function defaultCamera() {
  const narrow = plotEl.clientWidth < 600;
  const k = narrow ? 1.3 : 1;
  return { eye: { x: 1.1 * k, y: -2.3 * k, z: 0.85 * k }, center: { x: narrow ? 0.04 : -0.12, y: 0, z: -0.16 } };
}

const PLOT_CONFIG = {
  responsive: true,
  displaylogo: false,
  modeBarButtonsToRemove: ["resetCameraLastSave3d", "hoverClosest3d"],
  toImageButtonOptions: { format: "svg", filename: "cd69_response" },
};

// Readouts

function renderReadouts(p, dp) {
  document.getElementById("power").textContent = Number.isFinite(dp.slope) ? dp.slope.toFixed(2) : "n/a";
  const first = PEPTIDES[0];
  const last = PEPTIDES[PEPTIDES.length - 1];
  const s1 = localP15Slope(first.KD, p);
  const s2 = localP15Slope(last.KD, p);
  document.getElementById("local-slopes").textContent =
    `Slope of the straight-line fit of log P15 against log KD. ` +
    `The dashed line's own slope is ${Number.isFinite(s1) ? s1.toFixed(2) : "n/a"} at ${first.name} ` +
    `and ${Number.isFinite(s2) ? s2.toFixed(2) : "n/a"} at ${last.name}.`;
}

// Redraw at most once per frame, however fast the sliders move.

let pending = false;
function schedule() {
  if (pending) return;
  pending = true;
  requestAnimationFrame(() => {
    pending = false;
    render();
  });
}

function render() {
  const t = theme();
  const { data, dp } = traces(params, t);
  Plotly.react(plotEl, data, layout(t), PLOT_CONFIG);
  renderReadouts(params, dp);
  writeHash();
}

function resetToDefaults() {
  params = { ...DEFAULT_PARAMS };
  SLIDERS.forEach(syncControl);
  SLIDERS.forEach((s) => s.textInput.removeAttribute("aria-invalid"));
  schedule();
}

async function copyLink(button) {
  writeHash();
  const label = button.textContent;
  try {
    await navigator.clipboard.writeText(location.href);
    button.textContent = "Link copied";
  } catch {
    button.textContent = "Could not copy; use the address bar";
  }
  setTimeout(() => (button.textContent = label), 1800);
}

function init() {
  buildControls();
  buildGlossary();
  SLIDERS.forEach(syncControl);
  document.getElementById("reset").addEventListener("click", resetToDefaults);
  document.getElementById("reset-view").addEventListener("click", () =>
    Plotly.relayout(plotEl, { "scene.camera": defaultCamera() }));
  const copy = document.getElementById("copy-link");
  copy.addEventListener("click", () => copyLink(copy));
  window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", schedule);
  window.addEventListener("hashchange", () => {
    params = { ...DEFAULT_PARAMS, ...paramsFromHash() };
    SLIDERS.forEach(syncControl);
    schedule();
  });
  render();
}

init();
