# T cell activation explorer: notes for Claude

This repository is public, and so is this file.
Notes that must stay private (unpublished results, local paths, the upstream fit) go in `CLAUDE.local.md`, which git ignores.

## What this is

A static page, served by GitHub Pages at <https://allardjun.github.io/261023/>, that plots a closed-form model of T cell activation.
Eleven sliders set the parameters; a Plotly 3d plot shows the fraction of CD69+ cells against pMHC concentration for five peptides, placed by KD.
The page also reports the discrimination power.
Everything runs in the visitor's browser: no server, no build step, no bundler, no npm dependencies.

The model was ported from a Julia version in a separate, private repository.
The two projects are meant to diverge from here on.
The Julia version is only needed to regenerate `test/reference.json`.

## Wording rules for the page

- Call the fitted parameter set **"default values"** (or "the defaults"), never "fitted", "learned" or "best-fit".
  The fitting procedure is not peer-reviewed or published, so the page presents the values as a starting point, not a result.
  The README says so; the page itself no longer comments on where the defaults come from, at the owner's request.
- Write for a reader outside the field: short, concise sentences, without software or molecular-biology jargon.
  Explain a term the first time it appears (for example P15, KD, discrimination power).
- Cite sources by name and year in the text, with a DOI link where there is one.
  The page describes the model as "an extension of Lever et al. 2014" ([Nature Reviews Immunology](https://www.nature.com/articles/nri3728)); do not describe its unpublished results.

## Conventions

- Semantic line breaks everywhere, including HTML text, CSS and JS comments, and Markdown: one sentence per line, no hard wrapping.
- Atomic commits, each doing one thing, with a subject in the imperative mood.
- Colors live as CSS custom properties on `:root`, redefined for dark mode under `prefers-color-scheme` (guarded by `:root:not([data-theme="light"])`) and again under `:root[data-theme="dark"]`.
  `js/app.js` reads them at render time to theme the plot, so a new color goes in CSS first.
- Peptide colors are a reversed viridis, with a lighter stretch for the dark theme (`PEPTIDE_COLORS` in `js/app.js`).
- The layout must work at phone width (about 400 px) with a 16 px gutter and no sideways scrolling.
  Below 960 px the panel drops under the plot; below 600 px the plot is square.

## The model (`js/model.js`)

The pipeline, all in closed form:

1. Bound complexes `C_T`: the smaller root of the Lauffenburger quadratic, times `fA`.
   It is written as `2·fA·PT·TT / (b + sqrt(b² − 4·PT·TT))`, which equals the textbook `(b − sqrt(...))/2` but keeps its precision at low pMHC, where the textbook form subtracts two nearly equal numbers.
2. Proofreading: `C_N = (1 + g1·KD/kP)^(−nKP) · C_T`.
3. Switch: `X = C_N^nC / (C_N^nC + λX^nC)`.
4. Output: `O1 = O1max · X / (X_O1 + X)`.

Units: pMHC `x` and `KD` in µM, densities in µm⁻², rates in s⁻¹.
The code's `g1` is g̃₁ = g₁·kon³ᵈ, so the off-rate is `g1·KD`; the page labels it g̃₁.

Facts worth keeping in mind when changing the page:

- Seven parameters (`fA`, `alpha`, `TT`, `kon2d`, `g1`, `kP`, `lambdaX`) enter only through four combinations: x* = TT/α, y* = λX/(fA·TT), KD_Lauf = TT·kon2d/g̃₁ and KD_McK = kP/g̃₁.
  The panel lists these under "Combinations that matter".
  So many slider moves are equivalent, and that is expected, not a bug.
- `X` followed by `O1` is itself one Hill switch in `C_N`, with plateau `O1max/(1 + X_O1)` and threshold `λX·(X_O1/(1 + X_O1))^(1/nC)`.
  So `lambdaX`, `X_O1` and `O1max` trade off against one another.
- The response always rises as KD falls; there is no penalty for very tight binders in this model.

### Discrimination power

Following Pettmann et al., eLife 2021 ([doi:10.7554/eLife.67092](https://doi.org/10.7554/eLife.67092)):

- Potency P15 is the pMHC concentration that makes 15% of cells CD69+.
  It is an absolute level, not 15% of each peptide's own top.
- `log10Potency` finds it by bisection on log10 x in [−12, 4] (60 steps), and returns NaN when the response never reaches 15%.
- The discrimination power is the least-squares slope of log10 P15 against log10 KD over the five peptides, leaving out those that never reach 15%.
  With fewer than two left it is NaN, shown as "n/a".
- With the defaults it is 1.93.
- The dashed P15 line on the floor is curved: its local slope (`localP15Slope`, a central difference) is about 0.3 at KD = 1 µM, 1.0 at 9V, 2.8 at 4A, and 3.5 at 1000 µM.
- We use the panel fit, not the high-KD slope, on purpose.
  The fit is what Pettmann et al. measure, so the numbers compare directly.
  The local slope at the weak end blows up as a peptide's top response approaches 15%, which reflects the choice of threshold more than proofreading.
  The fitted slope does depend on the panel: theirs spans KD 7 to 1309 µM, ours 8 to 299 µM.

### Peptides

9V 7.97, 6V 22.8, 3Y 47.3, 4D 140, 4A 299 µM (Aleksic et al., Immunity 2010).
A sixth peptide, 6T (162 µM), exists in the data but was left out of the analysis this page follows.
The tested pMHC range is 1e-8 to 0.1 µM; the plot runs to 100 µM so every curve reaches its top, which is a prediction outside the data.

## Tests

`npm test` (which runs `node --test`) compares `js/model.js` with `test/reference.json`, for three parameter sets: the defaults, one that moves every parameter, and the nominal starting guess.
Tolerances, and why:

- Fraction CD69+ over a grid of x and KD: relative 1e-12.
- log10 P15: 1e-10, slope and intercept: 1e-9, since bisection stops at a width near 1e-17 in log10.
- Local slopes: 1e-6, since the central difference amplifies the bisection's rounding.

`null` in the JSON stands for NaN.
Do not change `node --test` to `node --test test/`: in Node 23 that form fails to find the files.
CI (`.github/workflows/test.yml`) runs the same on Node 22 for every push and pull request.

To regenerate the reference after an intended change to the model, run `julia test/make_reference.jl /path/to/model.jl` with the Julia model (see `CLAUDE.local.md` for the path).
If the JavaScript model diverges on purpose from the Julia one, the reference has to come from somewhere else, for example a frozen copy, so decide that explicitly rather than loosening tolerances.

## Plotly lessons

Plotly.js 2.35.2 (`plotly.js-dist-min`) is loaded from jsDelivr with a subresource-integrity hash.
To change the version, recompute the hash:

```sh
curl -sL https://cdn.jsdelivr.net/npm/plotly.js-dist-min@VERSION/plotly.min.js | openssl dgst -sha384 -binary | openssl base64 -A
```

Things that went wrong and how they were fixed:

- **Missing tick label on 3d log axes.**
  Plotly's automatic labels drew 10⁻¹ and 10¹ but not "1", and logged "error vectorizing text '1'".
  `powerTicks()` gives every label explicitly, as unicode superscripts.
- **Double-click does not reset a 3d camera**, unlike 2d plots.
  The "Reset view" button calls `Plotly.relayout(plotEl, {"scene.camera": defaultCamera()})`.
  The modebar's "reset camera to last save" button is removed, since it confused the two.
- **Keeping the user's rotation across redraws** needs `uirevision: "keep"` on both the layout and the scene; `Plotly.react` then leaves the camera alone.
- **Camera.**
  The default looks from the front right, slightly down, so pMHC runs along the front edge and KD recedes: eye (1.1, −2.3, 0.85), center (−0.12, 0, −0.16).
  Below 600 px the eye moves 1.3 times further out and the center shifts, so the whole box fits.
  Check both widths after any change to the axes or aspect ratio.
- **The zero tick on the z axis** collided with the corner ticks of the other two axes, so the z ticks are 0.25 to 1 and the floor marks zero.
- **Speed.**
  A full redraw (model plus `Plotly.react`) takes about 40 ms in Chrome on a laptop.
  Slider input is coalesced to at most one redraw per animation frame (`schedule()`).
  The model is under 4 ms of that (in Node: 0.4 ms for the curves, 3 ms for the P15 line, its 200 KDs times 60 bisection steps); Plotly's WebGL redraw is the rest.
  So speed work belongs in the plotting, not the model.

## Page behavior worth preserving

- Log sliders move in log10 space, with 500 steps; `nKP`, `nC`, `fA` and `O1max` move in plain steps.
- Each slider has a text box for an exact value; an out-of-range value gets `aria-invalid` and a red border, and is not applied.
- The address hash holds only the values that differ from the defaults, so shared links stay short; values from a hash are bounds-checked.
- The page re-renders when the system switches between light and dark.
- "Nominal" values were left off the page on purpose: they give almost no response over the plotted range, which would confuse a visitor.

## Checking changes in a browser

ES modules do not load from `file://`, so serve the folder with `python3 -m http.server 8000` and open <http://localhost:8000>.

- Check light and dark themes, a desktop width, and about 400 px wide (no sideways scroll).
- Browsers cache the old `app.js` aggressively during local testing; reload with the cache bypassed, or the change will seem not to work.
- Synthetic mouse drags from browser automation do not always rotate a WebGL canvas.
  Test rotation by hand, or check the camera through `plotEl.layout.scene.camera` after a relayout.
- Not yet tested: Safari on a real iPad, and Firefox.

## Deploying

GitHub Pages serves the `main` branch from the repository root (legacy build, no workflow), and `.nojekyll` makes it serve the files as they are.
A push to `main` goes live in about a minute.
Pages sends `Cache-Control: max-age=600`, so a visitor may see the old script for up to 10 minutes.

## Open questions

- No license yet, so the code is public but all rights reserved.
  Ask the owner which license before adding one.
- The default `fA` (0.00108) sits at the far left of its plain 0–1 slider, where it is hard to move finely.
  A log slider, for example 1e-4 to 1, would suit it better; this needs the owner's agreement, since plain steps for [0, 1] fractions was the original request.
