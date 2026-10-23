# CD69 response explorer

An interactive page for a kinetic proofreading model of T cell activation.
Sliders set the model's parameters, and a 3d plot shows the predicted fraction of CD69+ cells against pMHC concentration for five peptides of different binding strength.
The page also reports the discrimination power, as defined by Pettmann et al. ([eLife 2021](https://doi.org/10.7554/eLife.67092)).

Live page: <https://allardjun.github.io/261023/>

Everything runs in the visitor's browser.
There is no server and no build step: the page is static HTML, CSS and JavaScript, with [Plotly](https://plotly.com/javascript/) loaded from a CDN.

## Default values

The default parameter values come from one fit of the model to CD69 data.
The fitting method is not yet published, so the page presents them as a starting point, not as a result.

## Layout

- `index.html`: the page, including the model's equations.
- `css/style.css`: styles, with light and dark themes.
- `js/model.js`: the model, P15 and the discrimination power, with no dependencies.
- `js/app.js`: sliders, readouts, shareable links and the plot.
- `test/`: checks that `js/model.js` agrees with the Julia version of the model.

## Running locally

The page loads its scripts as ES modules, which browsers refuse to load from `file://`, so serve the folder:

```sh
python3 -m http.server 8000
```

Then open <http://localhost:8000>.

## Tests

`npm test` checks the JavaScript model against `test/reference.json` for three parameter sets, to a relative error of 1e-12.
It needs Node 20 or later and installs nothing.

The reference values come from the Julia version of the model, which lives in a separate repository.
To regenerate them after changing the model, run `julia test/make_reference.jl /path/to/model.jl`.

## Sources

- Peptide affinities: Aleksic et al., Immunity 2010.
- Discrimination power: Pettmann et al., eLife 2021, [doi:10.7554/eLife.67092](https://doi.org/10.7554/eLife.67092).
- The model is from a manuscript in preparation.
