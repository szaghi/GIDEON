# Verification

Every capability of GIDEON is backed by a verification case: a Fortran program in `src/tests/` that runs the
compiled kernels and asserts the result against something it shares no code with: the monomials a row must be exact
on, a smooth function and its exact derivatives, a mirror or symmetry identity, ADAM's own code compiled verbatim, the
host for the device. The generator adds a case of its own (C1): it checks the rows before any Fortran exists. The
programs are the test suite (`bash scripts/run_tests.sh`), and with `GIDEON_REPORT_DIR` set they also write the data
from which these pages are drawn (`scripts/run_report.sh`). Every table and figure of this section is generated from
that data: rerun the script after any change of the generator or the kernels.

Each page has the same sections: case and scope, models and numerics with links to the code, the program and how to
run it, results, findings. The tag of a card says what kind of check it is:

- **contract**: an identity that must hold exactly or within a stated tolerance: rows by their meaning, exactness,
  mirrors and symmetries, dispatch, equivalence with ADAM, host and device;
- **accuracy**: an order of accuracy, a convergence rate, a spectral property.

<div class="gideon-casts">
<figure class="gideon-cast" style="--ar: 2.25">
<video class="light-only" autoplay loop muted playsinline poster="/verification/casts/advection-light.png" aria-label="Linear advection of a wave packet, centred and upwind first derivatives"><source src="/verification/casts/advection-light.webm" type="video/webm"><source src="/verification/casts/advection-light.mp4" type="video/mp4"></video>
<video class="dark-only" autoplay loop muted playsinline poster="/verification/casts/advection-dark.png" aria-label="Linear advection of a wave packet, centred and upwind first derivatives"><source src="/verification/casts/advection-dark.webm" type="video/webm"><source src="/verification/casts/advection-dark.mp4" type="video/mp4"></video>
<figcaption>A wave packet advected over two periods, 128 cells, by GIDEON's compiled first derivatives: the centred kernels disperse (c2 lags behind the exact packet, c4 nearly keeps up), the upwind one dissipates. <a href="./spectral">Case A3</a>.</figcaption>
</figure>
</div>

## What a user of the kernels should know

The findings that change how the kernels are used, each documented on its page:

- **The orders of every kernel are observed in quadruple precision only.** In double precision the round-off floor
  hides the asymptotic range of the wide stencils; the R16P build resolves all of them ([A1](./orders)). The
  generator measures the rows themselves at 100 digits ([C1](./rows)).
- **Only the one- and two-cell upwind first derivatives are stable advection operators on their own.** From three
  cells on, every fully one-sided first derivative amplifies a band of wavenumbers of a wave moving towards it
  ([A3](./spectral)). The wider rows are exact and kept as building blocks.
- **A centred row does not dissipate.** The pair form makes the imaginary part of the modified wavenumber zero to
  round-off, in the compiled kernels as in the rows ([A3](./spectral)); a centred scheme needs its dissipation from
  elsewhere.
- **The FV centred cell-average derivatives have the rows of the FD node derivatives** of the same window, exactly:
  the averaging factor of a Fourier mode cancels. The kernels differ in arithmetic (a flux difference, a pair form),
  not in the operator ([C1](./rows)).
- **ADAM's left-side upwind face value and first derivative are wrong for two or more cells**, and ADAM's upwind names
  flip between the first derivative and the others; GIDEON fixes both and is bitwise equal to ADAM wherever ADAM is
  right ([C5](./adam)).
- **A run-time order costs a `select case` per call**, and an order without a kernel returns a NaN, never a number
  ([C4](./dispatch)).

## Contracts

What must hold exactly, or bit for bit: the rows, their compiled values, their identities, the dispatch, ADAM's
values, the device.

<div class="gideon-gallery gideon-catalog">

[![The windows](/guide/windows-upwind.png)<strong>C1 Rows by their meaning</strong><span>every row derived in exact rationals, its order at 100 digits, ADAM's tables cross-checked, the mirrors proven</span><em>contract</em>](./rows)

[<strong>C2 Exactness and leading error</strong><span>every compiled kernel exact on the monomials up to its degree, and off by the derived error on the next</span><em>contract</em>](./exactness)

[<strong>C3 Mirrors and symmetries</strong><span>left side = mirrored right side, odd data to zero, power-of-two scalings: bitwise</span><em>contract</em>](./properties)

[<strong>C4 Run-time dispatch</strong><span>every dispatcher calls its kernel, bitwise; NaN for an order without one</span><em>contract</em>](./dispatch)

[![Equivalence with ADAM](/verification/adam.png)<strong>C5 Equivalence with ADAM</strong><span>ADAM's FDV operators compiled verbatim: bitwise where ADAM is right, its defects A–F reproduced</span><em>contract</em>](./adam)

[<strong>C6 Host and device</strong><span>every kernel in an OpenACC loop on the GPU against the host</span><em>contract</em>](./device)

</div>

## Accuracy

The design order of every kernel and every 3-D operator, and how the operators treat each Fourier mode.

<div class="gideon-gallery gideon-catalog">

[![Orders of accuracy](/verification/orders.png)<strong>A1 Orders of accuracy</strong><span>e(h) on a generic smooth function, R4P, R8P and R16P: where round-off hides the order, and where it does not</span><em>accuracy</em>](./orders)

[<strong>A2 3-D operators</strong><span>gradient, divergence, curl, Laplacian, FD and FV, every order, exact on polynomials with unequal steps</span><em>accuracy</em>](./operators-3d)

[![Spectral properties](/verification/spectral-upwind.png)<strong>A3 Spectral properties</strong><span>dispersion and dissipation of the first and second derivatives, from the rows and from the compiled kernels; a wave packet advected</span><em>accuracy</em>](./spectral)

</div>

## Regression

The values themselves, pinned by golden files, the suite's power to catch planted defects, and the CI that runs it
all: [Regression](/regression/).

## Running it

```bash
fobis build --mode tests-gnu && bash scripts/run_tests.sh   # the suite, double precision
fobis build --mode tests-gnu-r16 && bash scripts/run_tests.sh --exe-dir exe-r16   # the orders, quadruple precision
scripts/run_report.sh [--gpu]                                 # every build, the data, tables, figures, casts
python3 scripts/gen_gideon.py --check                         # C1: the rows, before any Fortran
```

The programs and the page of each:

| Program | Page | Builds |
|---|---|---|
| `scripts/gen_gideon.py --check` | [C1 Rows by their meaning](./rows) | Python 3, exact rationals |
| `gideon_test_tables` | [C2 Exactness and leading error](./exactness) | every |
| `gideon_test_properties` | [C3 Mirrors and symmetries](./properties) | every |
| `gideon_test_dispatch` | [C4 Run-time dispatch](./dispatch) | every |
| `gideon_test_adam` | [C5 Equivalence with ADAM](./adam) | R8P only |
| `gideon_test_device` | [C6 Host and device](./device) | every (on the device with OpenACC) |
| `gideon_test_order` | [A1 Orders of accuracy](./orders) | every |
| `gideon_test_operators_3d` | [A2 3-D operators](./operators-3d) | every |
| `gideon_test_spectral`, `gideon_vv_advection` | [A3 Spectral properties](./spectral) | every |
| `gideon_test_golden` | [Golden values](/regression/golden) | every (one golden file per kind) |
| `gideon_test_kinds` | the kinds reach a client through `use gideon` | every |

## From the data to the pages

No number of this section is typed. `scripts/run_report.sh` builds the test programs (R8P, R16P, R4P, and with
`--gpu` the OpenACC build), runs each with `GIDEON_REPORT_DIR` set so that it writes its data (CSV files, the frames of
the advection), keeps every output, and then draws the pages from them with the scripts of `docs/scripts`:

| Script | Reads | Writes |
|---|---|---|
| `make_tables.py` | the generator's JSON, the test outputs and CSV files, the golden files, the recorded mutation run (`notes/mutation-run-*.txt`), the benchmark data (`notes/bench-data`) | `docs/verification/_tables`, `docs/regression/_tables` |
| `make_figures.py` | the same | every figure of `docs/public/guide`, `docs/public/verification`, `docs/public/regression`, light and dark |
| `make_animations.py` | the frames of `gideon_vv_advection` | `docs/public/verification/casts`: WebM, MP4 and a poster per theme (needs ffmpeg) |
| `make_cookbook_tables.py` | `src/lib/*.F90` only (`--check`: exit 1 when a table is stale) | the tables of the cookbook and of the [API overview](/api/) (docs/ford.md) |

The figures share two modules: `figure_theme.py`, the light (washi) and dark (smoked bronze) styles, SETH's and
ADAM's copy with GIDEON's dark tokens ([its identity](https://github.com/szaghi/GIDEON/blob/main/docs/brand/BRAND.md)),
which draws every figure a second time as `<name>-dark.png` for the dark theme; and `gideon_palette.py`, the grounds
and the colours of the data. The Python needs numpy and matplotlib (`GIDEON_DOCS_PYTHON`, else
`exe/.docs-venv/bin/python`, else `python3`).

What the pages say about the code is checked too: `scripts/check_docs.py` fails when a public entity, a module, a test
program, an example, a fobos mode, a script or a figure is named on no page, when a page names a `gideon_*` procedure
or an ADAM routine that does not exist, or when an image, a link, an include or a snippet does not resolve. The code of
the tutorials and the cookbook is compiled and run by `scripts/run_examples.sh`, each output compared with the one the
page shows.
