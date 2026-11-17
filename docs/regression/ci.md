# Continuous integration

Every push to any branch, and every pull request, runs the suite on GitHub Actions. The workflows are thin callers of
the reusable workflows of [fortran-toolchains](https://github.com/szaghi/fortran-toolchains), which build and test in
container images with the compilers installed; they are managed by `fobis scaffold` from the `[scaffold]` section of
`fobos`. A newer push to the same ref cancels the running checks.

| Workflow | When | What runs | Where |
|---|---|---|---|
| [`ci.yml`](https://github.com/szaghi/GIDEON/blob/main/.github/workflows/ci.yml) | push, pull request | `fobis fetch`, then `fobis rule --ex makecoverage-analysis`: the tests with coverage, the coverage badge, the R16P tests, the generator check, the examples, the docs check | `linux-gnu` image |
| [`matrix.yml`](https://github.com/szaghi/GIDEON/blob/main/.github/workflows/matrix.yml) | push, pull request | [`scripts/ci/matrix.sh`](https://github.com/szaghi/GIDEON/blob/main/scripts/ci/matrix.sh): `fobis fetch`, `fobis build --mode tests-$TOOLCHAIN_FAMILY`, `bash scripts/run_tests.sh` | `linux-gnu`, `linux-nvhpc`, `linux-oneapi` images |
| [`fpm.yml`](https://github.com/szaghi/GIDEON/blob/main/.github/workflows/fpm.yml) | push, pull request | `fpm test` and `fpm build --profile release` (the tests listed in `fpm.toml`) | `linux-gnu` image |
| [`docs.yml`](https://github.com/szaghi/GIDEON/blob/main/.github/workflows/docs.yml) | a successful CI run on the default branch, or by hand | the formal API and the VitePress site, with CI's coverage badge, deployed to GitHub Pages | `linux-gnu` image |
| [`release.yml`](https://github.com/szaghi/GIDEON/blob/main/.github/workflows/release.yml), [`install.yml`](https://github.com/szaghi/GIDEON/blob/main/.github/workflows/install.yml) | a `vX.Y.Z` tag | the GitHub release, then the release installed on a fresh runner as its users do (`scripts/install.sh`, fpm) | runner |

## The coverage run and the quadruple precision

`ci.yml` runs one FoBiS rule, the same that runs locally:

```ini
[rule-makecoverage-analysis]
rule_1 = fobis rule --ex makecoverage           # clean, build tests-gnu-coverage, run the suite
rule_2 = bash scripts/compute-coverage.sh       # gcov over the .gcda files: docs/public/coverage.json
rule_3 = fobis build --mode tests-gnu-r16       # then the tests in quadruple precision
rule_4 = bash scripts/run_tests.sh --exe-dir exe-r16
rule_5 = python3 scripts/gen_gideon.py --check  # the rows (C1) and the generated sources up to date
rule_6 = bash scripts/run_examples.sh           # the tutorials' and the cookbook's programs, R8P and R4P outputs
rule_7 = python3 scripts/check_docs.py          # the pages against the code
```

The coverage is the line coverage of the library and the tests summed by `gcov` (`scripts/compute-coverage.sh`); the
badge of the site reads the `coverage.json` of the last run on the default branch. The quadruple-precision build is in
the same rule because it is where the orders of every kernel are observed ([A1](/verification/orders)): a change
that breaks the order of a wide stencil passes every R8P test and fails here. The last three rules hold the
documentation to the code: a row or a generated source edited by hand, an example whose output changed, a page that
names a procedure the library does not have, or a public entity no page names, fails CI
([From the data to the pages](/verification/#from-the-data-to-the-pages)).

## The matrix and its floating-point flags

The matrix builds `tests-gnu`, `tests-nvf` or `tests-intel` (the family of the image) and runs the whole suite,
including the bitwise contracts: the mirrors and symmetries ([C3](/verification/properties)), the dispatch
([C4](/verification/dispatch)) and the bitwise equality with ADAM ([C5](/verification/adam)). Two compilers
need a flag for them to hold, and the test modes set it (`fobos`, `[common-variables]`):

```ini
# ifx's default fp-model (fast=1) reassociates inlined arithmetic: the pair form's exact symmetries and the bitwise
# equality with ADAM (7 curl values 1 ulp off) need the operations as written
$FP_INT  = -fp-model=precise
# nvfortran's -fast vectorises for the host it runs on: on another CPU (the CI runner) it reorders the longer sums and
# the bitwise equality with ADAM breaks at s = 4, 5. The test modes add -Kieee (strict IEEE); the libraries and the
# benchmarks keep -fast alone.
$IEEE_NVF = -Kieee
```

- **ifx, `-fp-model=precise`**, in every Intel mode, libraries included: GIDEON's results are defined by the order
  of the operations in the generated code. With the default model, ifx reassociates the pair form after inlining, and
  the bitwise contracts break.
- **nvfortran, `-Kieee`**, in the test modes only: the libraries and the benchmarks keep `-fast`, whose
  vectorisation depends on the CPU that compiles. On the CI runner, a CPU other than the developer's, the equality
  with ADAM broke at the widest stencils (the record in `fobos`, above). The tests use strict IEEE so that what they
  assert does not depend on the machine.
- **gfortran** needs no flag: `-O3` (the release profile) keeps the operations as written.

## What CI does not run

- the GPU: no runner has one. The device test runs locally (`tests-nvf-oac`, [C6](/verification/device));
- the R4P builds (`tests-*-r4`): locally, with their golden file ([Golden values](./golden));
- the mutation run, which rebuilds the suite once per mutant ([Mutation testing](./mutation));
- `scripts/run_report.sh`, which regenerates the data, tables, figures and casts of these pages: run it after a
  change of the generator or the kernels, and commit what it writes.
