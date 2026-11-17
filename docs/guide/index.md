# About

GIDEON (**G**enerated **I**nterpolations and **D**erivatives, **E**xact to **O**rder **N**) is the library of the
linear stencil operators of finite-difference and finite-volume schemes: derivatives of any order, interpolations and
extrapolations to faces, on point values and on cell averages, centred and upwind, and the conversions between
averages and point values.

It is the single source of truth for the FDV operators of ADAM (`adam_fdv_operators_library` and its FNL backend),
rebuilt with the conventions of SETH: one `pure`, device-callable kernel per operator and order, coefficients as
compile-time constants, generated in exact rational arithmetic.

The name comes from Gideon's fleece (Judges 6:36–40): a test, then the same test with the condition reversed, a wet
fleece on dry ground and a dry fleece on wet ground. A test and its mirrored control: verification first, and every
upwind operator generated as the mirror of its twin.

## Status

The first release (v0.1.0) carries ADAM's FD and FV operators, fixed (audit items A–F in `notes/adam-fdv-audit.md`),
with the upwind branch kept and extended to stencils of 9 cells and the sixth derivative: 243 kernels, 45 run-time
order dispatchers, 40 3-D operators and their 8 dispatchers. Bitwise equal to ADAM where ADAM is right, at 0.17–0.82
of its cost on an RTX 4070 and on one CPU core ([Performance](/guide/performance)). Tested with gfortran 16.2,
nvfortran 26.9 (also on the GPU, OpenACC) and ifx, in single and double precision, and with gfortran in quadruple.

## This documentation

| Section | What it holds |
|---|---|
| Guide | [installation](/guide/installation), a [quick start](/guide/quickstart), the [conventions](/guide/conventions) (names, windows, arithmetic), the [operators](/guide/operators), [precision](/guide/precision), the [API](/guide/api), [the generator](/guide/generator), [performance](/guide/performance), [sources](/guide/sources) |
| [Tutorials](/tutorials/) | complete programs, compiled and run, from a first derivative to the GPU |
| [Cookbook](/cookbook/) | recipes: choosing an operator, ghost cells, migrating from ADAM, averages and point values, a run-time order, other precisions, a new operator |
| [Verification](/verification/) | the contracts every kernel meets and the accuracy measured, one page per test |
| [Regression](/regression/) | golden values, mutation testing, continuous integration |
