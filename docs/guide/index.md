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
order dispatchers, 40 3-D operators and their 8 dispatchers. Phase 2 adds what ADAM's callers write by hand, on point
values: the face derivatives d1–d5 (the viscous fluxes), and the face value from one side (extrapolation) or
upwind-biased (interpolation): 288 kernels, 54 run-time order dispatchers. Bitwise equal to ADAM where ADAM is right, at 0.17–0.82
of its cost on an RTX 4070 and on one CPU core ([Performance](/guide/performance)). Tested with gfortran 16.2,
nvfortran 26.9 (also on the GPU, OpenACC) and ifx, in single and double precision, and with gfortran in quadruple.

## How to read this documentation

The sidebar is the reading order, and the **Next** link at the foot of every page follows it: read straight through,
or stop where you have what you need.

1. **Start here**: [installation](/guide/installation), a [quick start](/guide/quickstart) (one derivative, compiled
   and run), and the [conventions](/guide/conventions): how a kernel is named and which window of your array it
   reads. These four pages are all a first use needs.
2. **[Tutorials](/tutorials/)**: six complete programs, from a first derivative to the GPU, each compiled and run
   with the output on the page. Read them in order; each adds one idea.
3. **[Cookbook](/cookbook/)**: one recipe per task (choose an operator and its order, size the ghost cells, convert
   averages and point values, a fast run-time order, another precision, replace ADAM's calls, add an operator).
   Read the one you need.
4. **Reference**: every [operator](/guide/operators) with its window and order, the [API](/api/),
   [precision](/guide/precision), [the generator](/guide/generator) and [performance](/guide/performance); the
   module pages generated from the sources are under [API](/api/) in the top bar.
5. **[Verification](/verification/)** and **[Regression](/regression/)**: what is proved about every kernel and how,
   one page per test, and what keeps it from changing.
