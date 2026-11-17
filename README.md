# GIDEON

**G**enerated **I**nterpolations and **D**erivatives, **E**xact to **O**rder **N**: finite-difference and
finite-volume stencil operators for HPC, in modern Fortran.

[![CI](https://github.com/szaghi/GIDEON/actions/workflows/ci.yml/badge.svg)](https://github.com/szaghi/GIDEON/actions/workflows/ci.yml)
[![Matrix](https://github.com/szaghi/GIDEON/actions/workflows/matrix.yml/badge.svg)](https://github.com/szaghi/GIDEON/actions/workflows/matrix.yml)

GIDEON is the single source of truth for the FDV operators of ADAM: derivatives of any order, interpolations and
extrapolations to faces, on point values and cell averages, centred and upwind. Every coefficient is generated in
exact rational arithmetic and every operator is a `pure`, device-callable kernel with compile-time coefficients, the
conventions of SETH. Gideon's fleece, a test and its mirrored control, is the rule: verification first, and every
left-upwind operator generated as the mirror of its right twin. Tested with gfortran, nvfortran and ifx; OpenACC
builds with nvfortran.

**Inside**: 243 kernels (FD centred face and d1–d6; FV centred face, face d1–d5 and cell-average d1–d6; FV upwind
on right- or left-side stencils of 1 to 9 cells; the deconvolution between averages and point values), the gradient,
divergence, curl and Laplacian in 3-D, and a run-time order dispatcher for each. Bitwise equal to ADAM where ADAM is
right, fixed where it is not (`notes/adam-fdv-audit.md`), at 0.17–0.82 of its cost on a GPU and a CPU core.

```bash
fobis fetch && fobis build --mode tests-gnu && bash scripts/run_tests.sh
```

Documentation: https://szaghi.github.io/GIDEON.

## Authors

**[Stefano Zaghi](https://github.com/szaghi)** · stefano.zaghi@gmail.com
> *Chief Yak Shaver & Keeper of the Stencils*: tests everything twice, the second time upside down

**[Claude](https://claude.ai)** (Anthropic)
> *Tireless Rubber Duck & Professional Fleece Wetter*: AI pair programmer; derives the coefficients in exact rationals
> so humans don't have to, and trusts no fraction written in floating point

## License

GIDEON is available under your choice of license: GPL-3.0-or-later, BSD-2-Clause, BSD-3-Clause or MIT
(`LICENSE.gpl3.md`, `LICENSE.bsd-2.md`, `LICENSE.bsd-3.md`, `LICENSE.mit.md`).
