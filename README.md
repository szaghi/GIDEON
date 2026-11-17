<div align="center">

<picture>
  <source media="(prefers-color-scheme: dark)"  srcset="docs/public/logo-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="docs/public/logo.svg">
  <img src="docs/public/logo.svg" alt="GIDEON, an ensō touched once by its tangent, the G seal in gold" width="320">
</picture>

#### *the FD and FV stencil operators of ADAM, one pure kernel each, exact to order N, on every CPU and GPU*

### **G**enerated **I**nterpolations and **D**erivatives, **E**xact to **O**rder **N**

> Finite-difference and finite-volume stencil operators for HPC, in modern Fortran: derivatives of order 1 to 6,
> values and derivatives at faces, conversions between cell averages and point values, centred and upwind, and the
> 3-D gradient, divergence, curl and Laplacian. The single source of truth for ADAM's FDV operators: one `pure`,
> device-callable kernel per operator and order, coefficients generated in exact rational arithmetic and compiled as
> constants.

[![Version](https://img.shields.io/github/v/release/szaghi/GIDEON?label=version)](https://github.com/szaghi/GIDEON/releases)
[![CI](https://github.com/szaghi/GIDEON/actions/workflows/ci.yml/badge.svg)](https://github.com/szaghi/GIDEON/actions/workflows/ci.yml)
[![Matrix](https://github.com/szaghi/GIDEON/actions/workflows/matrix.yml/badge.svg)](https://github.com/szaghi/GIDEON/actions/workflows/matrix.yml)
[![Fortran](https://img.shields.io/badge/Fortran-modern-734f96?logo=fortran&logoColor=white)](https://fortran-lang.org/)
[![License: GPL v3](https://img.shields.io/badge/license-GPL%20v3-blue)](LICENSE.gpl3.md)
[![License: BSD-2](https://img.shields.io/badge/license-BSD--2--Clause-blue)](LICENSE.bsd-2.md)
[![License: BSD-3](https://img.shields.io/badge/license-BSD--3--Clause-blue)](LICENSE.bsd-3.md)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue)](LICENSE.mit.md)

<div>
<table>
<tr>
<td><b>🧮 Every operator, one kernel each</b><br><sub>243 kernels: FD centred face values and node derivatives d1–d6; FV centred face values, face derivatives d1–d5 and cell-average derivatives d1–d6; FV upwind on right- or left-side stencils of 1 to 9 cells; the deconvolution between averages and point values. 45 run-time order dispatchers on top. <a href="https://szaghi.github.io/GIDEON/guide/operators">Operators</a> · <a href="https://szaghi.github.io/GIDEON/guide/api">API</a></sub></td>
<td><b>⚡ Called from inside your loops</b><br><sub>Every kernel is <code>pure</code>, <code>!$acc routine seq</code> and <code>!$omp declare target</code>; the rows are compile-time constants, so nothing is set up or copied to the device. 0.17–0.82 of ADAM's cost, on an RTX 4070 and on one CPU core. <a href="https://szaghi.github.io/GIDEON/guide/performance">Performance</a></sub></td>
</tr>
<tr>
<td><b>🎯 Exact coefficients, bitwise ADAM</b><br><sub>Every row derived from what it must compute, in exact rationals, its order measured at 100 digits; ADAM's 93 non-zero rows agree. Bitwise equal to ADAM where ADAM is right, its defects (audit items A–F) fixed where it is not. <a href="https://szaghi.github.io/GIDEON/guide/generator">The generator</a> · <a href="https://szaghi.github.io/GIDEON/verification/adam">Equivalence with ADAM</a></sub></td>
<td><b>🪞 Every upwind operator with its mirror</b><br><sub>Gideon's fleece, a test and its reversed control: each left-side kernel reads the right-side row mirrored, bitwise. The 3-D gradient, divergence, curl and Laplacian, FD and FV, orders 2–10, over directional stencils. <a href="https://szaghi.github.io/GIDEON/guide/conventions">Conventions</a> · <a href="https://szaghi.github.io/GIDEON/verification/properties">Mirrors and symmetries</a></sub></td>
</tr>
<tr>
<td><b>🔬 Verification you can read</b><br><sub>Eleven programs, also the test suite: rows by their meaning, exactness to the leading error, mirrors, dispatch, ADAM, host against device, the order of all 243 kernels observed in quadruple precision, the 3-D operators, spectral properties with an advection animation, golden values. <a href="https://szaghi.github.io/GIDEON/verification/">Verification</a></sub></td>
<td><b>🧰 Three compilers, CPU and GPU</b><br><sub>gfortran, nvfortran and ifx, debug and optimised, plus the OpenACC build on the GPU; single, double and quadruple precision; built with FoBiS or fpm. A mutation run catches 22 of 22 planted defects; the tutorials' code is compiled and run in CI. <a href="https://szaghi.github.io/GIDEON/guide/installation">Installation</a> · <a href="https://szaghi.github.io/GIDEON/regression/mutation">Mutation testing</a></sub></td>
</tr>
</table>
</div>

**[Full documentation](https://szaghi.github.io/GIDEON/)**: [start here](https://szaghi.github.io/GIDEON/guide/) · [tutorials](https://szaghi.github.io/GIDEON/tutorials/) · [cookbook](https://szaghi.github.io/GIDEON/cookbook/) · [verification](https://szaghi.github.io/GIDEON/verification/) · [API reference](https://szaghi.github.io/GIDEON/api/)

</div>

---

## Quick start

```bash
git clone https://github.com/szaghi/GIDEON && cd GIDEON
fobis fetch                                                # PENF
fobis build --mode static-gnu                              # libgideon.a
fobis build --mode tests-gnu && bash scripts/run_tests.sh  # 11 programs
```

as a FoBiS dependency, in the `fobos` of your project (then `fobis fetch && fobis build`):

```ini
[dependencies]
deps_dir = src/third_party
# GIDEON's own dependency: fobis fetch is not transitive
PENF     = https://github.com/szaghi/PENF :: branch=master
GIDEON   = https://github.com/szaghi/GIDEON :: branch=main
```

or as an fpm dependency, in your `fpm.toml`:

```toml
[dependencies]
GIDEON = { git = "https://github.com/szaghi/GIDEON", branch = "main" }
```

Each kernel takes the window of data around one point as an array section, by offset, and returns one number:

```fortran
use gideon
real(R_P) :: u(-4:104), du, uf

call gideon_fd_d1_c6(u(i-3:i+3), ds, du)                   ! d1 at node i, order 6
call gideon_fv_face_ul5(u(i-4:i), uf)                      ! face i+1/2, 5 cells left
call gideon_fd_d1_c(order, u(i-order/2:i+order/2), ds, du) ! order at run time
```

The [quick start](https://szaghi.github.io/GIDEON/guide/quickstart) builds and runs a first program;
[conventions](https://szaghi.github.io/GIDEON/guide/conventions) explains the names, the windows and the kinds.

## Authors

**[Stefano Zaghi](https://github.com/szaghi)** · stefano.zaghi@gmail.com
> *Chief Yak Shaver & Keeper of the Stencils* — tests everything twice, the second time upside down

**[Claude](https://claude.ai)** (Anthropic)
> *Tireless Rubber Duck & Professional Fleece Wetter* — AI pair programmer; derives the coefficients in exact rationals
> so humans don't have to, and trusts no fraction written in floating point

## License

GIDEON is available under your choice of license:

| License | SPDX | File |
|---|---|---|
| GNU General Public License v3 | `GPL-3.0-or-later` | [LICENSE.gpl3.md](LICENSE.gpl3.md) |
| BSD 2-Clause | `BSD-2-Clause` | [LICENSE.bsd-2.md](LICENSE.bsd-2.md) |
| BSD 3-Clause | `BSD-3-Clause` | [LICENSE.bsd-3.md](LICENSE.bsd-3.md) |
| MIT | `MIT` | [LICENSE.mit.md](LICENSE.mit.md) |

© [Stefano Zaghi](https://github.com/szaghi)
