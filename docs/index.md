---
layout: home

hero:
  name: GIDEON
  text: Generated Interpolations and Derivatives, Exact to Order N
  tagline: The FD and FV stencil operators of ADAM, one pure kernel each, exact to order N, on every CPU and GPU. Finite-difference and finite-volume operators for HPC, in modern Fortran.
  image:
    light: /logo.svg
    dark: /logo-dark.svg
    alt: GIDEON, an ensō touched once by its tangent, the G seal in gold
  actions:
    - theme: brand
      text: Quick start
      link: /guide/quickstart
    - theme: alt
      text: Start here
      link: /guide/
    - theme: alt
      text: Tutorials
      link: /tutorials/
    - theme: alt
      text: Cookbook
      link: /cookbook/
    - theme: alt
      text: Verification
      link: /verification/
    - theme: alt
      text: API Reference
      link: /api/
    - theme: alt
      text: View on GitHub
      link: https://github.com/szaghi/GIDEON

features:
  - icon:
      light: /icons/stencils.svg
      dark: /icons/stencils-dark.svg
      width: 44
      height: 44
    title: Every operator, one kernel each
    details: 332 kernels. FD centred face values, node derivatives d1–d6 and face derivatives d1–d5 (the viscous fluxes); FD face values from one side (extrapolation) or upwind-biased (interpolation); FV centred face values, face derivatives d1–d5 and cell-average derivatives d1–d6; FV upwind on right- or left-side stencils of 1 to 9 cells, and upwind-biased; the deconvolution between averages and point values; flux-form dissipation, the upwind–central split, explicit filters. 64 run-time order dispatchers on top. 29 compact operators and their line solvers.
    link: /guide/operators
    linkText: Operators
  - icon:
      light: /icons/derivatives.svg
      dark: /icons/derivatives-dark.svg
      width: 44
      height: 44
    title: Called from inside your loops
    details: Every kernel is pure, acc routine seq and omp declare target; the rows are compile-time constants, so nothing is set up or copied to the device. 0.17–0.82 of ADAM's cost, on an RTX 4070 and on one CPU core.
    link: /guide/performance
    linkText: Performance
  - icon:
      light: /icons/exact.svg
      dark: /icons/exact-dark.svg
      width: 44
      height: 44
    title: Exact coefficients, bitwise ADAM
    details: Every row derived from what it must compute, in exact rationals, its order measured at 100 digits; ADAM's 93 non-zero rows agree. Bitwise equal to ADAM where ADAM is right, its defects (audit items A–F) fixed where it is not.
    link: /verification/adam
    linkText: Equivalence with ADAM
  - icon:
      light: /icons/faces.svg
      dark: /icons/faces-dark.svg
      width: 44
      height: 44
    title: Every upwind operator with its mirror
    details: Gideon's fleece, a test and its reversed control. Each left-side kernel reads the right-side row mirrored, bitwise. The 3-D gradient, divergence, curl and Laplacian, FD and FV, orders 2–10, over directional stencils.
    link: /verification/properties
    linkText: Mirrors and symmetries
  - icon:
      light: /icons/verified.svg
      dark: /icons/verified-dark.svg
      width: 44
      height: 44
    title: Verification you can read
    details: Sixteen programs, also the test suite. Rows by their meaning, exactness to the leading error, mirrors, dispatch, ADAM, host against device, the order of all 332 kernels and all 29 compact operators observed in quadruple precision, the 3-D operators, spectral properties with an advection animation, the upwind–central split, the line solvers against a dense solve, a compact–WENO coupling at a step, golden values.
    link: /verification/
    linkText: Verification
  - icon:
      light: /icons/kernels.svg
      dark: /icons/kernels-dark.svg
      width: 44
      height: 44
    title: Three compilers, CPU and GPU
    details: gfortran, nvfortran and ifx, debug and optimised, plus the OpenACC build on the GPU; single, double and quadruple precision; built with FoBiS or fpm. A mutation run catches 30 of 30 planted defects; the tutorials' code is compiled and run in CI.
    link: /guide/installation
    linkText: Installation
---

<script setup>
import { withBase } from 'vitepress'
import { VPTeamMembers } from 'vitepress/theme'

const authors = [
  {
    avatar: 'https://github.com/szaghi.png',
    name: 'Stefano Zaghi',
    title: 'Chief Yak Shaver & Keeper of the Stencils',
    desc: 'Tests everything twice, the second time upside down.',
    links: [{ icon: 'github', link: 'https://github.com/szaghi' }],
  },
  {
    avatar: withBase('/avatars/claude.svg'),
    name: 'Claude',
    org: 'Anthropic',
    orgLink: 'https://claude.ai',
    title: 'Tireless Rubber Duck & Professional Fleece Wetter',
    desc: "AI pair programmer; derives the coefficients in exact rationals so humans don't have to, and trusts no fraction written in floating point.",
  },
]
</script>

<div class="gideon-home">

## Quick start

<p class="gideon-lead">Build the library and run the verification programs, or pull GIDEON into your own project.</p>

::: code-group

```bash [Build and test]
git clone https://github.com/szaghi/GIDEON && cd GIDEON
fobis fetch                                                # PENF
fobis build --mode static-gnu                              # libgideon.a
fobis build --mode tests-gnu && bash scripts/run_tests.sh  # 16 programs
```

```ini [FoBiS dependency]
# fobos of your project; then: fobis fetch && fobis build
[dependencies]
deps_dir = src/third_party
# GIDEON's own dependency: fobis fetch is not transitive
PENF     = https://github.com/szaghi/PENF :: branch=master
GIDEON   = https://github.com/szaghi/GIDEON :: branch=main
```

```toml [fpm dependency]
# fpm.toml of your project
[dependencies]
GIDEON = { git = "https://github.com/szaghi/GIDEON", branch = "main" }
```

:::

<p class="gideon-lead">Each kernel takes the window of data around one point as an array section, by offset, and returns one number:</p>

```fortran
use gideon
real(R_P) :: u(-4:104), du, uf

call gideon_fd_d1_c6(u(i-3:i+3), ds, du)                   ! d1 at node i, order 6
call gideon_fv_face_ul5(u(i-4:i), uf)                      ! face i+1/2, 5 cells left
call gideon_fd_d1_c(order, u(i-order/2:i+order/2), ds, du) ! order at run time
```

<p class="gideon-lead">The <a href="./guide/quickstart">quick start</a> builds and runs a first program;
<a href="./guide/conventions">conventions</a> explains the names, the windows and the kinds.</p>

## Authors

</div>

<VPTeamMembers size="medium" :members="authors" />

<div class="gideon-home">

## License

<p class="gideon-lead">GIDEON is available under your choice of license.</p>

| License | SPDX | File |
|---|---|---|
| GNU General Public License v3 | `GPL-3.0-or-later` | [LICENSE.gpl3.md](https://github.com/szaghi/GIDEON/blob/main/LICENSE.gpl3.md) |
| BSD 2-Clause | `BSD-2-Clause` | [LICENSE.bsd-2.md](https://github.com/szaghi/GIDEON/blob/main/LICENSE.bsd-2.md) |
| BSD 3-Clause | `BSD-3-Clause` | [LICENSE.bsd-3.md](https://github.com/szaghi/GIDEON/blob/main/LICENSE.bsd-3.md) |
| MIT | `MIT` | [LICENSE.mit.md](https://github.com/szaghi/GIDEON/blob/main/LICENSE.mit.md) |

<p class="gideon-lead">© <a href="https://github.com/szaghi">Stefano Zaghi</a></p>

</div>
