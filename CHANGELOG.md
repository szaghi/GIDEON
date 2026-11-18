# Changelog

All notable changes to this project are documented here.
Versions follow [Semantic Versioning](https://semver.org/).
Format follows [Keep a Changelog](https://keepachangelog.com/).

## [0.3.0] — 2026-11-18
### Added
- **fd**: Face derivatives, one-sided and biased face values from point values

- **dissipation**: Face differences, the upwind-central split and explicit filters

- **solvers**: Tridiagonal and pentadiagonal line solvers, bounded and periodic

- **compact**: Compact first and second derivatives, the band and the right side derived together

- **compact**: Compact filters, their parameter alpha_f at run time

- **compact**: Compact face values and derivatives, centred and upwind-biased

- **solvers**: Masked line solves, the rows of a compact-WENO hybrid replaced

- **solvers**: Factor once, solve many: the constant-band line solves without division


### Documentation
- Phase 2 documented, a recipe for face values from point values

- **api**: One API reference, the generated one, introduced by the overview

- Roadmap #2 documented, the compact operators, the split, the solvers

- The factored solves documented, their mutants, all 56 killed


### Performance
- **bench**: Compact operators against explicit ones of equal order, budget missed

- **bench**: The factored compact solves, budget met


## [0.2.0] — 2026-11-17
### Documentation
- The shell of the refactored documentation, and its audit

- **tutorials**: Six tutorials whose code is compiled and run, the output the program's own

- **verification**: Case pages, figures and casts from report data; the regression section

- **cookbook**: Recipes as compiled programs, and the guide completed

- The pages checked against the code, in CI

- One sidebar in reading order, a static redirect for the moved verification page

- README and landing page in SETH's style

- **brand**: The Ensō to sessen identity, its assets and its generator

- **theme**: The Kogane identity in the site: tokens, hero, navbar, favicons, icons

- **figures**: The dark figures and casts on the smoked-bronze ground

- **readme**: The logo at the head of the README, as SETH's


## [0.1.0] — 2026-11-17
### Added
- Scaffold GIDEON, the FD/FV stencil operators library for HPC

- **generator**: Every stencil row derived by its meaning, in exact arithmetic

- **kernels**: The FD and FV operators of ADAM, generated from the derived rows

- **dispatch**: Run-time order dispatchers, an order without a kernel returns NaN

- **operators**: Gradient, divergence, curl and Laplacian over directional stencils


### Documentation
- **readme**: The authors, and the license in its own section.

- Conventions, operators, API, verification, performance and sources


### Performance
- **bench**: GIDEON's kernels against ADAM's, on the GPU and one CPU core



