# Changelog

All notable changes to this project are documented here.
Versions follow [Semantic Versioning](https://semver.org/).
Format follows [Keep a Changelog](https://keepachangelog.com/).

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



