# GIDEON documentation audit (2026-11-17, before the ADAM wiring)

What the code offers, what the documentation said about it at v0.1.0, and what the refactor adds. The refactor follows
the owner's scientific Fortran projects: SETH (verification case pages, figures and animations from report data,
generated tables, sumi-style theme), SELAH (tutorials whose code is compiled examples with captured output), ADAM
(tutorials, cookbook, concepts/reference split).

## Features (the public surface)

| Feature | Where | Documented at v0.1.0 | Gap |
|---|---|---|---|
| 243 kernels: FD centred face, d1-d6; FV centred face, face d1-d5, cell d1-d6; FV upwind (right/left, 1-9 cells) face, face d1-d5, cell d1-d6; deconvolution | `gideon_fd_centred`, `gideon_fv_centred`, `gideon_fv_upwind`, `gideon_deconvolution` | generated table, conventions | no figure of a window; no worked use; no spectral (dispersion/dissipation) view |
| 45 run-time order dispatchers, NaN for an order without a kernel | `gideon_dispatch` | conventions, API | its cost (1.09-1.18 of the kernel) only in performance; no recipe for hoisting |
| 40 3-D operators (gradient, divergence, curl, Laplacian; FD, FV; c2-c10) + 8 dispatchers | `gideon_operators_3d` | API, operators | no tutorial; the gather of directional stencils not shown |
| The rows as named constants | `gideon_tables` | not documented | reference page |
| Kinds R4P/R8P/R16P, `GIDEON_NAN` | `gideon_kind_parameters` | conventions (short) | precision page: what each kind verifies, R16P host-only |
| The generator: rows by meaning, ADAM cross-check, 100-digit orders, emission | `scripts/gen_gideon.py` | operators (short) | developer page: how to add an operator |
| Tests: tables, order, properties, dispatch, operators_3d, adam, golden, device, kinds | `src/tests/` | one verification page | case pages with the house sections, figures, generated tables |
| Golden files v1 (R8P, R4P, R16P), mutation testing, CI matrix | `src/tests/golden`, `scripts/mutation_test.py`, `.github/workflows` | verification page (lines) | regression section |
| Benchmark against ADAM | `src/bench/`, `notes/bench.md` | performance page (copy of the note) | figures |
| Build: FoBiS modes (gnu/nvf/intel, debug, R4P, R16P, OpenACC, bench), fpm, install.sh | `fobos`, `fpm.toml`, `scripts/install.sh` | installation (short) | every mode, the flags (-fp-model=precise, -Kieee) and why |
| ADAM migration map | `docs/guide/api.md` | table | cookbook recipe with a before/after |

## Sections after the refactor

Guide (about, installation, quickstart, conventions, operators, precision, API, generator, performance, sources),
Tutorials (compiled examples), Cookbook (recipes), Verification (case pages per gate, figures, animations), Regression
(golden, mutation, CI), API reference (formal).

## Closure (end of the night)

Every gap above has its page: windows figures (conventions), the spectral view (A3, with the advection cast), the
dispatcher's cost and its hoisting (cookbook run-time-order), the 3-D operators and the stencil gather (tutorial
operators-3d), the rows (API, generated table), precision (guide), the generator and a new operator (guide, cookbook),
case pages C1-C6 and A1-A3, regression (golden, mutation, CI), performance figure, installation with every mode, and
the ADAM migration recipe. Coverage is now a check, not this table: `scripts/check_docs.py` (in CI's
makecoverage-analysis rule, with `scripts/gen_gideon.py --check` and `scripts/run_examples.sh`).

What no check proves:
- the ADAM "before" calls of cookbook/from-adam: the routine names are checked against ADAM's verbatim copy
  (check_docs, adam), the windows and arguments are not compiled (their equivalence is gideon_test_adam's);
- the prose: that a sentence describes what the code does is read, not tested;
- the GPU tutorial (T6) builds outside the fobos modes (no examples-nvf-oac mode): its output is the host run's.
