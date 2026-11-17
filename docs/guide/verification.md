# Verification

What each check proves, where it lives, and what it measured on 2026-11-17 (gfortran 16.2 and nvfortran 26.9, debug
and optimised, R8P and R4P; the device test on an RTX 4070).

## In the generator (`scripts/gen_gideon.py --check`)

| Check | Result |
|---|---|
| Every row solved exactly from its moment conditions; its exact degree measured beyond them | 243 rows |
| ADAM's committed tables parsed and compared with the derived rows | 93 of 93 equal |
| Left-side upwind rows solved independently and compared with the mirrored right ones | all mirrors |
| Order of every row on a generic smooth function, 100-digit decimals, h = 2⁻¹¹ and 2⁻¹² (bound: design − 0.25) | 243 of 243, smallest margin −0.002 |
| The emitted Fortran rows read back and compared with the derived ones | equal |

## In the test suite (`scripts/run_tests.sh`)

| Program | Gate | What it checks |
|---|---|---|
| `gideon_test_tables` | G2 | every compiled kernel exact on xʲ up to its degree, and on x^(K+1) its exact leading error, within 64 eps, in every kind |
| `gideon_test_order` | G3 | observed orders on the generic function where the precision resolves a plateau (R8P: 96 kernels; R4P: 11; the rest covered by the generator's check) |
| `gideon_test_properties` | G4 | the mirror identity, the (anti)symmetries, constant data, the scalings: all bitwise |
| `gideon_test_dispatch` | — | every dispatcher equals its kernel bitwise; NaN for an order without a kernel |
| `gideon_test_operators_3d` | — | the 3-D operators exact on random polynomial fields of their degree, FD and FV, every order |
| `gideon_test_adam` | G6 | ADAM's module compiled verbatim: 279 values bitwise, 20 within 64 eps, ADAM's defects reproduced (8), its zero rows (43) |
| `gideon_test_golden` | G7 | 1944 values against `src/tests/golden/gideon_golden_v1.txt` (R4P: `_r4_v1`) |
| `gideon_test_device` | G8 | all 243 kernels on the GPU against the host: 15.9 M values, 76 % bitwise, worst 4.3e-16 |
| debug modes | G5 | every test under `-fcheck=bounds`: a wrong window is an error |

## Mutation testing (`scripts/mutation_test.py`)

22 mutants (a coefficient of each row family, a pair sign, a window offset, the division by the step, the mirror and
its sign, the order of a flux difference, a dispatcher's case and its NaN, the 3-D wiring): 22 killed
(`notes/mutation-run-2026-11-17.txt`).

## Why the orders are measured in the generator

In double precision the widest stencils (upwind, 7 to 9 cells) leave their pre-asymptotic range only where round-off
already dominates. The ratio of truncation to round-off depends on the wavenumber times the step alone, so no choice
of test function changes it. The generator measures every row in 100-digit arithmetic. The Fortran test measures the
compiled kernels where a plateau exists, and the exactness test proves the compiled rows in every kind.
