# Verification

What each check proves, where it lives, and what it measured on 2026-11-17 (gfortran 16.2, nvfortran 26.9 and ifx
2026.2, debug and optimised, R8P and R4P; gfortran in R16P, the quadruple-precision verification build; the device test
on an RTX 4070). The Intel modes compile with `-fp-model=precise`: ifx's default model reassociates inlined arithmetic,
which costs the bitwise equality with ADAM (7 curl values 1 ulp off) and the exactness of the pair form's symmetries.

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
| `gideon_test_order` | G3 | observed orders on the generic function where the precision resolves a plateau: R16P all 243 kernels; R8P 96; R4P 11 |
| `gideon_test_properties` | G4 | the mirror identity, the (anti)symmetries, constant data, the scalings: all bitwise |
| `gideon_test_dispatch` | — | every dispatcher equals its kernel bitwise; NaN for an order without a kernel |
| `gideon_test_operators_3d` | — | the 3-D operators exact on random polynomial fields of their degree, FD and FV, every order |
| `gideon_test_adam` | G6 | ADAM's module compiled verbatim: 279 values bitwise, 20 within 64 eps, ADAM's defects reproduced (8), its zero rows (43) |
| `gideon_test_golden` | G7 | 1944 values against `src/tests/golden/gideon_golden_v1.txt` (R4P: `_r4_v1`, R16P: `_r16_v1`) |
| `gideon_test_device` | G8 | all 243 kernels on the GPU against the host: 15.9 M values, 76 % bitwise, worst 4.3e-16 |
| debug modes | G5 | every test under `-fcheck=bounds`: a wrong window is an error |

## Mutation testing (`scripts/mutation_test.py`)

22 mutants (a coefficient of each row family, a pair sign, a window offset, the division by the step, the mirror and
its sign, the order of a flux difference, a dispatcher's case and its NaN, the 3-D wiring): 22 killed
(`notes/mutation-run-2026-11-17.txt`).

## Why the orders are observed in quadruple precision

In double precision the widest stencils (upwind, 7 to 9 cells) leave their pre-asymptotic range only where round-off
already dominates. The ratio of truncation to round-off depends on the wavenumber times the step alone, so no choice
of test function changes it. The quadruple-precision build (`tests-gnu-r16`, host only, run by CI after the coverage)
moves the round-off floor from 1e-16 to 1e-34: there the order test observes all 243 compiled kernels. The generator
also measures every row in 100-digit arithmetic, and the exactness test proves the compiled rows in every kind.
