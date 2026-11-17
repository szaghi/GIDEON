# C5 Equivalence with ADAM

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| C5 | Contracts | C | CPU, R8P only (gfortran; nvfortran); the CI matrix | yes | `gideon_test_adam` |

## Case and scope

GIDEON replaces ADAM's FDV operators. Where ADAM is right, a GIDEON kernel must return ADAM's value, bitwise, so
that ADAM's results do not move when it switches; where ADAM is wrong, GIDEON must be right, and the difference must
be the documented defect, not another one. The audit of ADAM's operators (`notes/adam-fdv-audit.md`, ADAM `ff879596`)
found six items:

| Item | ADAM's defect | How GIDEON avoids it | Checked here |
|---|---|---|---|
| A | `compute_derivative1_fv_lupwind` and `compute_reconstruction_r_fv_lupwind`, $s \ge 2$: the left-side row is mirrored twice; first order instead of order $s$ | left rows are the right rows read mirrored, never stored ([C3](./properties)) | ADAM misses a linear (face value) or quadratic (d1) function that GIDEON reproduces exactly |
| B | `compute_derivative1_fv_rupwind`: its window `q(0:)` over `[0:1+s]` and ADAM's own test disagree; the test reads one element out of bounds | explicit-shape windows by offset, one convention per family | the right-side d1 on ADAM's documented window, cell at `q(1)`: bitwise |
| C | rows ADAM has not got (FD d3–d6, FV face d2–d5 and cell d3–d6 at small $s$; upwind d6) return 0 | absent kernels; a dispatcher returns NaN ([C4](./dispatch)) | ADAM returns exactly 0 |
| D | ADAM's FDV tests print errors, they do not assert | every GIDEON test asserts a frozen bound | (the test suite itself) |
| E | ADAM's generator validates its tables against the committed tables, not against their meaning | rows derived by meaning ([C1](./rows)) | (C1) |
| F | the upwind names flip: ADAM's `rupwind` d1 reads the right side, its `rupwind` d2–d5 the left | `ur`, `ul` name the side the stencil lies on, for every operator | each ADAM upwind operator mapped to the GIDEON side it really reads |

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| ADAM's code | `adam_fdv_operators_library.F90` and `adam_fnl_fdv_operators_library.F90`, byte copies of `notes/adam-reference/`, compiled into the R8P test builds only | [`src/tests/`](https://github.com/szaghi/GIDEON/tree/main/src/tests) |
| Data | random numbers in $[-1, 1)$, a fixed seed; $\Delta s = 0.37$, $\Delta xyz = (0.37, 0.53, 0.29)$; every half-stencil $s = 1 \dots 5$ | [`gideon_test_adam.F90`](https://github.com/szaghi/GIDEON/blob/main/src/tests/gideon_test_adam.F90) |
| Index maps | ADAM's node or cell sits at `q(1)`, GIDEON's at offset 0; faces at the same indices | |
| Bitwise | `a == g` | `same` |
| Within tolerance | $\vert a - g\vert \le 64\,\varepsilon \cdot 10^3 \sum\vert q\vert/\Delta s^n$: ADAM sums the right-side upwind rows farthest cell first, GIDEON nearest first | `close` |
| Defect | GIDEON within $640\,\varepsilon$ of the exact value, ADAM off by more than $10^{-3}$ | `defect` |
| Absent | ADAM returns exactly 0 | `zero` |
| 3-D operators | ADAM's FNL (directional stencils) and CPU (cube) signatures against GIDEON's directional one | `three_d` |

## Program

```bash
./exe/gideon_test_adam
```

## Results

![Comparisons with ADAM by operator group and class](/verification/adam.png)

*Every comparison of the test by ADAM operator group: bitwise, within the summation-order tolerance, ADAM's defect
reproduced (GIDEON exact), ADAM's zero row.*

<!--@include: ./_tables/adam-groups.md-->

<!--@include: ./_tables/adam-builds.md-->

## Findings

- GIDEON is bitwise equal to ADAM on every centred operator, the deconvolution, the 3-D operators in both of ADAM's
  signatures, and the upwind operators that ADAM sums in GIDEON's order: switching ADAM to GIDEON does not move a
  bit of its centred results.
- The right-side upwind d2–d5 and face derivatives differ by the summation order only, within the tolerance; the
  largest differences are in the table.
- ADAM's left-side face value and d1 are wrong for $s \ge 2$ (A), at every $s$ the test runs; for $s = 1$ the double
  mirror is the identity and they are right (bitwise).
- ADAM's right-side d1 is right on its documented window `[0:1+s]` with the cell at `q(1)`; B is a mismatch between
  that window and ADAM's own test, not a wrong row.
- The mapping of F is in the table: ADAM's `rupwind` d2–d5 and `rupwind_d` are GIDEON's left side, its `lupwind` the
  right side ([Replace ADAM's FDV calls](/cookbook/from-adam)).
