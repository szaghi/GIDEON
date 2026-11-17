# C2 Exactness and leading error

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| C2 | Contracts | C | CPU (gfortran R4P, R8P, R16P; nvfortran); the CI matrix | no | `gideon_test_tables` |

## Case and scope

C1 checks the rows as fractions; this case checks the compiled kernels, in the working kind, against the same
meaning. A kernel of exact degree $K$ must return the exact target on every monomial $x^0 \dots x^K$, and on
$x^{K+1}$ it must miss by exactly the leading error the generator derived, not by more and not by less. The second
half catches what exactness alone would not: a row that is exact but wrong beyond its degree (a different stencil
of the same order), or a compiled coefficient off by more than its rounding.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Data | the monomial $x^j$ on the kernel's window ($h = 1$, $x_0 = 0$): point values, or exact cell averages $(b^{j+1} - a^{j+1})/(j+1)$ | [`gideon_test_functions.F90`](https://github.com/szaghi/GIDEON/blob/main/src/tests/gideon_test_functions.F90) `data_monomial` |
| Target | the n-th derivative at 0, or its average over cell 0 | `target_monomial` |
| Kernels | all of them, through the registry's `apply` (one `select case` over the 288 kernels) | [`gideon_test_registry.F90`](https://github.com/szaghi/GIDEON/blob/main/src/tests/gideon_test_registry.F90) (generated) |
| Asserted, $j \le K$ | $\vert k(x^j) - T(x^j)\vert \le 64\,\varepsilon \sum_k \vert c_k\vert \max(\vert q\vert, 1)$ | [`gideon_test_tables.F90`](https://github.com/szaghi/GIDEON/blob/main/src/tests/gideon_test_tables.F90) |
| Asserted, $j = K + 1$ | $\vert (k - T) - E\vert \le$ the same tolerance, $E$ the exact leading error (`kerror`, derived by `leading_error` in the generator) | |

The tolerance, 64 ε of the scale, was frozen before the kernels were written: the round-off of a dot product of a
few terms, with margin. In R4P it also covers the two coefficients whose literal rounds before the division
([Precision](/guide/precision)).

## Program

```bash
./exe/gideon_test_tables                 # R8P
./exe-r4/gideon_test_tables              # R4P
./exe-r16/gideon_test_tables             # R16P
```

## Results

<!--@include: ./_tables/exactness-builds.md-->

The checks are $K + 2$ per kernel ($x^0 \dots x^{K+1}$); the count is recomputed from the generator's rows when this
table is written, and must match the program's.

<details>
<summary>The exact leading error of every kernel</summary>

<!--@include: ./_tables/exactness-leading.md-->

</details>

## Findings

- Every kernel is exact to its degree and misses the next monomial by the derived amount, in every kind: the
  compiled rows are the derived ones, neither less nor more exact.
- The leading error is the constant of the order: at step $h$ the error of a kernel on a smooth $u$ is
  $E\,h^{p}\,u^{(K+1)}/(K+1)!$ to leading order. Two kernels of the same derivative and order compare by
  $\vert E\vert$, and the one-sided rows pay for their window ([Choose an operator](/cookbook/choose)):

<!--@include: ./_tables/exactness-constants.md-->
