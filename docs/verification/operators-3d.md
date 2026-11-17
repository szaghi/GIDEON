# A2 3-D operators

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| A2 | Accuracy | P | CPU (gfortran R4P, R8P, R16P; nvfortran); the CI matrix | no | `gideon_test_operators_3d` |

## Case and scope

The gradient, divergence, curl and Laplacian are sums of 1-D kernels over directional stencils, gathered by the
caller (ADAM's FNL signature: `vy_x` is the y component along x). What can go wrong is not a row, already verified,
but the wiring: a component on the wrong axis, a step of another direction, a sign of the curl. This case checks each
operator, FD and FV, at every order, through its run-time dispatcher, against exact derivatives of random
polynomials, with unequal steps so that a swapped step shows. Their equality with ADAM's 3-D operators is in
[C5](./adam).

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Fields | three random polynomials of total degree $p$ (the order), coefficients in $[-1, 1]$, around $(0.3, -0.2, 0.1)$; steps $(0.5, 0.25, 0.75)$ | [`gideon_test_operators_3d.F90`](https://github.com/szaghi/GIDEON/blob/main/src/tests/gideon_test_operators_3d.F90) |
| FD data | point values along each axis | `line` |
| FV data | averages over the 1-D cells along the axis of each stencil (the cells the FV operators difference); targets the averages of the derivatives over that cell | `line`, `dmean` |
| Exactness | a first derivative of order $p$ is exact on degree $p$, a second derivative of order $p$ on degree $p + 1$: both on these fields | |
| Asserted | error $\le 10^4\,\varepsilon\,(1 + \max\vert{\rm data}\vert)/\min(h)^n$, frozen | `check` |
| NaN | order 3 (no kernel) returns NaN from the gradient and the Laplacian, FD and FV | |
| Operators | below, one dispatcher per operator and data | [`gideon_operators_3d.F90`](https://github.com/szaghi/GIDEON/blob/main/src/lib/gideon_operators_3d.F90) (generated) |

<!--@include: ./_tables/operators-3d-list.md-->

## Program

```bash
./exe/gideon_test_operators_3d
```

## Results

<!--@include: ./_tables/operators-3d-builds.md-->

## Findings

- Every operator is exact on polynomials of its degree, with unequal steps, in every kind: the stencils and the steps
  are wired to the right axes.
- The test is exactness, not convergence: the 3-D operators add no arithmetic beyond the 1-D kernels' (a sum per
  component), so their orders are those of [A1](./orders). Three of the mutants of
  [Mutation testing](/regression/mutation) (a curl component, a divergence step, a Laplacian direction) are killed here
  and by C5.
