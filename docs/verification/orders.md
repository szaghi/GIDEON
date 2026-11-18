# A1 Orders of accuracy

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| A1 | Accuracy | P | CPU (gfortran R4P, R8P, R16P; nvfortran R8P); the CI matrix and the CI's quadruple-precision run | yes | `gideon_test_order` |

## Case and scope

Every compiled kernel must converge at its design order on a smooth function. C1 measures the rows at 100 digits;
this case measures the kernels as compiled, in the working kind, where round-off competes with truncation. The
question is twofold: does every kernel reach its order where the arithmetic can show it, and where can it? The
answer to the second is the reason the orders are verified in quadruple precision.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Function | $g(x) = \sin(1.3x + 0.4) + 0.5\cos(0.7x - 0.2)$ at $x_0 = 0.37$: two incommensurate waves at a generic phase, no derivative vanishes at the point, so no leading term cancels by accident | [`gideon_test_functions.F90`](https://github.com/szaghi/GIDEON/blob/main/src/tests/gideon_test_functions.F90) |
| Data, target | point values, or exact cell averages (differences of the primitive over $h$); the exact derivative, or its exact cell average | `data_generic`, `target_generic` |
| Steps | $h = 2^{-1} \dots 2^{-20}$, as long as $1/h^{n+1}$ fits the kind | [`gideon_test_order.F90`](https://github.com/szaghi/GIDEON/blob/main/src/tests/gideon_test_order.F90) |
| Round-off floor | $\varepsilon\,[\sum_k\vert c_k\vert(\max\vert q\vert + 2/h)/h^n + \vert T\vert + 2/h]$, the cell-average terms only for averaged data or target: an upper bound of the rounding of the data, the dot product and the target | |
| Resolved step | error above 10 floors, and the window shorter than a twelfth of the shorter wavelength | |
| Observed order | $p = \log_2(e(h)/e(h/2))$ between resolved steps; a plateau is three consecutive resolved steps whose two orders differ by at most 0.25 | |
| Asserted (frozen) | on the finest plateau both orders $\ge$ design $- 0.25$; a kernel with no plateau is *unresolved in this precision*, not failed | |

## Program

```bash
./exe/gideon_test_order          # R8P
./exe-r16/gideon_test_order      # R16P: every kernel measured
GIDEON_REPORT_DIR=dir ./exe-r16/gideon_test_order   # also order-r128.csv, order-summary-r128.csv in dir
```

## Results

![e(h) of four kernels in R8P and R16P with their round-off floors](/verification/convergence.png)

*Four kernels, R8P (umber) and R16P (blue): filled markers are resolved steps, hollow ones are not; the dotted lines
are the round-off floors. In R8P the floor of `fd_d1_c10`, `fv_face_ur9` and `fv_d6_ur4` meets the truncation error
before three resolved steps line up: no plateau. In R16P the floor drops with ε and every kernel has
one.*

<!--@include: ./_tables/orders-builds.md-->

By family: how many kernels each kind measures, and the smallest margin of the quadruple-precision orders:

<!--@include: ./_tables/orders-families.md-->

![Observed minus design order of every kernel](/verification/orders.png)

*The smaller of the two observed orders on the finest plateau, minus the design order, for every kernel: all 332 in
R16P, the measured ones in R8P. The dashed line is the frozen bound.*

<details>
<summary>Every kernel: the two orders of its finest plateau, R8P and R16P, and the generator's 100-digit order</summary>

<!--@include: ./_tables/orders-kernels.md-->

</details>

## Findings

- Every kernel reaches its design order: in R16P all of them, on a plateau, within the bound; in R8P and R4P every
  kernel that has a plateau.
- Double precision cannot verify the wide stencils. Whether a step is resolved depends on the wavenumber times $h$
  only, and a wide stencil leaves its pre-asymptotic range only at steps where its round-off floor,
  $\propto \varepsilon/h^n$, is already above its truncation error. A different function does not help; a wider
  kind does. This is why `tests-gnu-r16` is a test mode and why CI runs it.
- The R8P plateaus scatter more than the R16P ones: they sit at coarser steps, where the next term of the error is
  not yet negligible. The bound, design $- 0.25$, holds for them all.
- R4P resolves only a handful of kernels: a single-precision build is verified by C2, C3 and the golden values,
  not by its orders ([Single and quadruple precision](/cookbook/precision)).
