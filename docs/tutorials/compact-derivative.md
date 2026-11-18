# T7 A compact derivative on a periodic line

**Goal.** The first derivative of $u = e^{\sin x}$ on a periodic line with the compact (Padé) operators of order 4,
6, 8 and 10, against the explicit centred operator of order 6, and the orders measured by the program itself.

**Concepts.** A compact operator is a band solved along the line, its right side an ordinary kernel
([Compact operators](/guide/compact)): `gideon_fd_d1_p6` returns the right side at one node, from the window
`q(-2:2)`; the band `FD_D1_P6_LEFT` ($[1/3, 1, 1/3]$, stored as pairs: index $j$ multiplies the unknowns at $\pm j$) is
solved by `gideon_tri_periodic_const`. [Line operators](/guide/conventions#line-operators).

## The program

<<< @/examples/src/compact_derivative.F90{fortran}

- **The right side.** One kernel call per node, exactly as for an explicit derivative ([T1](./fd-derivative)): the
  window by offset, the step an argument, the ghost nodes the periodic image. The widest right side, `p10`, reads
  `q(-3:3)`: three ghosts on each side.
- **The solve.** The array of right sides becomes the derivative in place: `gideon_tri_periodic_const(n, a, b, c, d,
  w, z)` solves $a\,x_{i-1} + b\,x_i + c\,x_{i+1} = d_i$ with indices modulo $n$, the band passed as scalars (the
  `_const` variant), two scratch arrays the caller's. Orders 8 and 10 have pentadiagonal bands:
  `gideon_penta_periodic_const`, four scratch arrays.
- **Nothing else changes.** The kernels and the solvers are `pure` and device-callable, so the same two calls run one
  line per thread on a GPU ([T6](./gpu), [C8](/verification/solvers)).

## The output

```text
<!--@include: ../examples/output/compact_derivative.txt-->
```

## What it shows

- **The orders.** 4.01, 6.02, 8.02 at 128 nodes; `p10` measures 10.05 at 64 and then meets round-off (its error at 128
  nodes is $2.6 \cdot 10^{-14}$). The orders of every compact operator are measured in quadruple precision by
  [A4](/verification/lines).
- **The accuracy at equal order.** `p6` reads five points where the explicit `c6` reads seven, and its error is about
  1/15 of `c6`'s on the same grid ($1.7 \cdot 10^{-9}$ against $2.5 \cdot 10^{-8}$ at 128 nodes): the band buys the
  resolution of the short waves. `p4`, with three points, is more accurate than
  `c6` on the coarsest grid (16 nodes, where the short waves dominate), less on the finer ones, where order wins.
- **The cost.** A solve per line. Here each grid has one line, so `_const` is the natural call; when many lines share
  the band and the length (every line of a grid direction), factor the band once with `gideon_tri_periodic_factor`
  and solve each line with `gideon_tri_periodic_solve`, multiplications only: `p6` then costs about 3 times `c6` per
  point on a CPU core, 5 times on a GPU ([Compact operators](/guide/compact), [Performance](/guide/performance)).

**Next.** A compact filter after each step of a central scheme: [Filtering a central scheme](/cookbook/filters);
compact rows coupled with a nonlinear scheme at a step: [Compact-WENO rows](/cookbook/compact-weno).
