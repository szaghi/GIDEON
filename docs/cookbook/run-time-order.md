# A run-time order, fast

**Problem.** The order of the operators is an input of the solver (PRISM's `fdv_order`), not a compile-time constant,
and the loops over the grid must not pay for that choice at every point.

## What the dispatcher costs

Every 1-D operator has a run-time order dispatcher, `gideon_<data>_<op>_<bias>(order, q, [ds,] r)`: one
`select case(order)` that calls the kernel of that order, and returns `GIDEON_NAN` for an order without one. It is
device code like the kernels (`!$acc routine seq`), uniform across the points when the order is, so it does not
diverge on a GPU. It is not free. Measured against the kernel (`gideon_fd_d1_c(6, ...)` against `gideon_fd_d1_c6`,
$2^{24}$ points, the median of four runs, `notes/bench.md`):

| | RTX 4070, nvfortran OpenACC | one core, gfortran | one core, nvfortran |
|---|---|---|---|
| dispatcher / kernel | 1.093 (1.093–1.104) | 1.183 (1.134–1.220) | 1.115 (1.102–1.136) |

The `select case` is evaluated at every call: no compiler hoisted it out of the loop, although the order is
loop-invariant, and on the CPU it also blocks the vectorisation of the loop. The dispatcher is still cheaper than
ADAM's run-time-`s` path, about 0.45 to 0.65 of it (1.093 × 0.532 on the GPU, 1.183 × 0.425 with gfortran), but
the 9 to 18 % are avoidable.

## Recipe: validate once, choose once

1. **Validate the order before the loop.** Call the dispatcher once on a probe window of that order; a NaN result
   means no kernel exists (`ieee_is_nan` from `ieee_arithmetic`). Inside a device loop there is no way to stop, so a
   bad order caught there is a field of NaN.
2. **Hoist the choice.** One `select case` outside, one loop per order inside, each calling its kernel: the loop body
   is the kernel's, with its compile-time row and its vectorisation.
3. **Keep the dispatcher** where the loop body is not the hot spot, or where one loop per order would multiply a long
   body: it is the same arithmetic, bitwise.

A procedure pointer chosen once is not a substitute: the kernels of different orders have different window shapes, so
they share no interface, and a procedure pointer is not device code.

## Code

<<< @/examples/src/cookbook_run_time_order.F90{fortran}

```text
<!--@include: ../examples/output/cookbook_run_time_order.txt-->
```

Order 3 has no centred kernel and is refused before any loop runs; every other order gives the same values through
the dispatcher and through the hoisted kernels, to the last bit (`gideon_test_dispatch` asserts the same for all 288
kernels and their 54 dispatchers). On a GPU the hoisted form puts the `!$acc parallel loop` on each loop inside the
cases, so every case launches its own kernel with its own compile-time row ([On the GPU](/tutorials/gpu)).

**Reference.** The measurement: [Performance](/guide/performance); the dispatchers' NaN contract:
[Run-time dispatch](/verification/dispatch); every dispatcher with its window: [API](/guide/api).
