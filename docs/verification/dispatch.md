# C4 Run-time dispatch

| ID | Group | Type | Backends | Figure | Program |
|---|---|---|---|---|---|
| C4 | Contracts | C | CPU (gfortran R4P, R8P, R16P; nvfortran); the CI matrix | no | `gideon_test_dispatch` |

## Case and scope

A solver that reads its order from input (ADAM's PRISM, `fdv_order`) calls a dispatcher,
`gideon_<data>_<op>_<bias>(order, q, [ds,] r)`, which selects the kernel by a `select case` on `order`. Two things
can go wrong: a case calls the kernel of another order or operator, and an order with no kernel returns a number.
ADAM returned 0 for an order it had no row for (audit item C), a plausible value that a caller cannot tell from a
result. A GIDEON dispatcher returns `GIDEON_NAN`, a quiet NaN, for such an order: `error stop` is not device code.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Dispatchers | one per 1-D operator family and bias, generated with their windows as expressions of `order` (each checked by the generator against every kernel of the group) | [`gideon_dispatch.F90`](https://github.com/szaghi/GIDEON/blob/main/src/lib/gideon_dispatch.F90) (generated), `window_expr` in [`gen_gideon.py`](https://github.com/szaghi/GIDEON/blob/main/scripts/gen_gideon.py) |
| Equality | every kernel against its dispatcher at the kernel's order, on the same random window, $\Delta s = 0.37$: bitwise | [`gideon_test_dispatch.F90`](https://github.com/szaghi/GIDEON/blob/main/src/tests/gideon_test_dispatch.F90), registry `apply_dispatch` |
| NaN | every dispatcher at an order without a kernel whose window stays inside the test's data (the largest order plus 2 for a centred family, plus 1 for an upwind one): the result must be a NaN (`r /= r`) | registry `apply_dispatch_bad` |
| 3-D operators | their dispatchers, and the NaN at order 3, are tested by [A2](./operators-3d) | |

The window of a dispatcher grows with `order`, so a caller passes `q(i+lo(order):i+hi(order))` and the dispatcher
hands the whole section to the kernel. Every dispatcher, its window and the order of its NaN check:

<!--@include: ./_tables/dispatch-counts.md-->

<details>
<summary>Every dispatcher</summary>

<!--@include: ./_tables/dispatch-list.md-->

</details>

## Program

```bash
./exe/gideon_test_dispatch
```

## Results

<!--@include: ./_tables/properties-builds.md-->

## Findings

- Every dispatcher returns its kernel's value bitwise, and a NaN for an order without a kernel, in every kind.
- The NaN is a quiet NaN, a constant built by `transfer` from its bit pattern in each kind
  (`gideon_kind_parameters`): no floating-point operation produces it, so a build that traps invalid operations does
  not stop on it. A caller detects it with `r /= r` or `ieee_is_nan`; a NaN that reaches a solution field propagates
  instead of hiding.
- The dispatch is not free: the `select case` is evaluated at every call and the compilers do not hoist it out of a
  loop with a loop-invariant order. Its cost against the kernel is measured in [Performance](/guide/performance); with
  an order known at compile time, call the kernel ([A run-time order, fast](/cookbook/run-time-order)).
