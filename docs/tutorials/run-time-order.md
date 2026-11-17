# T4 The order at run time

**Goal.** Choose the order when the program runs: what a dispatcher does, where its window bounds come from, what it
returns for an order without a kernel, what it costs against the kernel, and how to avoid that cost.

**Concepts.** [Run-time order](/guide/conventions#run-time-order) and the [API](/guide/api#run-time-order); the
dispatcher of every kernel is in the last column of the table in [Operators](/guide/operators#every-kernel). The
measured cost: [Performance](/guide/performance#the-run-time-dispatcher-budget-missed). T1 for the kernels.

`gideon_fd_d1_c(order, q, ds, d)` is one `select case` on `order` that calls `gideon_fd_d1_c<order>`. Its dummy is
declared with bounds that are expressions of `order`, so the caller passes the whole window of that order and the
kernel receives it with the offsets it expects:

| Dispatcher | Window | Orders with a kernel |
|---|---|---|
| `gideon_fd_d1_c` | `q(-(order/2):order/2)` | 2, 4, 6, 8, 10 |
| `gideon_fd_d3_c` | `q(-(order/2+1):order/2+1)` | 2, 4, 6, 8 |
| `gideon_fd_face_c` | `q(1-(order/2):order/2)` | 2, 4, 6, 8, 10 |
| `gideon_fv_face_ur` | `q(1:order)` | 1 to 9 |
| `gideon_fv_face_ul` | `q(-(order-1):0)` | 1 to 9 |
| `gideon_fv_d1_ur` | `q(0:order)` | 1 to 9 |
| `gideon_fv_d1_ul` | `q(-(order):0)` | 1 to 9 |

## The program

<<< @/examples/src/run_time_order.F90{fortran}

- **Orders 1 to 12.** The same loop for every order: `s = order/2`, the section `u(i-s:i+s)`. Six ghost nodes cover
  the widest window tried, order 12.
- **An order without a kernel.** The dispatcher writes `GIDEON_NAN`, a quiet NaN, into the result. The program tests
  it with `ieee_is_nan`; `d /= d` is the same test without the module, but compilers may fold it away under
  value-unsafe optimisations (`-ffast-math`).
- **The dispatcher is the kernel.** Order 6 through `gideon_fd_d1_c` and through `gideon_fd_d1_c6`: compared bit for
  bit.
- **The choice made once.** A `select case` outside the loop, one loop per order with the kernel inside: the same
  values, and no branch per point.

## Build and run

```bash
fobis build --mode examples-gnu && ./exe-examples/run_time_order
```

or against the library:

```bash
fobis build --mode static-gnu
gfortran -I lib/gnu/static/mod docs/examples/src/run_time_order.F90 lib/gnu/static/libgideon.a -o run_time_order
./run_time_order
```

## Output

<<< @/examples/output/run_time_order.txt{text}

## What to notice

- **The window follows the order, the result follows the kernel.** Orders 2 to 10 give the errors of T1 at 64 nodes
  (6.512E-03 for order 2, 1.577E-06 for order 6), the kernels' own numbers.
- **An odd order has a window but no kernel.** Order 3 passes `q(-1:1)`, a legal window; the `select case` has no
  branch for it and returns `GIDEON_NAN`, as for order 12. No `error stop`: the dispatchers are device code, where it
  is not allowed, and ADAM's silent 0 (audit item C) looked like a valid derivative. A NaN propagates through every
  value computed from it, so one check on the result of a step finds it.
- **The dispatcher costs a branch per call.** Measured against `gideon_fd_d1_c6` on $2^{24}$ points: 1.093 times the
  kernel's time on an RTX 4070, 1.183 with gfortran on one core, 1.115 with nvfortran
  ([Performance](/guide/performance#the-run-time-dispatcher-budget-missed)). The compilers do not hoist the
  loop-invariant `select case` out of the loop, and on the CPU it blocks the loop's vectorisation. The last block of
  the program is the remedy: the same numbers, the branch taken once. The timings are not printed: they depend on the
  machine (exercise 1).

## Exercises

1. Time the dispatcher loop and the hoisted one on $2^{22}$ nodes with `system_clock` (not `cpu_time`), five
   repetitions each, and compare the medians.
2. Upwind faces at run time: `gideon_fv_face_ul(order, u(i-order+1:i), f)` for orders 1 to 10, on cell averages
   (T3). Which order returns `GIDEON_NAN`?
3. Read the order from the command line (`get_command_argument`), the reason a run-time order exists.
4. `gideon_fd_d3_c` reads a wider window than `gideon_fd_d1_c` of the same order. Write its section for node `i`
   from the table above and check it against `gideon_fd_d3_c6`.
