# T1 A derivative on a periodic grid

**Goal.** The first derivative of $u = e^{\sin x}$ at the nodes of a periodic grid with the centred kernels of order
2, 4, 6 and 8, and the orders measured by the program itself.

**Concepts.** [Names](/guide/conventions#names): `gideon_fd_d1_c6` is the first derivative (`d1`) at a node, from
point values (`fd`), centred (`c`), order 6. [Windows](/guide/conventions#windows): the kernel's dummy is `q(-3:3)`,
the node at offset 0; the caller passes the section `u(i-3:i+3)`. Ghost nodes: [Ghost cells for a
window](/cookbook/ghosts).

## The program

<<< @/examples/src/fd_derivative.F90{fortran}

- **The grid.** $n$ nodes $x_i = i h$, $i = 0, \dots, n-1$, $h = 2\pi/n$. The period holds $n$ distinct nodes: $x = 2\pi$
  is node 0 again.
- **Ghost nodes.** The widest window is `q(-4:4)` (`gideon_fd_d1_c8`): node 0 reads `u(-4:4)`, node $n-1$ reads
  `u(n-5:n+3)`. So `u` carries `G = 4` nodes beyond each end, filled with the periodic image before the loop. The
  kernels never see a boundary, only a window.
- **The call.** `call gideon_fd_d1_c6(u(i-3:i+3), h, du(i,3))`: the section is the window, and node `i` lands at
  offset 0 because the dummy is declared `q(-3:3)`. The step is an argument; the kernel divides by it (decision D4,
  [Arithmetic](/guide/conventions#arithmetic)).
- **The measurement.** The error is the maximum over the nodes against $u' = \cos x\, e^{\sin x}$; the observed order
  is $\log_2(e_{n/2}/e_n)$, the step halved.

## Build and run

```bash
fobis build --mode examples-gnu && ./exe-examples/fd_derivative
```

or against the library:

```bash
fobis build --mode static-gnu
gfortran -I lib/gnu/static/mod docs/examples/src/fd_derivative.F90 lib/gnu/static/libgideon.a -o fd_derivative
./fd_derivative
```

## Output

<<< @/examples/output/fd_derivative.txt{text}

## What to notice

- Every column reaches its design order: 2.00, 3.99, 5.97 and 7.91 between 64 and 128 nodes. Between 16 and 32 nodes
  the orders are lower (7.23 for `c8`): $e^{\sin x}$ carries every harmonic, and 16 nodes per period do not yet resolve
  it. An order is a statement about $h \to 0$; it shows once the grid resolves the function.
- At 128 nodes the eighth-order error is $2.0 \times 10^{-10}$, eight million times below the second-order one, for a
  window of 9 nodes instead of 3. It is still a truncation error: the round-off of this kernel is about
  $\varepsilon \max|u| \sum|c_m| / h \approx 3 \times 10^{-14}$ here, and it grows as $1/h$ (exercise 3).
- The ghost count is the half-width of the widest window used, not a constant of the library. Each kernel's window is
  in the table of [Operators](/guide/operators#every-kernel).

::: warning The window's length is the caller's
`gideon_fd_d1_c8` declares `q(-4:4)`: it reads 9 values from the first element of the section it receives. gfortran
rejects a section of constant wrong length at compile time ("Actual argument contains too few elements"), but a
section with run-time bounds is not checked, even with `-fcheck=all` (gfortran 16.2): `u(i-3:i+3)` passed to
`gideon_fd_d1_c8` returns the derivative at node `i+1`. Write the section from the kernel's window, `u(i-S:i+S)` for
`q(-S:S)`, as the program does.
:::

## Exercises

1. Replace `gideon_fd_d1_c*` by `gideon_fd_d2_c*` (the same windows for orders 2 to 8) and measure against
   $u'' = (\cos^2 x - \sin x)\, e^{\sin x}$.
2. Add the tenth-order kernel, `gideon_fd_d1_c10(u(i-5:i+5), h, d)`. What else must change in the program?
3. Set `NGRID = 8` (up to 2048 nodes). The `c8` error stops near $10^{-13}$ at 512 nodes, then grows: the observed
   order becomes $-1$, round-off amplified by the division by the step.
4. The value at the midpoint $x_{i+1/2}$: `gideon_fd_face_c6(u(i-2:i+3), f)`, window `q(-2:3)`. Measure its order.
