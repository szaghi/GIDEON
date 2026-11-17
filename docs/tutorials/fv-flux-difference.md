# T2 Finite volumes, conservatively

**Goal.** From the cell averages of $u = \sin x + \tfrac12\cos 3x$ on a periodic grid: the sixth-order values at the
faces, the cell average of $u'$ as their difference, and the telescoping sum that makes a finite-volume update
conservative.

**Concepts.** [Windows](/guide/conventions#windows): an `fv` face operator reads the cell $[o-1, o]$ at offset $o$
and returns the value at the face $x = 0$, between cells 0 and 1; an `fv` cell operator reads the cell
$[o - \tfrac12, o + \tfrac12]$ and returns the average over cell 0. [Arithmetic](/guide/conventions#arithmetic), the
flux difference: an `fv` cell derivative is $(g_{+1/2} - g_{-1/2})/h^n$, $g$ the face row. Averages against point
values: [Averages and point values](/cookbook/deconvolution).

Integrate $u_t + f(u)_x = 0$ over cell $i$:

$$ \frac{d\bar u_i}{dt} = -\frac{f_{i+1/2} - f_{i-1/2}}{h}, \qquad
   \sum_i h\,\frac{d\bar u_i}{dt} = -\sum_i \left(f_{i+1/2} - f_{i-1/2}\right) = 0 \ \text{(periodic)}. $$

The cell average of $u_x$ is exactly $(u_{i+1/2} - u_{i-1/2})/h$: a face value of order $p$ gives a cell derivative
of order $p$, and every face value enters the sum twice, with opposite signs.

## The program

<<< @/examples/src/fv_flux_difference.F90{fortran}

- **Exact averages.** Cell $i$ is $[ih, (i+1)h]$, and its average is a difference of the primitive `prim`:
  $\bar u_i = (P(x_{i+1/2}) - P(x_{i-1/2}))/h$. `xf(i)` is the face $i+\tfrac12$, the right face of cell $i$.
- **Face values.** `gideon_fv_face_c6(q(i-2:i+3), f(i))`: window `q(-2:3)`, cells $i-2, \dots, i+3$, the face between
  cells $i$ and $i+1$. One value per face, `f(-1) = f(n-1)` closing the period.
- **Cell derivative.** `gideon_fv_d1_c6(q(i-3:i+3), h, d(i))`: window `q(-3:3)`, the union of the windows of the two
  faces of cell $i$. Three ghost cells on each side cover both calls.
- **Two checks.** How many cells have `d(i)` equal, bit for bit, to `(f(i) - f(i-1))/h`; and whether the sum of the
  flux differences is zero to round-off, $|\sum_i (f_{i+1/2} - f_{i-1/2})| \le n\,\varepsilon \max|f|$.

## Build and run

```bash
fobis build --mode examples-gnu && ./exe-examples/fv_flux_difference
```

or against the library:

```bash
fobis build --mode static-gnu
gfortran -I lib/gnu/static/mod docs/examples/src/fv_flux_difference.F90 lib/gnu/static/libgideon.a -o fv_flux_difference
./fv_flux_difference
```

## Output

<<< @/examples/output/fv_flux_difference.txt{text}

## What to notice

- **The cell kernel is the flux difference, bitwise**, in every cell of every grid. `gideon_fv_d1_c6` evaluates the
  two face rows exactly as `gideon_fv_face_c6` does, the same coefficients in the same order, then divides their
  difference by the step. A scheme may compute one value per face and difference them ($n$ face evaluations), or call
  the cell kernel per cell ($2n$): the numbers are the same, the cost is the caller's choice.
- **Conservation is to round-off, not exact zero.** The face value that cell $i$ uses on its right is the one cell
  $i+1$ uses on its left, so the sum telescopes; but each difference $f_{i+1/2} - f_{i-1/2}$ is rounded, and a sum of
  rounded differences is not zero. What a conservative scheme needs is the first property: one number per face, shared
  by the two cells.
- **Both errors are of order 6**, 5.98 at 128 cells. The cell-derivative error is about three times the face error:
  the leading face error is carried by the $\cos 3x$ mode (its sixth derivative is $3^6/2$ against 1 for $\sin x$), and
  differencing it across a cell differentiates it, which multiplies by 3.

## Exercises

1. Fourth order: `gideon_fv_face_c4(q(i-1:i+2), f(i))` and `gideon_fv_d1_c4(q(i-2:i+2), h, d(i))`, with `G = 2`.
2. Make it a scheme: one forward-Euler step $\bar u_i \leftarrow \bar u_i - \Delta t\,(f_{i+1/2} - f_{i-1/2})/h$ of
   $u_t + u_x = 0$, and compare $\sum_i \bar u_i$ before and after.
3. Feed point values at the cell centres instead of averages: the face error falls to order 2, the size of
   $\bar u - u(x_i) = h^2 u''/24$. `gideon_fd_average_c<order>` turns point values into averages
   ([Averages and point values](/cookbook/deconvolution)).
4. The second derivative: `gideon_fv_face_d1_c6(q(i-2:i+3), h, g(i))` is the face value of $u'$, and
   `gideon_fv_d2_c6(q(i-3:i+3), h, d2(i))` the average of $u''$ over the cell. Compare `d2(i)` with
   `(g(i) - g(i-1))/h`: equal to round-off, not bitwise, because the cell kernel divides $g_{+1/2} - g_{-1/2}$ by
   `h*h` once and the face kernel divides each $g$ by `h`. With `h` replaced by `1._R_P` in both calls the equality is
   bitwise again.
