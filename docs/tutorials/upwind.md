# T3 Upwind faces and their mirror

**Goal.** The third-order face value from the cells on either side of a face; the mirror identity that makes every
left-side kernel the reflection of its right-side twin, checked bit for bit; the upwind face values of the advection
$u_t + a u_x = 0$ for both signs of $a$.

**Concepts.** [Names](/guide/conventions#names): the bias names **where the stencil lies**, `ur` on the right of the
face, `ul` on the left. [Windows](/guide/conventions#windows): an upwind face from $S$ cells reads `q(1:S)` (right)
or `q(1-S:0)` (left); an upwind cell operator reads `q(0:S)` or `q(-S:0)`. Mirrors:
[Arithmetic](/guide/conventions#arithmetic). T2 for the face values as fluxes.

The flow carries information from upwind: for $a > 0$ the face $i+\tfrac12$ takes its value from the cells on its
left, `gideon_fv_face_ul3(u(i-2:i), f)`; for $a < 0$ from the cells on its right, `gideon_fv_face_ur3(u(i+1:i+3), f)`.
The order-3 row is $\tfrac{11}{6}, -\tfrac76, \tfrac13$, the nearest cell first: the stencil lies wholly on one side
of the face ($i-2, i-1, i$ for `ul3`, the leftmost candidate stencil of WENO5), not the upwind-biased $i-1, i, i+1$.

## The mirror

Reflect the data about the face, $r_m = q_{1-m}$: the left-side kernel on $q$ is the right-side kernel on $r$. For a
cell operator the reflection is about the cell centre, $r_m = q_{-m}$, and the $n$-th derivative changes sign as
$(-1)^n$. GIDEON stores no left-side row: each `ul` kernel indexes the `ur` row backwards, so the identity holds by
construction, and it holds bitwise because the reflected kernel performs the same operations on the same numbers.

## The program

<<< @/examples/src/upwind.F90{fortran}

- **The mirror check.** `q(i:i-2:-1)` is the reflected window, passed as it is: a section with a negative stride is
  copied into a contiguous temporary for the explicit-shape dummy. For `d1` the comparison is with `-dr`, for `d2`
  with `+dr`. The data are arbitrary: the identity is algebraic, not a property of smooth fields.
- **The residual.** $R_i = -a\,(f_{i+1/2} - f_{i-1/2})/h$ with $f$ the upwind face value, against the exact
  $-a\,(u(x_{i+1/2}) - u(x_{i-1/2}))/h$. Cell averages of $e^{\sin x}$ come from an 8-point Gauss-Legendre rule per
  cell, exact far below the errors measured.
- **Ghosts.** For $a > 0$ the faces reach 2 cells to the left of cell 0; for $a < 0$ 3 cells to the right of cell
  $n-1$. `G = 3` covers both.

## Build and run

```bash
fobis build --mode examples-gnu && ./exe-examples/upwind
```

or against the library:

```bash
fobis build --mode static-gnu
gfortran -I lib/gnu/static/mod docs/examples/src/upwind.F90 lib/gnu/static/libgideon.a -o upwind
./upwind
```

## Output

<<< @/examples/output/upwind.txt{text}

## What to notice

- **The mirror is exact**: 1000 of 1000 windows for the face value, the first and the second derivative, with the
  sign $(-1)^n$. This is what [C3 Mirrors and symmetries](/verification/properties) verifies for every left-side
  kernel.
- **Order 3 from three cells**, 2.98 at 128 cells, for either sign of $a$.
- **The two columns are equal.** $e^{\sin x}$ is symmetric about $x = \pi/2$, and the reflection $x \to \pi - x$ maps
  faces onto faces ($\pi$ is a multiple of $h$ for even $n$). Reflected, the flow to the right becomes the flow to the
  left, and the mirror identity carries the errors across: the same numbers to the digits printed.
- **Accuracy says nothing about stability.** The downwind face value (`ur3` for $a > 0$) has an error of the same
  order, and the upwind one is not stable by being upwind. The flux difference of `ul3` is
  $(11\bar u_i - 18\bar u_{i-1} + 9\bar u_{i-2} - 2\bar u_{i-3})/(6h)$; a mode $e^{\mathrm{i}\theta x/h}$ of
  $u_t + a u_x = 0$ grows at the rate $\tfrac{a}{h}(1-\cos\theta)^2(4\cos\theta - 1)/3$, positive for
  $\theta < \arccos\tfrac14 \approx 1.318$. Of the left-side first derivatives only those of 1 and 2 cells damp every
  wavenumber; the wider one-sided rows amplify a band, and they are building blocks (the candidate stencils a WENO
  reconstruction weighs), not schemes on their own ([A3 Spectral properties](/verification/spectral)).

## Exercises

1. Fifth order: `gideon_fv_face_ul5(u(i-4:i), f(i))` for $a > 0$, `gideon_fv_face_ur5(u(i+1:i+5), f(i))` for
   $a < 0$, and `G = 5`.
2. The cell kernel is the flux difference (T2), upwind too: compare `gideon_fv_d1_ul3(u(i-3:i), h, d)` with
   `(f(i) - f(i-1))/h` for $a > 0$, bit for bit.
3. The symbol, from the compiled kernels: with `h = 1`, call `gideon_fv_d1_ul3` on the window
   $q_o = \cos(o\theta)$, $o = -3, \dots, 0$; the result is $-1$ times the growth rate of the mode $\theta$ (for
   $a = 1$). Scan $\theta \in (0, \pi]$, find the band where it is negative and compare its edge with
   $\arccos\tfrac14$. Repeat for `gideon_fv_d1_ul2`: no such band.
4. The third derivative: `gideon_fv_d3_ul3` against `gideon_fv_d3_ur3` on the reflected window. Which sign?
