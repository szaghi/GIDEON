# Ghost cells for a window

**Problem.** Size the halo of an array (the ghost cells or nodes beyond the domain, filled by a boundary condition or
by a neighbour's data) so that every call of a kernel reads data that exists.

## The rule

Every kernel declares its window as an explicit-shape dummy `q(LO:HI)`, by offset from the point of the result
([Conventions](/guide/conventions)). The caller passes `q(i+LO:i+HI)`, so the reach of a loop follows from the window:

- **A node or cell operator** (`fd_dN`, `fv_dN`, `fv_point`, `fd_average`) at the points $i = 1 \dots N$ reads
  $1+\mathrm{LO} \dots N+\mathrm{HI}$: $\max(0, -\mathrm{LO})$ ghosts on the left, $\max(0, \mathrm{HI})$ on the right.
- **A face operator** (`face`, `face_dN`) at the face between $i$ and $i+1$ reads the same offsets. The $N+1$ faces of
  cells $1 \dots N$ are $i = 0 \dots N$, so the reach is $\mathrm{LO} \dots N+\mathrm{HI}$: $\max(0, 1-\mathrm{LO})$
  ghosts on the left, $\max(0, \mathrm{HI})$ on the right.
- **An upwind pair** (`_ur` and `_ul` of the same order, the side chosen per point by the flow) needs the larger of
  the two on each side: the same number on both.
- **A 3-D operator** reads one stencil per axis through the point, `(-S:S)`: $S$ ghosts per axis, and none on the
  edges or corners of a block.

Two different mistakes, two different outcomes. A halo too thin makes the last sections of a loop leave the array:
a checked build (`-fcheck=bounds`, `-C`, `-check all`) stops there with a bounds error. A section of the wrong
**length** is not caught: the explicit-shape dummy takes the section's elements in order, so `u(i-3:i+3)` passed to a
kernel whose window is `q(-4:4)` is read as if it started at offset −4, and the result is that of node `i+1`, silently,
even with `-fcheck=all` (gfortran 16.2); only a constant-length section of the wrong size fails, at compile time. The
window bounds of each kernel are the caller's responsibility: write the section from the kernel's own `LO` and `HI`, as
the code below does with named constants.

## The table

Ghosts per side for every operator and order, computed by `docs/scripts/make_cookbook_tables.py` from the `q(LO:HI)`
each compiled kernel declares (`src/lib/*.F90`) with the rule above; `left / right` where the sides differ:

<!--@include: ./_tables/ghosts.md-->

Two consequences. A centred node or cell kernel of order $p$ needs $p/2$ ghosts for d1 and d2, one more for d3 and d4,
two more for d5 and d6. A one-sided kernel needs none on the side away from its stencil, but the pair a solver uses
needs the full width on both.

## Code

A node loop with `gideon_fd_d1_c6` (window `q(-3:3)`), and a face loop of the upwind pair `gideon_fv_face_ul3`
(`q(-2:0)`) and `gideon_fv_face_ur3` (`q(1:3)`) over the faces $0 \dots N$, on a periodic grid; the halo is declared
from the window bounds, and the ghosts are copies from the other end of the domain:

<<< @/examples/src/cookbook_ghosts.F90{fortran}

```text
<!--@include: ../examples/output/cookbook_ghosts.txt-->
```

Both loops need 3 ghosts on each side, as the table says (`fd_d1_c` at order 6, `fv_face_ur + _ul` at order 3). The
flux difference of the left-side faces, `(f(i+1/2) − f(i−1/2))/h`, equals the cell kernel `gideon_fv_d1_ul3` in every
cell, bitwise: the cell kernel computes the two face sums the same way and in the same order.

**Reference.** The windows: [Conventions](/guide/conventions), every kernel's in [Operators](/guide/operators);
a wrong window (ADAM's audit item B): [Replace ADAM's FDV calls](/cookbook/from-adam).
