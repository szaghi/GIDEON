# Face values and derivatives from point values

**Problem.** A finite-difference solver holds point values $u_i$ at the nodes and needs quantities at the faces
$i + 1/2$ between them: the gradient of a viscous flux, the face value an upwind flux reads, or the value at a
boundary face that no node sits on. Written by hand, each is a row of fractions copied from a paper (FLUME writes its
viscous row $\big(27(W_1 - W_0) - (W_2 - W_{-1})\big)/24\,\Delta x$ inline); GIDEON derives them by their meaning, like
every other row.

## The kernels

All read point values at $x = o - 1/2$ around the face $x = 0$ (offset 0 is node $i$, offset 1 node $i+1$), so the
section for the face $i + 1/2$ is `u(i+LO:i+HI)`:

| Kernel | Returns | Window | Orders |
|---|---|---|---|
| `gideon_fd_face_c<p>` | the value at the face, centred | `q(1-S:S)` | 2–10 |
| `gideon_fd_face_dN_c<p>` | the N-th derivative at the face, centred (N = 1 … 5) | `q(1-S:S)` | d1 2–10; d2, d3 2–8; d4, d5 2–6 |
| `gideon_fd_face_bl<p>` | the value at the face, upwind-biased: $k+1$ points on its left, $k$ on its right | `q(-k:k)` | 3, 5, 7, 9 |
| `gideon_fd_face_br<p>` | the same, one point more on the right | `q(1-k:k+1)` | 3, 5, 7, 9 |
| `gideon_fd_face_ur<p>` | the value at the face from the $S$ points on its right only (extrapolation) | `q(1:S)` | 1–9 |
| `gideon_fd_face_ul<p>` | the same from the points on its left | `q(1-S:0)` | 1–9 |

`gideon_fd_face_d1_c4` is FLUME's row, $\tfrac{9}{8}(q_1 - q_0) - \tfrac{1}{24}(q_2 - q_{-1})$, divided by the step;
`gideon_fd_face_bl3` is $(-q_{-1} + 6 q_0 + 3 q_1)/8$ and `gideon_fd_face_bl5`
$(3 q_{-2} - 20 q_{-1} + 90 q_0 + 60 q_1 - 5 q_2)/128$, the linear upwind interpolation of point values that a WENO
scheme on point values starts from. The left-side and left-biased kernels are the mirrors of the right ones
($o \to 1 - o$), proven by the generator and tested bitwise ([Mirrors and symmetries](/verification/properties)).
Each has a run-time order dispatcher: `gideon_fd_face_d1_c(order, ...)`, `gideon_fd_face_bl(order, ...)`, …
([API](/api/#run-time-order-dispatchers)).

## Which one

- **A viscous or diffusive flux** $\mu\, \partial u/\partial x$ at the face: `gideon_fd_face_d1_c<p>`, the same order as
  the convective part. It reads the nodes on both sides, so it is centred by construction; the second and higher face
  derivatives serve hyperviscosity and the fourth-order dissipation terms.
- **A convective face value**: centred (`gideon_fd_face_c<p>`) where no upwinding is wanted, as in a split or
  energy-preserving form; upwind-biased (`bl` when the flow goes to the right, `br` to the left) where the scheme needs
  the dissipation of the upwind side. A biased stencil of $2k+1$ points is order $2k+1$, one less than the centred
  stencil of $2k+2$ points, and its error is larger at equal width (below): the price of the bias, which damps the
  shortest waves instead of letting them disperse. That trade is the same one the FV upwind kernels make
  ([Spectral properties](/verification/spectral)).
- **A boundary face** where only interior nodes exist: `gideon_fd_face_ur<p>` at a left boundary (the nodes lie to its
  right), `ul` at a right one. One-sided extrapolation amplifies the data's noise as its width grows: the sum of the
  absolute weights is 2 for `ur2`, 7/2 for `ur3`, 6 for `ur4` ($(35, -35, 21, -5)/16$) and about 110 for `ur9`,
  against 5/4 for the centred `c4`. A closure takes the lowest order that keeps the scheme's global order, as for any
  one-sided stencil.

## Code

The three uses, each measured by halving the step: the viscous derivative at every face of a periodic grid; centred and
upwind-biased face values; the value at the boundary face $x = 0$ from the nodes $x = h/2, 3h/2, \dots$ inside:

<<< @/examples/src/cookbook_face_from_points.F90{fortran}

```text
<!--@include: ../examples/output/cookbook_face_from_points.txt-->
```

Every kernel shows its order. The five-point `bl5` is less accurate than the four-point centred `c4` on this smooth
field, and `bl3` than both: the upwind bias costs accuracy where nothing needs damping, and buys stability where
something does. The ghost cells each kernel needs follow from its window ([Ghost cells for a window](/cookbook/ghosts)).

**Reference.** Roadmap issue #1, phase 2 (2.1, 2.2). The rows are derived in exact rationals by the generator
([The generator](/guide/generator)); their orders are observed in quadruple precision with every other kernel
([Orders of accuracy](/verification/orders)).
