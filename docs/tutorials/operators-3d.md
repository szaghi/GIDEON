# T5 Gradient, divergence, curl

**Goal.** The gradient, Laplacian, divergence and curl at every node of a periodic $n^3$ grid of point values: gather
one contiguous stencil per direction (and per component) out of the 3-D arrays, call the operators, compare with the
exact values, and observe orders 2, 4 and 6.

**Concepts.** [3-D operators](/guide/api#_3-d-operators) and their place among the [families](/guide/operators#families):
ADAM's FNL signature, one stencil per direction and component, gathered by the caller. T4 for the run-time order;
[A2 3-D operators](/verification/operators-3d) for their verification.

The fields:

$$ \varphi = \sin x \sin y \sin z, \qquad \nabla\varphi = (\cos x \sin y \sin z,\ \sin x \cos y \sin z,\
   \sin x \sin y \cos z), \qquad \nabla^2\varphi = -3\varphi, $$

$$ \mathbf v = (\sin x \cos y,\ \sin y \cos z,\ \sin z \cos x), \qquad
   \nabla\cdot\mathbf v = \cos x\cos y + \cos y\cos z + \cos z\cos x, \qquad
   \nabla\times\mathbf v = (\sin y \sin z,\ \sin z \sin x,\ \sin x \sin y). $$

## The gather

| Operator | Stencils | From the 3-D arrays at node `(i,j,k)` |
|---|---|---|
| gradient, Laplacian | `qx`, `qy`, `qz`: the field along each axis | `phi(i-s:i+s,j,k)`, `phi(i,j-s:j+s,k)`, `phi(i,j,k-s:k+s)` |
| divergence | `vx_x`, `vy_y`, `vz_z`: each component along its own axis | `vx(i-s:i+s,j,k)`, `vy(i,j-s:j+s,k)`, `vz(i,j,k-s:k+s)` |
| curl | `vy_x`, `vz_x`, `vx_y`, `vz_y`, `vx_z`, `vy_z`: the off-axis components | `vy(i-s:i+s,j,k)`, `vz(i-s:i+s,j,k)`, ... |

`vy_x` is the y component along x. A stencil along x is contiguous in memory; along y and z it is strided. The
operators take explicit-shape dummies (`qx(-S:S)`), so a strided section passed directly would be copied into a
temporary by the compiler at every call, and in device code such a copy is what ADAM's FNL signature exists to avoid.
The caller copies into small local arrays instead: the data movement is explicit, and one gather of $\varphi$ serves
both the gradient and the Laplacian.

## The program

<<< @/examples/src/operators_3d.F90{fortran}

- **Ghosts.** The fields are evaluated on the ghost nodes too, as a halo exchange would deliver them; `G = 3` is the
  half-width of the order-6 stencils.
- **Run-time order.** The local arrays are declared `(-G:G)` once; for order `order` the program fills and passes
  the sections `qx(-s:s)`, `s = order/2`, to the dispatchers `gideon_fd_<op>_c(order, ...)`: the windows of T4.
- **Errors.** The maximum over the nodes, and over the components for the vector results.

## Build and run

```bash
fobis build --mode examples-gnu && ./exe-examples/operators_3d
```

or against the library:

```bash
fobis build --mode static-gnu
gfortran -I lib/gnu/static/mod docs/examples/src/operators_3d.F90 lib/gnu/static/libgideon.a -o operators_3d
./operators_3d
```

## Output

<<< @/examples/output/operators_3d.txt{text}

## What to notice

- **Every operator reaches its order**: 2.00, 4.00 and 5.99 between 32 and 64 nodes per direction, the orders of the
  1-D kernels they sum (each 3-D operator is the 1-D kernels in ADAM's order, x then y then z).
- **The divergence error is three times the gradient's**: each of its three terms carries the error of one first
  derivative of a sine, and the three peak together where $\cos x = \cos y = \cos z = 1$. The gradient and the curl
  have equal errors: each of their components is one first derivative of a sine times a product of two others.
- **The Laplacian at order 6 is more accurate than the gradient** (4.8E-09 against 6.4E-09 at $n = 64$), at order 2
  less (2.4E-03 against 1.6E-03): the leading constants of the d2 and d1 rows differ, the orders do not.

## Exercises

1. Finite volumes: replace the fields by their cell averages (a product of 1-D averages; the average of $\sin$ over
   $[x - h/2, x + h/2]$ is $\sin x \cdot \sin(h/2)/(h/2)$) and call `gideon_fv_gradient_c`, ...; the exact answer is
   the cell average of the derivative.
2. Order known at compile time: replace the dispatchers by `gideon_fd_gradient_c6` and friends, and the sections by
   the whole local arrays declared `(-3:3)`.
3. Discrete identities: compute $\nabla\varphi$ at every node (ghosts included, by periodicity), then its curl. The
   result is round-off, not $O(h^p)$: centred first derivatives along different axes commute on a uniform periodic
   grid. Does $\nabla\cdot(\nabla\times\mathbf v)$ vanish the same way?
