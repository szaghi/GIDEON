# Averages and point values

**Problem.** A finite-volume solver holds cell averages $\bar u_i$, but some of its work needs point values: a
nonlinear flux, a source term, an initial condition given as a function, a diagnostic compared with point data. The
two differ by $h^2 u''/24 + O(h^4)$; ignoring the difference makes the step second order whatever the order of the
reconstruction around it.

## The kernels

| Kernel | From | To | Window | Order |
|---|---|---|---|---|
| `gideon_fv_point_c<p>` | cell averages | the point value at the centre of cell 0 | `q(-S:S)` | $p = 2S+2$: 4, 6, 8, 10, 12 |
| `gideon_fd_average_c<p>` | point values at the cell centres | the average over cell 0 | `q(-S:S)` | the same |

Both are symmetric rows in pair form, derived by their meaning: `gideon_fv_point_c4` is exact on every cubic,
$u_i = \tfrac{13}{12}\bar u_i - \tfrac{1}{24}(\bar u_{i+1} + \bar u_{i-1})$, i.e. $\bar u_i - \delta^2 \bar u_i/24$;
`gideon_fd_average_c4` is $\tfrac{11}{12}u_i + \tfrac{1}{24}(u_{i+1} + u_{i-1})$. They are each other's inverse to
their order, not exactly: a round trip returns the data plus an $O(h^p)$ error. ADAM's `compute_avg_to_point` and
`compute_point_to_avg` are the same rows (bitwise, [Replace ADAM's FDV calls](/cookbook/from-adam)).

## When an FV solver needs them

- **A nonlinear flux.** The average of $f(u)$ over a cell is not $f(\bar u)$: for $f = u^2/2$ they differ by
  $h^2 (u')^2/24 + O(h^4)$. The route to order 4 and above: deconvolve the averages to point values, evaluate $f$ at
  the points, convolve back to the average of $f$. Each step costs one
  window, so the ghost layer is two half-widths (here 2 for order 4, 4 for order 6).
- **Initial conditions.** A solution given as a function is sampled at the centres and averaged with
  `gideon_fd_average_c<p>`, at the order of the scheme; where the cell integral is known in closed form, that is exact.
- **Diagnostics and output.** A plot, a probe, an error norm against a point-wise reference: `gideon_fv_point_c<p>`
  first, or the comparison carries the $O(h^2)$ difference.

In more dimensions the average over a box is the product of the three 1-D averages, so the 1-D kernels applied along
each axis in turn give the box conversion at the same order.

## Code: the average of $u^2/2$

The cell averages of $u = \sin 2\pi x$ on a periodic grid, and the average of $f = u^2/2$ over each cell three ways:
$f(\bar u)$; deconvolution, $f$ at the points and convolution at order 4; the same at order 6:

<<< @/examples/src/cookbook_deconvolution.F90{fortran}

```text
<!--@include: ../examples/output/cookbook_deconvolution.txt-->
```

$f(\bar u)$ converges at order 2, whatever the scheme around it; the deconvolve-evaluate-convolve route at the order of
its kernels. The point-value array `fp4` is computed on cells $-1 \dots N+2$, one half-width beyond the cells whose
average is wanted, because the convolution needs it there.

**Reference.** The expansion $\bar u = u + h^2 u''/24 + \dots$, as ADAM cites it from P. McCorquodale, P. Colella,
*A high-order finite-volume method for conservation laws on locally refined grids*, Comm. App. Math. Comp. Sci. 6
(2011) 1–25 ([Sources](/guide/sources)); GIDEON derives the rows from the cell averages of monomials, not from the
expansion. The rows' exactness and order:
[Exactness and leading error](/verification/exactness), [Orders of accuracy](/verification/orders).
