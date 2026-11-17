# Choose an operator and its order

Four choices fix a kernel, `gideon_<data>_<op>_<bias><order>`: the data you hold, the quantity you want, the side of
the stencil, the order of accuracy. Each is settled by the scheme, not by taste; this page says how, and what each
choice costs in window width.

## The decision table

| You hold | You want | Kernel family | Window |
|---|---|---|---|
| point values at nodes | a derivative at a node | `gideon_fd_dN_c<p>` | centred, `q(-S:S)` |
| point values at nodes | the value midway between two nodes | `gideon_fd_face_c<p>`, `_bl<p>`, `_br<p>` | centred `q(1-S:S)`, upwind-biased `q(-k:k)` or `q(1-k:k+1)` |
| point values at nodes | a derivative midway between two nodes (a viscous flux) | `gideon_fd_face_dN_c<p>` | `q(1-S:S)` |
| point values inside a boundary | the value at the boundary face | `gideon_fd_face_ur<p>`, `_ul<p>` | one side, `q(1:p)` or `q(1-p:0)` |
| cell averages | the value at a face (the flux a conservative scheme differences) | `gideon_fv_face_c<p>`, `_ur<p>`, `_ul<p>` | centred `q(1-S:S)`, one side `q(1:p)` or `q(1-p:0)` |
| cell averages | a derivative at a face (a viscous flux) | `gideon_fv_face_dN_c<p>`, `_ur<p>`, `_ul<p>` | as the face value, wider by the derivative |
| cell averages | the average of a derivative over a cell | `gideon_fv_dN_c<p>`, `_ur<p>`, `_ul<p>` | `q(-S:S)`, `q(0:W)`, `q(-W:0)` |
| cell averages | the point value at the cell centre | `gideon_fv_point_c<p>` | `q(-S:S)` |
| point values | the cell average | `gideon_fd_average_c<p>` | `q(-S:S)` |
| one stencil per axis | gradient, divergence, curl, Laplacian | `gideon_{fd,fv}_<op>_c<p>` | `(-S:S)` per axis |

A cell kernel `gideon_fv_dN_*` is the flux difference of the face kernel of the derivative below,
`(g(+1/2) − g(−1/2))/ds^N`, computed inside: the same value, to the last bit, as differencing the face values yourself
([Ghost cells](/cookbook/ghosts) checks it). Call the face kernel when the flux needs work between the reconstruction
and the difference (a Riemann solver, an upwind choice per face, a flux function); call the cell kernel when it does not.

## FD or FV: the data decide

A kernel's row is exact on the data it was derived for. Point values and cell averages of the same function differ
by $h^2 u''/24 + O(h^4)$, so feeding averages to an FD kernel, or point values to an FV one, adds an $O(h^2)$ error
whatever the kernel's order. The program below takes the first derivative of $\sin 2\pi x$ with sixth-order kernels
three ways:

<<< @/examples/src/cookbook_choose.F90{fortran}

```text
<!--@include: ../examples/output/cookbook_choose.txt-->
```

The FD kernel on point values and the FV kernel on averages converge at order 6 (each against its own target: the
derivative at the centre, the average of the derivative over the cell). The FD kernel on averages stalls at order 2.
To move between the two kinds of data at high order, use the deconvolution kernels
([Averages and point values](/cookbook/deconvolution)).

## Centred or upwind

- **Centred** (`c`): the pair form, $c_m (q_m \pm q_{-m})$, makes a first derivative antisymmetric and a second one
  symmetric exactly: no numerical dissipation (the imaginary part of the modified wavenumber is zero to round-off,
  [Spectral properties](/verification/spectral)). The default for smooth fields, viscous and diffusive terms, the 3-D
  operators, and the linear part of a hybrid scheme.
- **Upwind** (`ur`, `ul`): the stencil on one side only. Name it by where the stencil lies: for a flow to the right
  (positive speed) the upwind side of a face is its **left**, so `gideon_fv_face_ul<p>`; for a flow to the left,
  `gideon_fv_face_ur<p>`. A solver switches between the two per face on the sign of the speed. The FV stencils are
  fully one-sided, and only the narrowest are stable advection operators on their own. The modified wavenumber of the left-side cell first derivative, the flux difference of the left-side faces,
  decides it (`gideon_test_spectral`, [Spectral properties](/verification/spectral)): `gideon_fv_d1_ul1` and `_ul2`
  damp every wavenumber of a right-moving wave; `gideon_fv_d1_ul3` to `_ul9` amplify a band. For `_ul3` the real part
  of its symbol is $-(\cos\xi - 1)^2 (4\cos\xi - 1)/3$, so it amplifies every wavenumber $0 < \xi < \arccos(1/4)
  \approx 0.42\pi$; the wider ones amplify more, and higher (the right-side twins mirror all this for a left-moving
  wave). The one-sided rows of 3 to 9 cells are building blocks, as the candidate stencils of a WENO reconstruction
  are: exact to their order, combined by a scheme that supplies the stability (kept by owner decision). GIDEON
  provides them; it does not make them a scheme.
- **Upwind-biased** (`bl`, `br`, point values only): $k + 1$ points on the upwind side of the face and $k$ on the
  other, order $2k + 1$, named by the side holding the extra point (`gideon_fd_face_bl<p>` for a flow to the right).
  The linear upwind interpolation of point values; on point values the fully one-sided `ur`/`ul` face values are
  extrapolations, for boundary faces ([Face values from point values](/cookbook/face-from-points)).

## The derivative and the order are two numbers

The name carries the **order of accuracy** `p`; the derivative `N` is in the operator (`d1`…`d6`). For a fixed window,
a higher derivative leaves fewer degrees of freedom for accuracy: a centred window of $2S+1$ nodes gives order $2S$ for
d1 and d2, $2S-2$ for d3 and d4, $2S-4$ for d5 and d6. So `gideon_fd_d1_c6` reads 7 nodes and `gideon_fd_d3_c6` reads 9.
The orders are measured by the generator from the row's exact degree, not assumed from the window
([The generator](/guide/generator)). The table below lists, for every operator, its orders and the points of the
window of each (`order: points`):

<!--@include: ./_tables/orders.md-->

## Which order

The error of a kernel of order $p$ at step $h$ is $C\,h^p\,|u^{(N+p)}|$, with $C$ the exact leading-error constant
the generator derives (`kerror` in the test registry). Higher $p$ pays off while $u$ is resolved and its error is above
the round-off floor, about $\varepsilon \sum_k |c_k| \max|q| / h^N$:

- **Resolution.** For a mode of wavenumber $k$, the gain from one order to the next is a factor of order $(kh)^2$:
  at four points per wavelength ($kh = \pi/2$) a wider stencil gains little. Order 6 to 10 pays off for well resolved
  turbulence or wave propagation; at under-resolved scales the dispersion error (the modified wavenumber) decides,
  [Spectral properties](/verification/spectral).
- **Round-off.** The floor grows like $h^{-N}$: a d4 or d6 at a small step reaches it early, in R4P first
  ([Single and quadruple precision](/cookbook/precision)).
- **Ghosts and cost.** The window sets the halo ([Ghost cells](/cookbook/ghosts)) and the flops: one multiply-add per
  pair for a centred kernel, per cell for an upwind one.
- **Run-time choice.** When the order is an input, one dispatcher covers every order of an operator, at a measured
  9 to 18 % over the kernel ([A run-time order, fast](/cookbook/run-time-order)).
