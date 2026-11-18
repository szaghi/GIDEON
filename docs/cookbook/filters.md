# Filtering a central scheme

**Problem.** A central scheme leaves the odd-even mode $(-1)^i$, the shortest wave of the grid, undamped: it has no
derivative to see it with, so round-off and nonlinear products feed it and nothing removes it. A filter applied after
each step removes it, and should leave the resolved waves alone, to the order of the scheme.

## The kernels

| Kernel | What | Window |
|---|---|---|
| `gideon_fd_filter_c<2k>`, `gideon_fv_filter_c<2k>` | explicit centred filter of order $2k$: exact below degree $2k$, the odd-even mode removed (the binomial filter $u - (-1)^k \delta^{2k} u/4^k$) | `q(-k:k)` |
| `gideon_fd_filter_p<2k>` | the right side of the compact filter of order $2k$, $\alpha_f$ an argument; the band $[\alpha_f, 1, \alpha_f]$ solved along the line | `q(-k:k)` |

The explicit filter is one kernel call per point. The compact one is a right side per point and a solve:

```fortran
do i=1, n
   call gideon_fd_filter_p8(u(i-4:i+4), alpha, uf(i))
enddo
call gideon_tri_periodic_const(n, alpha, 1._R_P, alpha, uf, w1, w2)   ! |alpha| < 1/2
```

$\alpha_f$ sets the sharpness: larger values leave more of the resolved band untouched (Gaitonde and Visbal use
0.3–0.49); at $1/2$ the band is singular at the odd-even mode.

## Code

Point values of $\sin x + \cos(2x)/2$ plus an odd-even mode of amplitude 0.1 on a periodic line, every filter of
orders 2–10, explicit and compact ($\alpha_f = 0.45$); the error on the smooth wave on 16 and 32 points, and what is
left of the odd-even mode:

<<< @/examples/src/cookbook_filter.F90{fortran}

```text
<!--@include: ../examples/output/cookbook_filter.txt-->
```

- Every filter removes the odd-even mode to round-off: that is its definition (the right symbol is zero at $\xi = \pi$).
- Every filter keeps the smooth wave to its order; the compact one, at equal order, is 16 to 18 times more accurate on
  this wave, at every order alike: its transfer function stays near 1 over more of the band.
- Filter after every step of the central scheme, on each line of each direction. The compact filter costs a solve per
  line; the explicit one is a single pass.

**Reference.** Roadmap issue #2, phase 2 (2.1, 2.2). The rows: [Operators](/guide/operators#explicit-filters),
[Compact operators](/guide/compact); the transfer function in $[0, 1]$ and the odd-even mode removed are checked for
every compact filter by [A4](/verification/lines).
