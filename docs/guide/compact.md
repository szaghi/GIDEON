# Compact operators and line solvers

A compact (implicit, Padé) operator couples the unknowns along a line:

$$ \sum_{j=-l}^{l} \alpha_j\, u^{(n)}_{i+j} = \frac{1}{h^n} \sum_{m} w_m\, q_{i+m} . $$

The band $\alpha$ on the left (tridiagonal, $l = 1$, or pentadiagonal, $l = 2$, with $\alpha_0 = 1$) buys accuracy and
resolution that an explicit stencil of the same width cannot have: `gideon_fd_d1_p6` reaches order 6 with a right side
of five points where the explicit `gideon_fd_d1_c6` needs seven, and its error is about 1/15 of the explicit one's
on the same grid ([T7](/tutorials/compact-derivative)). The price is a solve along every line.

GIDEON splits a compact operator into three things (roadmap #2, decision D2), so that every rule of the explicit
kernels still holds:

| Part | What | Where |
|---|---|---|
| the right side | an ordinary per-point kernel: window by offset, `pure`, device-callable, its row a compile-time constant | `gideon_compact`, the row `<NAME>_RIGHT` |
| the band | named constants, pair form (index $j$ multiplies the unknowns at $\pm j$) or by offset when not symmetric | `gideon_tables`, `<NAME>_LEFT` |
| the solve | one line per call, the scratch the caller's arrays | `gideon_line_solvers` |

The caller computes the right side point by point with the kernel, then solves the band on the line:

```fortran
do i=0, n-1                                   ! the right side, one kernel call per node
   call gideon_fd_d1_p6(u(i-2:i+2), h, du(i))
enddo
call gideon_tri_periodic_const(n, FD_D1_P6_LEFT(1), FD_D1_P6_LEFT(0), FD_D1_P6_LEFT(1), du, w1, w2)
```

## The operators

<!--@include: ./_compact-table.md-->

- **Node derivatives** (`gideon_fd_d1_p<o>`, `gideon_fd_d2_p<o>`): Lele's centred first and second derivatives,
  tridiagonal orders 4 and 6, pentadiagonal 8 and 10. The band and the right side are solved together for the highest
  order the window allows: the derived rows are Lele's (1992) as the generator recalls them (a cross-check, not
  evidence; [Sources](/guide/sources)).
- **Compact filters** (`gideon_fd_filter_p<2N>`): Gaitonde and Visbal's form,
  $\alpha_f \tilde u_{i-1} + \tilde u_i + \alpha_f \tilde u_{i+1} = \sum_m w_m(\alpha_f)\, u_{i+m}$, exact on the
  polynomials of degree below $2N$ and zero on the odd-even mode. The conditions are linear in $\alpha_f$, so the
  weights are `_RIGHT + alpha_f * _ALPHA`, the kernel takes `alpha` as an argument and the caller passes the band
  `[alpha, 1, alpha]` to the solver; $|\alpha_f| < 1/2$ (at $1/2$ the band is singular at the odd-even mode).
- **Face values and derivatives from point values** (`gideon_fd_face_p<o>`, `gideon_fd_face_d1_p<o>`): the unknowns at
  the faces, the data at the nodes on either side (staggered use: the viscous flux of a compact scheme).
- **Face values from cell averages**: centred (`gideon_fv_face_p<o>`) and upwind-biased (`gideon_fv_face_pbr<o>`,
  `gideon_fv_face_pbl<o>`, $o = 5, 7$), a tridiagonal band whose two off-diagonals are both solved for, $2k+1$ cells
  with one more on the upwind side. `gideon_fv_face_pbl5`, band $(1/2, 1, 1/6)$ on cells $(1/18, 19/18, 5/9)$, is
  Pirozzoli's (2002) compact upwind face value: the linear branch of his compact–WENO hybrid. `pbl` reads `pbr`'s row
  mirrored, bitwise, and keeps its own (mirrored) band.

No run-time order dispatcher exists for the compact operators: the band changes with the order, so a dispatcher of the
right side alone would be half an operator.

## The solvers

| Routine | Band | Line |
|---|---|---|
| `gideon_tri_bounded_{const,var,masked}` | $a_i x_{i-1} + b_i x_i + c_i x_{i+1}$ | bounded: Thomas' algorithm |
| `gideon_tri_periodic_{const,var,masked}` | the same, indices modulo $n$ ($n \ge 3$) | periodic |
| `gideon_penta_bounded_{const,var,masked}` | $e_i x_{i-2} + a_i x_{i-1} + b_i x_i + c_i x_{i+1} + f_i x_{i+2}$ | bounded: elimination without pivoting |
| `gideon_penta_periodic_{const,var,masked}` | the same, indices modulo $n$ ($n \ge 5$) | periodic |

Every routine, the private helpers included (they are generated in the same three variants):

<!--@include: ./_solvers-table.md-->

- `_const`: the band a constant (scalars): the compact operators. `_var`: one coefficient per row (arrays). `_masked`:
  the constant band, with the rows where `mask(i)` replaced by the identity: their right side is the value the caller
  gives there (below). The three are emitted from one template by the generator, so they cannot drift: on the same
  band they agree bitwise ([C8](/verification/solvers)).
- A periodic line is solved as its interior (a bounded band, for the right side and for the columns of the border
  unknowns, with one factor) followed by the border's 1×1 or 2×2 Schur complement, solved by elimination: no
  Sherman–Morrison shift to choose, and a border row replaced by the identity returns its value exactly.
- Elimination without pivoting is stable for diagonally dominant bands and for symmetric positive definite ones, and
  every GIDEON band is one or the other. `gideon_fd_d1_p10`'s band $[1/20, 1/2, 1, 1/2, 1/20]$ is not diagonally
  dominant ($2\alpha + 2\beta = 1.1$); its symbol $1 + \cos\xi + \cos(2\xi)/10 \ge 1/10$ makes it positive definite.
- A coefficient that reaches outside a bounded line ($a_1$, $c_n$, …) is never read.
- Device use: one line per thread, `!$acc parallel loop` over the lines, the scratch a column of the caller's arrays
  (`w(:,l)`), the inner loops `!$acc loop seq`. Nothing is allocated.

## Compact-WENO coupling

In a hybrid compact–WENO/TENO scheme the face values come from the compact operator where the solution is smooth and
from a nonlinear scheme where a sensor flags a discontinuity (Pirozzoli 2002). For the solve, a flagged face drops its
compact row: its row becomes the identity and its right side the nonlinear value. That is `_masked`:

```fortran
do i=1, n
   mask(i) = sensor(i) > threshold                 ! the caller's switch
   if (mask(i)) then
      f(i) = weno_value(i)                        ! the caller's nonlinear value
   else
      call gideon_fv_face_pbl5(u(i-1:i+1), f(i))  ! the compact row's right side
   endif
enddo
call gideon_tri_periodic_masked(n, FV_FACE_PBL5_LEFT(-1), FV_FACE_PBL5_LEFT(0), FV_FACE_PBL5_LEFT(1), mask, f, w, z)
```

GIDEON provides the operator and the solve; the sensor, the threshold and the nonlinear value are the caller's (SETH's:
decision D7). The flagged faces return their given values exactly, every face keeps one value (the flux form stays
conservative), and away from the flagged region the compact order is kept: the band's inverse decays geometrically,
so the influence of a step falls below the truncation error a few dozen cells away ([A5](/verification/coupling)).
A worked example: [Compact-WENO rows](/cookbook/compact-weno).

## Not here (yet)

- Boundary closures for bounded lines (D4: periodic lines only in this roadmap; the bounded solvers exist).
- Optimised (dispersion-relation-preserving, selective) filters and rows: their coefficients are not derivable by
  meaning (D3).
- A multi-dimensional compact operator: apply along one axis at a time, the caller's loop over lines.
