# API

`use gideon` gives the kinds (`I4P`, `R4P`, `R8P`, `R_P` from PENF), `GIDEON_NAN`, the rows (`gideon_tables`) and every
operator. The full list, with windows and orders, is on [Operators](./operators).

## 1-D kernels

```fortran
pure subroutine gideon_<data>_<op>_<bias><order>(q, [ds,] r)
real(R_P), intent(in)  :: q(LO:HI) ! the window, by offset (Conventions)
real(R_P), intent(in)  :: ds       ! the grid step (derivatives only)
real(R_P), intent(out) :: r        ! the result
```

```fortran
use gideon
real(R_P) :: u(-2:102), du, uf
integer   :: i
! sixth-order first derivative at node i, point values
call gideon_fd_d1_c6(u(i-3:i+3), 0.01_R_P, du)
! fifth-order value at the face i+1/2 from the five cells on its left (cell averages)
call gideon_fv_face_ul5(u(i-4:i), uf)
```

## Run-time order

```fortran
pure subroutine gideon_<data>_<op>_<bias>(order, q, [ds,] r)
integer(I4P), intent(in) :: order
real(R_P),    intent(in) :: q(lo(order):hi(order)) ! the window of that order
```

```fortran
call gideon_fd_d1_c(order, u(i-order/2:i+order/2), ds, du) ! GIDEON_NAN if no kernel has that order
```

The dispatch costs 9–22 % over the kernel (see [Performance](./performance)): with an order known at compile time,
call the kernel.

## 3-D operators

```fortran
call gideon_fd_gradient_c6(qx, qy, qz, dxyz, g)                       ! qx(-3:3): the field along x, ...
call gideon_fd_divergence_c6(vx_x, vy_y, vz_z, dxyz, div)             ! each component along its axis
call gideon_fd_curl_c6(vy_x, vz_x, vx_y, vz_y, vx_z, vy_z, dxyz, c)   ! the off-axis components
call gideon_fd_laplacian_c6(qx, qy, qz, dxyz, lap)
```

`gideon_fv_*` take cell averages along each axis; `gideon_<data>_<op>_c(order, ...)` dispatch at run time.

## From ADAM

| ADAM | GIDEON |
|---|---|
| `compute_derivativeN_fd_centered(s, ds, q(1-s:1+s), d)` | `gideon_fd_dN_c<order>(q(-S:S), ds, d)`, the node at offset 0 |
| `compute_reconstruction_r_fd_centered(s, q(1-s:s), f)` | `gideon_fd_face_c<2s>(q(1-s:s), f)` |
| `compute_derivativeN_fv_centered(s, ds, q, d)` | `gideon_fv_dN_c<order>` |
| `compute_reconstruction_r_fv_centered[_dN](s, q, f)` | `gideon_fv_face_c<2s>`, `gideon_fv_face_dN_c<order>(q, 1, f)` (ADAM's is unscaled) |
| `compute_avg_to_point`, `compute_point_to_avg` | `gideon_fv_point_c<2s+2>`, `gideon_fd_average_c<2s+2>` |
| `compute_derivative1_fv_rupwind`, `compute_derivative{2..5}_fv_lupwind` | `gideon_fv_dN_ur<order>` (right side) |
| `compute_derivative1_fv_lupwind`, `compute_derivative{2..5}_fv_rupwind` | `gideon_fv_dN_ul<order>` (left side) |
| `compute_{gradient,divergence,curl,laplacian}_{fd,fv}_centered[_dev]` | `gideon_{fd,fv}_{gradient,divergence,curl,laplacian}_c<2s>` |

The orders and the bitwise or tolerance-level equivalence of each pair are in `notes/adam-fdv-audit.md`, §5.
