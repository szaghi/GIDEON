program cookbook_from_adam
!< Cookbook, replace ADAM's FDV calls: the GIDEON call that replaces each ADAM family, on the data of a cubic, where
!< every operator below is exact: each line prints the value and the exact one. ADAM's half-stencil s is a run-time
!< integer here (s = 3), so the dispatchers are called; a compile-time s calls the kernel of that order. The ADAM call
!< each line replaces is in the comment above it (docs/cookbook/from-adam.md).
use gideon
implicit none
integer(I4P), parameter :: S = 3, I = 10                 ! ADAM's half-stencil; the node or cell of the calls
real(R_P),    parameter :: H = 0.1_R_P                   ! grid step
real(R_P)               :: u(0:20), a(0:20)              ! point values at x = j h; averages over [x - h/2, x + h/2]
real(R_P)               :: g(3), x, r
integer(I4P)            :: j

do j=0, 20
   x = j*H
   u(j) = p(x)
   a(j) = x**3 + x*H**2/4 - 2*x                          ! the average of p over the cell centred at x
enddo
x = I*H

! compute_derivative1_fd_centered(s=s, ds=ds, q=u(i-s:i+s), dq_ds=d): the same section, order 2s
call gideon_fd_d1_c(2*S, u(I-S:I+S), H, r)
call show('fd d1 at node i', r, 3*x**2 - 2)
! compute_derivative3_fd_centered(s=s, ds=ds, q=u(i-s:i+s), d3q_ds3=d): order 2s - 2 for d3, d4; 2s - 4 for d5, d6
call gideon_fd_d3_c(2*S - 2, u(I-S:I+S), H, r)
call show('fd d3 at node i', r, 6._R_P)
! compute_reconstruction_r_fd_centered(s=s, q=u(i+1-s:i+s), qr=f): the face i+1/2
call gideon_fd_face_c(2*S, u(I+1-S:I+S), r)
call show('fd value at face i+1/2', r, p(x + H/2))
! compute_derivative2_fv_centered(s=s, ds=ds, q=a(i-s:i+s), d2q_ds2=d): the cell average of d2
call gideon_fv_d2_c(2*S, a(I-S:I+S), H, r)
call show('fv d2 over cell i', r, 6*x)
! compute_reconstruction_r_fv_centered_d1(s=s, q=a(i+1-s:i+s), qr=f): ADAM's is unscaled, d1 times h; ds = 1 keeps it
call gideon_fv_face_d1_c(2*S, a(I+1-S:I+S), 1._R_P, r)
call show('fv d1*h at face i+1/2', r, (3*(x + H/2)**2 - 2)*H)
! compute_avg_to_point(s=s, q=a(i-s:i+s), qp=v): order 2s + 2
call gideon_fv_point_c(2*S + 2, a(I-S:I+S), r)
call show('fv point value at cell i', r, p(x))
! compute_derivative1_fv_rupwind(s=s, ds=ds, q=a(i-1:i+s), dq_ds=d): ADAM's window starts one cell early, unread
call gideon_fv_d1_ur(S, a(I:I+S), H, r)
call show('fv d1 over cell i, right side', r, (p(x + H/2) - p(x - H/2))/H)
! compute_derivative1_fv_lupwind(s=s, ds=ds, q=a(i-s:i), dq_ds=d): wrong in ADAM for s >= 2 (audit A)
call gideon_fv_d1_ul(S, a(I-S:I), H, r)
call show('fv d1 over cell i, left side', r, (p(x + H/2) - p(x - H/2))/H)
! compute_derivative2_fv_rupwind(s=s, ds=ds, q=a(i-s:i), d2q_ds2=d): ADAM's d2-d5 "rupwind" read the left side (F)
call gideon_fv_d2_ul(S - 1, a(I-S:I), H, r)
call show('fv d2 over cell i, left side', r, 6*x)
! compute_derivative2_fv_lupwind(s=s, ds=ds, q=a(i:i+s), d2q_ds2=d): and "lupwind" the right side
call gideon_fv_d2_ur(S - 1, a(I:I+S), H, r)
call show('fv d2 over cell i, right side', r, 6*x)
! compute_gradient_fd_centered(s=s, dxyz=dxyz, q=f(i-s:i+s,j-s:j+s,k-s:k+s), gradient=g): one stencil per axis
call gideon_fd_gradient_c(2*S, u(I-S:I+S), u(I-S:I+S), u(I-S:I+S), [H, H, H], g)
call show('fd gradient, x component', g(1), 3*x**2 - 2)

contains
   pure function p(x) result(v)
   !< The cubic of the data.
   real(R_P), intent(in) :: x
   real(R_P)             :: v
   v = x**3 - 2*x
   endfunction p

   subroutine show(what, r, exact)
   !< One line: the operator, its value, the exact one.
   character(*), intent(in) :: what
   real(R_P),    intent(in) :: r, exact
   print '(A,T32,2F14.8)', what, r, exact
   endsubroutine show
endprogram cookbook_from_adam
