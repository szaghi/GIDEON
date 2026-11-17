!< Audit check (2026-11-17): convergence of ADAM's FV upwind first derivatives on cell averages of sin(2 pi x).
!< At fixed x0 = 0.3, h = 1/50 .. 1/400, the error of each kernel against the exact cell average of the derivative on
!< cells j-1, j, j+1. Windows as declared by the kernels: right-upwind q(0:) [the caller passes qa(j:j+s), as ADAM's
!< test_fdv_operators_trigonometric does], left-upwind q(-s:0) = qa(j-s:j). Not part of the GIDEON build: compile with
!< PENF and adam_fdv_operators_library.F90 (this directory), e.g. scripts in notes/adam-fdv-audit.md.
program check_upwind_d1
use penf
use adam_fdv_operators_library
implicit none
real(R8P)    :: qa(-40:1000), h, d, pi, e(-1:1), x0
integer      :: s, j, m, k

pi = acos(-1._R8P)
x0 = 0.3_R8P
do s=1, 5
   print '(A,I0)', 's = ', s
   do k=0, 3
      h = 1._R8P/(50*2**k)
      j = nint(x0/h)
      do m=-40, 1000
         qa(m) = -(cos(2*pi*(m+0.5_R8P)*h) - cos(2*pi*(m-0.5_R8P)*h))/(2*pi*h) ! averages of sin(2 pi x)
      enddo
      do m=-1, 1
         e(m) = (sin(2*pi*(j+m+0.5_R8P)*h) - sin(2*pi*(j+m-0.5_R8P)*h))/h      ! averages of its derivative
      enddo
      call compute_derivative1_fv_rupwind(s=s, ds=h, q=qa(j:j+s), dq_ds=d)
      write(*,'(A,ES9.2,A,3ES10.2)', advance='no') '  h=', h, '  RU err vs j-1,j,j+1:', abs(d - e)
      call compute_derivative1_fv_lupwind(s=s, ds=h, q=qa(j-s:j), dq_ds=d)
      write(*,'(A,3ES10.2)') '  | LU:', abs(d - e)
   enddo
enddo
endprogram check_upwind_d1
