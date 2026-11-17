program cookbook_choose
!< Cookbook, choose an operator: the first derivative of u = sin(2 pi x) on a periodic grid of n cells, sixth-order
!< centred kernels. The FD kernel on point values and the FV kernel on cell averages reach order 6; the FD kernel fed
!< with cell averages (the wrong data) stops at order 2, whatever its own order.
use gideon
implicit none
integer(I4P), parameter :: NGRID = 4                 ! grids of 16, 32, 64, 128 cells
real(R_P),    parameter :: PI = acos(-1._R_P)
real(R_P)               :: h, x, d, err(3,NGRID)     ! step, a centre, a derivative, max error per case and grid
real(R_P)               :: p(-3:3), a(-3:3)          ! point values and averages around one cell, by offset
integer(I4P)            :: n, i, k, o

err = 0._R_P
do k=1, NGRID
   n = 8 * 2**k
   h = 1._R_P / n
   do i=1, n
      x = (i - 0.5_R_P)*h                            ! centre of cell i (periodic: no ghosts needed, u is known)
      do o=-3, 3
         p(o) = sin(2*PI*(x + o*h))                                        ! point values at the centres
         a(o) = (cos(2*PI*(x + (o - 0.5_R_P)*h)) - cos(2*PI*(x + (o + 0.5_R_P)*h)))/(2*PI*h) ! cell averages
      enddo
      ! FD kernel on point values: the derivative at the centre
      call gideon_fd_d1_c6(p, h, d)
      err(1,k) = max(err(1,k), abs(d - 2*PI*cos(2*PI*x)))
      ! FV kernel on cell averages: the average of the derivative over the cell
      call gideon_fv_d1_c6(a, h, d)
      err(2,k) = max(err(2,k), abs(d - (sin(2*PI*(x + h/2)) - sin(2*PI*(x - h/2)))/h))
      ! FD kernel on cell averages, read as point values: the wrong data
      call gideon_fd_d1_c6(a, h, d)
      err(3,k) = max(err(3,k), abs(d - 2*PI*cos(2*PI*x)))
   enddo
enddo

print '(A)', 'd/dx sin(2 pi x), sixth-order centred kernels, max error over the cells'
print '(A)', '     n  fd on points  fv on averages  fd on averages'
do k=1, NGRID
   print '(I6,ES13.3,ES16.3,ES16.3)', 8 * 2**k, err(:,k)
enddo
print '(A)', 'observed orders, log2(e(n/2) / e(n))'
do k=2, NGRID
   print '(I6,F13.2,F16.2,F16.2)', 8 * 2**k, log(err(:,k-1) / err(:,k)) / log(2._R_P)
enddo
endprogram cookbook_choose
