program cookbook_deconvolution
!< Cookbook, averages and point values: the cell average of a nonlinear flux f(u) = u^2/2 from the cell averages of
!< u = sin(2 pi x). f(<u>) is second order whatever the reconstruction that follows; deconvolving to point values
!< (gideon_fv_point_c<p>), evaluating f there and convolving back (gideon_fd_average_c<p>) is of order p.
use gideon
implicit none
integer(I4P), parameter :: NGRID = 4, G = 4          ! grids of 16..128 cells; ghosts: two windows of half-width 2
real(R_P),    parameter :: PI = acos(-1._R_P)
real(R_P), allocatable  :: a(:), fp4(:), fp6(:)      ! averages of u; f at the point values, orders 4 and 6
real(R_P)               :: h, v, fex, err(3,NGRID)
integer(I4P)            :: n, i, k

err = 0._R_P
do k=1, NGRID
   n = 8 * 2**k
   h = 1._R_P / n
   allocate(a(1-G:n+G), fp4(-1:n+2), fp6(-1:n+2))
   do i=1-G, n+G                                     ! cell i spans [(i - 1) h, i h]; periodic, so the ghosts are known
      a(i) = (cos(2*PI*(i - 1)*h) - cos(2*PI*i*h))/(2*PI*h)
   enddo
   do i=-1, n+2                                      ! point values at the centres, then f there
      call gideon_fv_point_c4(a(i-1:i+1), v) ; fp4(i) = v**2/2
      call gideon_fv_point_c6(a(i-2:i+2), v) ; fp6(i) = v**2/2
   enddo
   do i=1, n
      ! exact <f> = (1 - <cos(4 pi x)>)/4 over the cell
      fex = (1 - (sin(4*PI*i*h) - sin(4*PI*(i - 1)*h))/(4*PI*h))/4
      err(1,k) = max(err(1,k), abs(a(i)**2/2 - fex))  ! f(<u>)
      call gideon_fd_average_c4(fp4(i-1:i+1), v) ; err(2,k) = max(err(2,k), abs(v - fex))
      call gideon_fd_average_c6(fp6(i-2:i+2), v) ; err(3,k) = max(err(3,k), abs(v - fex))
   enddo
   deallocate(a, fp4, fp6)
enddo

print '(A)', 'cell average of u^2/2, u = sin(2 pi x), max error over the cells'
print '(A)', '     n     f(<u>)   order 4   order 6'
do k=1, NGRID
   print '(I6,3ES10.2)', 8 * 2**k, err(:,k)
enddo
print '(A)', 'observed orders, log2(e(n/2) / e(n))'
do k=2, NGRID
   print '(I6,3F10.2)', 8 * 2**k, log(err(:,k-1) / err(:,k)) / log(2._R_P)
enddo
endprogram cookbook_deconvolution
