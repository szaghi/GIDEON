program compact_derivative
!< Tutorial T7: the first derivative of u = exp(sin x) on a periodic line, compact (Pade) operators of order 4 to 10
!< against the explicit centred one of order 6: the right side point by point, then one solve of the band.
use gideon
implicit none
integer(I4P), parameter :: G = 3                 ! ghost nodes on each side: the half-width of the widest right side
integer(I4P), parameter :: NGRID = 4             ! number of grids: 16, 32, 64, 128 nodes
real(R_P),    parameter :: PI = acos(-1._R_P)
real(R_P), allocatable  :: x(:), u(:), du(:,:), w(:,:) ! nodes, point values (with ghosts), derivatives, scratch
real(R_P)               :: h, err(5,NGRID)
integer(I4P)            :: n, i, k, o

do k=1, NGRID
   n = 8 * 2**k
   h = 2*PI / n
   allocate(x(-G:n-1+G), u(-G:n-1+G), du(0:n-1,5), w(0:n-1,4))
   x = [(i*h, i=-G, n-1+G)]
   u(0:n-1) = exp(sin(x(0:n-1)))
   u(-G:-1) = u(n-G:n-1) ; u(n:n-1+G) = u(0:G-1)  ! periodic ghosts
   do i=0, n-1                                    ! the right sides, one kernel call per node
      call gideon_fd_d1_p4(u(i-1:i+1), h, du(i,1))
      call gideon_fd_d1_p6(u(i-2:i+2), h, du(i,2))
      call gideon_fd_d1_p8(u(i-2:i+2), h, du(i,3))
      call gideon_fd_d1_p10(u(i-3:i+3), h, du(i,4))
      call gideon_fd_d1_c6(u(i-3:i+3), h, du(i,5)) ! the explicit operator of order 6, for comparison
   enddo
   ! the bands: tridiagonal for orders 4 and 6, pentadiagonal for 8 and 10 (their pairs: index j at offsets +-j)
   call gideon_tri_periodic_const(n, FD_D1_P4_LEFT(1), FD_D1_P4_LEFT(0), FD_D1_P4_LEFT(1), du(:,1), w(:,1), w(:,2))
   call gideon_tri_periodic_const(n, FD_D1_P6_LEFT(1), FD_D1_P6_LEFT(0), FD_D1_P6_LEFT(1), du(:,2), w(:,1), w(:,2))
   call gideon_penta_periodic_const(n, FD_D1_P8_LEFT(2), FD_D1_P8_LEFT(1), FD_D1_P8_LEFT(0), FD_D1_P8_LEFT(1), &
                                    FD_D1_P8_LEFT(2), du(:,3), w(:,1), w(:,2), w(:,3), w(:,4))
   call gideon_penta_periodic_const(n, FD_D1_P10_LEFT(2), FD_D1_P10_LEFT(1), FD_D1_P10_LEFT(0), FD_D1_P10_LEFT(1), &
                                    FD_D1_P10_LEFT(2), du(:,4), w(:,1), w(:,2), w(:,3), w(:,4))
   do o=1, 5
      err(o,k) = maxval(abs(du(:,o) - cos(x(0:n-1))*u(0:n-1)))
   enddo
   deallocate(x, u, du, w)
enddo

print '(A)', 'd/dx exp(sin x) on [0, 2 pi), max error over the nodes'
print '(A)', '     n        p4         p6         p8        p10         c6'
do k=1, NGRID
   print '(I6,5ES11.3)', 8 * 2**k, err(:,k)
enddo
print '(A)', 'observed orders, log2(e(n/2) / e(n))'
do k=2, NGRID
   print '(I6,5F11.2)', 8 * 2**k, log(err(:,k-1) / err(:,k)) / log(2._R_P)
enddo
endprogram compact_derivative
