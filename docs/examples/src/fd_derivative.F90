program fd_derivative
!< Tutorial T1: the first derivative of u = exp(sin x) on a periodic grid, centred kernels of order 2 to 8, and the
!< orders observed by halving the step.
use gideon
implicit none
integer(I4P), parameter :: G = 4                 ! ghost nodes on each side: the half-width of the widest window
integer(I4P), parameter :: NGRID = 4             ! number of grids: 16, 32, 64, 128 nodes
real(R_P),    parameter :: PI = acos(-1._R_P)
real(R_P), allocatable  :: x(:), u(:), du(:,:)   ! nodes, point values (with ghosts), derivatives (one column per order)
real(R_P)               :: h, err(4,NGRID)       ! step, max error per order and grid
integer(I4P)            :: n, i, k, o

do k=1, NGRID
   n = 8 * 2**k
   h = 2*PI / n
   allocate(x(-G:n-1+G), u(-G:n-1+G), du(0:n-1,4))
   x = [(i*h, i=-G, n-1+G)]
   u(0:n-1) = exp(sin(x(0:n-1)))                 ! the n distinct nodes of the period [0, 2 pi)
   u(-G:-1) = u(n-G:n-1)                          ! ghosts on the left: the last G nodes
   u(n:n-1+G) = u(0:G-1)                          ! ghosts on the right: the first G nodes
   do i=0, n-1
      call gideon_fd_d1_c2(u(i-1:i+1), h, du(i,1)) ! window q(-1:1): node i at offset 0
      call gideon_fd_d1_c4(u(i-2:i+2), h, du(i,2)) ! window q(-2:2)
      call gideon_fd_d1_c6(u(i-3:i+3), h, du(i,3)) ! window q(-3:3)
      call gideon_fd_d1_c8(u(i-4:i+4), h, du(i,4)) ! window q(-4:4)
   enddo
   do o=1, 4
      err(o,k) = maxval(abs(du(:,o) - cos(x(0:n-1))*u(0:n-1)))
   enddo
   deallocate(x, u, du)
enddo

print '(A)', 'd/dx exp(sin x) on [0, 2 pi), max error over the nodes'
print '(A)', '     n       c2         c4         c6         c8'
do k=1, NGRID
   print '(I6,4ES11.3)', 8 * 2**k, err(:,k)
enddo
print '(A)', 'observed orders, log2(e(n/2) / e(n))'
do k=2, NGRID
   print '(I6,4F11.2)', 8 * 2**k, log(err(:,k-1) / err(:,k)) / log(2._R_P)
enddo
endprogram fd_derivative
