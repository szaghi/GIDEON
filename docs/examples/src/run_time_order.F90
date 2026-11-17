program run_time_order
!< Tutorial T4: the order chosen when the program runs. The dispatcher gideon_fd_d1_c over orders 1 to 12, its window
!< bounds as expressions of the order, the NaN of an order without a kernel, and the same values from a choice made
!< once, outside the loop.
use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
use gideon
implicit none
integer(I4P), parameter :: N = 64                ! nodes
integer(I4P), parameter :: G = 6                 ! ghost nodes: half of the largest order tried
real(R_P),    parameter :: PI = acos(-1._R_P)
real(R_P)               :: x(-G:N-1+G), u(-G:N-1+G), du(0:N-1), dk(0:N-1), h
integer(I4P)            :: order, s, i

h = 2*PI / N
x = [(i*h, i=-G, N-1+G)]
u(0:N-1) = exp(sin(x(0:N-1)))
u(-G:-1) = u(N-G:N-1)
u(N:N-1+G) = u(0:G-1)

print '(A)', 'gideon_fd_d1_c(order, u(i-s:i+s), h, d), s = order/2, d/dx exp(sin x), 64 nodes'
print '(A)', 'order  window        max error'
do order=1, 12
   s = order / 2
   do i=0, N-1
      call gideon_fd_d1_c(order, u(i-s:i+s), h, du(i)) ! the window's bounds follow the order
   enddo
   if (any(ieee_is_nan(du))) then
      print '(I5,A,I0,A,I0,A,T22,A)', order, '  q(', -s, ':', s, ')', 'no kernel: GIDEON_NAN'
   else
      print '(I5,A,I0,A,I0,A,T21,ES10.3)', order, '  q(', -s, ':', s, ')', maxval(abs(du - cos(x(0:N-1))*u(0:N-1)))
   endif
enddo

! the dispatcher calls the kernel: the same values, bitwise
order = 6
do i=0, N-1
   call gideon_fd_d1_c(order, u(i-3:i+3), h, du(i))
   call gideon_fd_d1_c6(u(i-3:i+3), h, dk(i))
enddo
print '(A,I0,A,I0)', 'gideon_fd_d1_c(6, ...) == gideon_fd_d1_c6: ', count(du == dk), ' of ', N

! the choice made once: one loop per order, the kernel inside (no select case per point)
select case(order)
case(2)
   do i=0, N-1 ; call gideon_fd_d1_c2(u(i-1:i+1), h, dk(i)) ; enddo
case(4)
   do i=0, N-1 ; call gideon_fd_d1_c4(u(i-2:i+2), h, dk(i)) ; enddo
case(6)
   do i=0, N-1 ; call gideon_fd_d1_c6(u(i-3:i+3), h, dk(i)) ; enddo
case default
   error stop 'run_time_order: order not supported by this program'
endselect
print '(A,I0,A,I0)', 'select case outside the loop == dispatcher:  ', count(du == dk), ' of ', N
endprogram run_time_order
