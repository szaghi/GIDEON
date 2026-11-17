program cookbook_run_time_order
!< Cookbook, a run-time order: the order read at run time, validated once, then either dispatched at every point
!< (gideon_fd_d1_c, one select case per call) or chosen once outside the loop (one loop per order, each calling its
!< kernel). Both give the same values, bitwise; the second pays no dispatch inside the loop.
use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
use gideon
implicit none
integer(I4P), parameter :: N = 64, G = 5            ! nodes 0..N-1, periodic; ghosts: the widest window, order 10
integer(I4P), parameter :: ORDERS(6) = [2, 3, 4, 6, 8, 10]
real(R_P),    parameter :: PI = acos(-1._R_P)
real(R_P)               :: u(-G:N-1+G), dd(0:N-1), dk(0:N-1), probe(-G:G), t, h
integer(I4P)            :: order, i, k

h = 2*PI / N
u(0:N-1) = [(exp(sin(i*h)), i=0, N-1)]
u(-G:-1) = u(N-G:N-1) ; u(N:N-1+G) = u(0:G-1)
probe = 0._R_P
do k=1, size(ORDERS)
   order = ORDERS(k)                                 ! in a solver: read from the input
   ! validate once: an order without a kernel returns GIDEON_NAN (a quiet NaN), never a number
   call gideon_fd_d1_c(order, probe(-(order/2):order/2), 1._R_P, t)
   if (ieee_is_nan(t)) then
      print '(A,I0,A)', 'order ', order, ': no kernel (GIDEON_NAN), refused before the loop'
      cycle
   endif
   ! (a) the dispatcher at every node
   do i=0, N-1
      call gideon_fd_d1_c(order, u(i-order/2:i+order/2), h, dd(i))
   enddo
   ! (b) the choice hoisted out of the loop: one loop per order, the kernel inside
   select case(order)
   case(2)  ; do i=0, N-1 ; call gideon_fd_d1_c2(u(i-1:i+1), h, dk(i)) ; enddo
   case(4)  ; do i=0, N-1 ; call gideon_fd_d1_c4(u(i-2:i+2), h, dk(i)) ; enddo
   case(6)  ; do i=0, N-1 ; call gideon_fd_d1_c6(u(i-3:i+3), h, dk(i)) ; enddo
   case(8)  ; do i=0, N-1 ; call gideon_fd_d1_c8(u(i-4:i+4), h, dk(i)) ; enddo
   case(10) ; do i=0, N-1 ; call gideon_fd_d1_c10(u(i-5:i+5), h, dk(i)) ; enddo
   endselect
   print '(A,I0,A,L1,A,ES10.3)', 'order ', order, ': dispatched = hoisted bitwise ', all(dd == dk), &
                                 ', max error ', maxval(abs(dk - [(cos(i*h)*exp(sin(i*h)), i=0, N-1)]))
enddo
endprogram cookbook_run_time_order
