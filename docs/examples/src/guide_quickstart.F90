program guide_quickstart
!< Quick start: the first derivative of sin x at the nodes of a periodic grid, a sixth-order centred kernel, then the
!< same through the run-time order dispatcher.
use gideon
implicit none
integer(I4P), parameter :: N = 32, S = 3             ! nodes 0..N-1; the half-width of the window q(-3:3)
real(R_P),    parameter :: PI = acos(-1._R_P)
real(R_P)               :: u(-S:N-1+S), du(0:N-1), dr(0:N-1), h
integer(I4P)            :: i, order

h = 2*PI / N
u = [(sin(i*h), i=-S, N-1+S)]                         ! the ghosts by periodicity: sin is known everywhere
do i=0, N-1
   call gideon_fd_d1_c6(u(i-3:i+3), h, du(i))        ! node i at offset 0 of the window
enddo
order = 6                                            ! known only at run time, in general
do i=0, N-1
   call gideon_fd_d1_c(order, u(i-order/2:i+order/2), h, dr(i))
enddo
print '(A,ES10.3)', 'gideon_fd_d1_c6, max error: ', maxval(abs(du - [(cos(i*h), i=0, N-1)]))
print '(A,L1)',     'gideon_fd_d1_c(6, ...) gives the same values, bitwise: ', all(dr == du)
endprogram guide_quickstart
