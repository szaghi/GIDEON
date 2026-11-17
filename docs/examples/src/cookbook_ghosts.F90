program cookbook_ghosts
!< Cookbook, ghost cells: the halo of an array derived from the windows of the kernels it feeds, filled periodically.
!< A node loop (gideon_fd_d1_c6, window q(-3:3)) and a face loop of upwind values (gideon_fv_face_ul3, q(-2:0), and
!< gideon_fv_face_ur3, q(1:3)), whose flux difference is the cell kernel gideon_fv_d1_ul3, bitwise.
use gideon
implicit none
integer(I4P), parameter :: N = 32                ! nodes or cells 1..N
real(R_P),    parameter :: PI = acos(-1._R_P)
! node loop: the window of gideon_fd_d1_c6 is q(-3:3); node i reads i-3..i+3, so nodes 1..N read -2..N+3
integer(I4P), parameter :: D1_LO = -3, D1_HI = 3
real(R_P)               :: u(1+D1_LO:N+D1_HI), du(N)
! face loop: face i+1/2 (i = 0..N) reads cells i-2..i (left side) or i+1..i+3 (right side): cells -2..N+3
integer(I4P), parameter :: UL_LO = -2, UR_HI = 3
real(R_P)               :: a(UL_LO:N+UR_HI), fl(0:N), fr(0:N), dl
real(R_P)               :: h, x
integer(I4P)            :: i, same

h = 1._R_P / N
! point values at the nodes x = (i - 1) h, then the ghosts by periodicity
do i=1, N
   u(i) = sin(2*PI*(i - 1)*h)
enddo
u(1+D1_LO:0) = u(N+1+D1_LO:N)                    ! left ghosts: the last -D1_LO nodes
u(N+1:N+D1_HI) = u(1:D1_HI)                      ! right ghosts: the first D1_HI nodes
do i=1, N
   call gideon_fd_d1_c6(u(i+D1_LO:i+D1_HI), h, du(i))
enddo
print '(A,I0,A,I0,A,ES10.3)', 'nodes: halo ', -D1_LO, ' + ', D1_HI, ', max error of d1 ', &
                              maxval(abs(du - [(2*PI*cos(2*PI*(i - 1)*h), i=1, N)]))

! cell averages over [(i - 1) h, i h], then the ghosts by periodicity
do i=1, N
   a(i) = (cos(2*PI*(i - 1)*h) - cos(2*PI*i*h))/(2*PI*h)
enddo
a(UL_LO:0) = a(N+UL_LO:N)
a(N+1:N+UR_HI) = a(1:UR_HI)
do i=0, N
   call gideon_fv_face_ul3(a(i-2:i), fl(i))      ! the upwind value for a flow to the right
   call gideon_fv_face_ur3(a(i+1:i+3), fr(i))    ! the upwind value for a flow to the left
enddo
x = 0._R_P
same = 0
do i=1, N
   x = max(x, abs(fl(i) - sin(2*PI*i*h)), abs(fr(i) - sin(2*PI*i*h)))
   call gideon_fv_d1_ul3(a(i-3:i), h, dl)        ! the cell kernel reads cells i-3..i: the same halo on the left
   if (dl == (fl(i) - fl(i-1))/h) same = same + 1
enddo
print '(A,I0,A,I0,A,ES10.3)', 'faces: halo ', -UL_LO + 1, ' + ', UR_HI, ', max error of the face values ', x
print '(A,I0,A,I0)', 'flux difference of the left-side faces = gideon_fv_d1_ul3, bitwise, in cells: ', same, ' of ', N
endprogram cookbook_ghosts
