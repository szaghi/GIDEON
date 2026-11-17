program upwind
!< Tutorial T3: upwind stencils. The mirror identity between left- and right-side kernels, checked bitwise, and the
!< upwind face values of the advection u_t + a u_x = 0 for both signs of a.
use gideon
implicit none
integer(I4P), parameter :: M = 1000              ! windows for the mirror check
integer(I4P), parameter :: G = 3                 ! ghost cells on each side: the reach of the order-3 faces
integer(I4P), parameter :: NGRID = 4             ! number of grids: 16, 32, 64, 128 cells
real(R_P),    parameter :: PI = acos(-1._R_P)
real(R_P)               :: q(-10:M+10)           ! arbitrary data for the mirror check
real(R_P)               :: fl, fr, dl, dr        ! left- and right-side results
real(R_P), allocatable  :: xf(:), u(:), f(:), res(:) ! faces, cell averages (with ghosts), face values, residuals
real(R_P)               :: h, a, err(2,NGRID)
integer(I4P)            :: i, k, n, s, nface, nd1, nd2

! the mirror: a left-side kernel on q equals its right-side twin on q reversed, times (-1)^n for the n-th derivative
q = [(sin(0.7_R_P*i + 0.3_R_P) + 0.01_R_P*i, i=-10, M+10)]
nface = 0 ; nd1 = 0 ; nd2 = 0
do i=1, M
   call gideon_fv_face_ul3(q(i-2:i), fl)          ! the face i+1/2 from cells i-2, i-1, i
   call gideon_fv_face_ur3(q(i:i-2:-1), fr)       ! the same cells, reversed: r(m) = q(i+1-m), m = 1..3
   if (fl == fr) nface = nface + 1
   call gideon_fv_d1_ul3(q(i-3:i), 1._R_P, dl)   ! d1 averaged over cell i, cells i-3..i
   call gideon_fv_d1_ur3(q(i:i-3:-1), 1._R_P, dr) ! reversed: r(m) = q(i-m), m = 0..3
   if (dl == -dr) nd1 = nd1 + 1
   call gideon_fv_d2_ul3(q(i-4:i), 1._R_P, dl)
   call gideon_fv_d2_ur3(q(i:i-4:-1), 1._R_P, dr)
   if (dl == dr) nd2 = nd2 + 1
enddo
print '(A)', 'mirror identity, bitwise, over 1000 windows'
print '(A,I0)', '  face value:  ul3(q) ==  ur3(reversed q)  ', nface
print '(A,I0)', '  d1 average:  ul3(q) == -ur3(reversed q)  ', nd1
print '(A,I0)', '  d2 average:  ul3(q) ==  ur3(reversed q)  ', nd2

! advection: the residual -a (f(i+1/2) - f(i-1/2))/h, f the face value from the upwind side
print '(A)', 'residual of u_t + a u_x = 0, u = exp(sin x), order-3 upwind faces, max error'
print '(A)', '     n     a = +1      a = -1'
do k=1, NGRID
   n = 8 * 2**k
   h = 2*PI / n
   allocate(xf(-1:n-1), u(-G:n-1+G), f(-1:n-1), res(0:n-1))
   xf = [(i*h + h, i=-1, n-1)]                    ! xf(i): the face i+1/2, right of cell i
   u(0:n-1) = cell_average(xf(-1:n-2), xf(0:n-1))
   u(-G:-1) = u(n-G:n-1)                          ! periodic ghost cells
   u(n:n-1+G) = u(0:G-1)
   do s=1, 2
      a = merge(1._R_P, -1._R_P, s == 1)
      do i=0, n-1
         if (a > 0._R_P) then
            call gideon_fv_face_ul3(u(i-2:i), f(i)) ! the flow comes from the left: cells i-2..i
         else
            call gideon_fv_face_ur3(u(i+1:i+3), f(i)) ! the flow comes from the right: cells i+1..i+3
         endif
      enddo
      f(-1) = f(n-1)
      res = -a * (f(0:n-1) - f(-1:n-2)) / h
      err(s,k) = maxval(abs(res + a * (exp(sin(xf(0:n-1))) - exp(sin(xf(-1:n-2)))) / h))
   enddo
   print '(I6,2ES12.3)', n, err(:,k)
   deallocate(xf, u, f, res)
enddo
print '(A)', 'observed orders, log2(e(n/2) / e(n))'
do k=2, NGRID
   print '(I6,2F12.2)', 8 * 2**k, log(err(:,k-1) / err(:,k)) / log(2._R_P)
enddo
contains
   elemental function cell_average(xa, xb) result(avg)
   !< The average of exp(sin x) over [xa, xb]: 8-point Gauss-Legendre per cell, exact far below the errors above.
   real(R_P), intent(in) :: xa, xb
   real(R_P)             :: avg
   real(R_P), parameter  :: GX(4) = [0.1834346424956498_R_P, 0.5255324099163290_R_P, &
                                     0.7966664774136267_R_P, 0.9602898564975363_R_P]
   real(R_P), parameter  :: GW(4) = [0.3626837833783620_R_P, 0.3137066458778873_R_P, &
                                     0.2223810344533745_R_P, 0.1012285362903763_R_P]
   real(R_P)             :: c, r
   c = (xa + xb) / 2 ; r = (xb - xa) / 2
   avg = sum(GW * (exp(sin(c - r*GX)) + exp(sin(c + r*GX)))) / 2
   endfunction cell_average
endprogram upwind
