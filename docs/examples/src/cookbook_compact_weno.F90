program cookbook_compact_weno
!< Cookbook: compact-WENO rows. Face values of the cell averages of a smooth bump and a square wave on a periodic line:
!< the explicit upwind-biased bl5, the compact upwind-biased pbl5 on every face, and pbl5 with the faces a sensor flags
!< replaced by a given value (here the first-order upwind one, the stand-in for a WENO/TENO value): one masked solve.
use gideon
implicit none
integer(I4P), parameter :: N = 100, G = 2
real(R_P)               :: u(1-G:N+G), f(0:N,3), w(N,2), x, h, nu(0:N+1), exact
logical                 :: mask(N)
integer(I4P)            :: i, s
character(*), parameter :: NAMES(3) = [character(24) :: 'explicit bl5', 'compact pbl5', 'compact pbl5, flagged']

h = 1._R_P/N
do i=1, N
   x = (i - 0.5_R_P)*h
   u(i) = merge(1._R_P, 0._R_P, x > 0.6_R_P .and. x < 0.85_R_P) + bump_average(x - h/2, x + h/2)
enddo
u(1-G:0) = u(N-G+1:N) ; u(N+1:N+G) = u(1:G)
do i=0, N+1 ! the sensor: the second difference of the averages over the line's scale (a jump: 1; the bump: < 0.06)
   nu(i) = abs(u(i+1) - 2*u(i) + u(i-1))/maxval(abs(u(1:N)))
enddo
do i=1, N   ! face i is between cells i and i+1 (face N: between N and 1)
   call gideon_fv_face_bl5(u(i-2:i+2), f(i,1))
   call gideon_fv_face_pbl5(u(i-1:i+1), f(i,2))
   mask(i) = max(nu(i), nu(i+1)) > 0.25_R_P
   f(i,3) = merge(u(i), f(i,2), mask(i))  ! a flagged face: its row the identity, its right side the value given
enddo
call gideon_tri_periodic_const(N, FV_FACE_PBL5_LEFT(-1), FV_FACE_PBL5_LEFT(0), FV_FACE_PBL5_LEFT(1), f(1:N,2), &
                               w(:,1), w(:,2))
call gideon_tri_periodic_masked(N, FV_FACE_PBL5_LEFT(-1), FV_FACE_PBL5_LEFT(0), FV_FACE_PBL5_LEFT(1), mask, f(1:N,3), &
                                w(:,1), w(:,2))
print '(A,I0,A)', 'face values of a bump and a square wave on 100 cells; ', count(mask), ' faces flagged by the sensor'
print '(A)', 'scheme                     overshoot  undershoot  bump error'
do s=1, 3
   exact = 0._R_P
   do i=5, 45 ! the faces of the bump, far from the square wave
      exact = max(exact, abs(f(i,s) - exp(-((i*h - 0.25_R_P)/0.06_R_P)**2)))
   enddo
   print '(A24,3F12.5)', NAMES(s), max(0._R_P, maxval(f(1:N,s)) - 1), max(0._R_P, -minval(f(1:N,s))), exact
enddo
contains
   pure function bump_average(a, b) result(v)
   !< The average over [a, b] of exp(-((x - 1/4)/0.06)^2), a 5-point Gauss-Legendre rule.
   real(R_P), intent(in) :: a, b
   real(R_P)             :: v
   real(R_P), parameter  :: XG(5) = [-0.9061798459386640_R_P, -0.5384693101056831_R_P, 0._R_P, &
                                      0.5384693101056831_R_P, 0.9061798459386640_R_P]
   real(R_P), parameter  :: WG(5) = [0.2369268850561891_R_P, 0.4786286704993665_R_P, 0.5688888888888889_R_P, &
                                      0.4786286704993665_R_P, 0.2369268850561891_R_P]
   v = sum(WG*exp(-(((a + b)/2 + (b - a)/2*XG - 0.25_R_P)/0.06_R_P)**2))/2
   endfunction bump_average
endprogram cookbook_compact_weno
