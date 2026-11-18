program cookbook_dissipation
!< Cookbook: dissipation dialled by a sensor. A smooth bump and a square wave advected once around a periodic line
!< (u_t + u_x = 0, cell averages, SSP-RK3) with four face values: centred of order 6; the centred one plus the dissipation part of the
!< upwind-biased bl5 (that is bl5 itself); the same with the dissipation scaled by a sensor; and Jameson-Schmidt-Turkel's
!< blend of the second and fourth face differences on the centred value of order 2.
use gideon
implicit none
integer(I4P), parameter :: N = 100, G = 3       ! cells, ghost cells
real(R_P),    parameter :: CFL = 0.4_R_P
character(*), parameter :: NAMES(4) = [character(16) :: 'central c6', 'biased bl5', 'sensor c6 + s d', 'JST on c2']
real(R_P)    :: u(1-G:N+G), u0(N), u1(N), un(N), x, h, dt, mass0
integer(I4P) :: scheme, i, step, nsteps

h = 1._R_P/N
dt = CFL*h
nsteps = nint(1._R_P/dt)
dt = 1._R_P/nsteps                                ! exactly one period
do i=1, N ! a Gaussian bump on [0, 1/2), a square wave on [0.6, 0.85]: cell averages (the bump's by a 5-point Gauss rule)
   x = (i - 0.5_R_P)*h
   u0(i) = merge(1._R_P, 0._R_P, x > 0.6_R_P .and. x < 0.85_R_P) + bump_average(x - h/2, x + h/2)
enddo
print '(A)', 'a bump and a square wave advected once around 100 cells, CFL 0.4, SSP-RK3'
print '(A)', 'scheme           bump error  square overshoot  undershoot  mass conserved'
do scheme=1, 4
   u(1:N) = u0
   mass0 = sum(u0)
   do step=1, nsteps                              ! SSP-RK3 (Shu-Osher)
      u1 = u(1:N)
      call update(scheme, u, dt, un)
      u(1:N) = un
      call update(scheme, u, dt, un)
      u(1:N) = 0.75_R_P*u1 + 0.25_R_P*un
      call update(scheme, u, dt, un)
      u(1:N) = u1/3 + 2*un/3
   enddo
   print '(A16,F11.5,F18.5,F12.5,L16)', NAMES(scheme), sum(abs(u(1:N/2) - u0(1:N/2)))*h, &
      maxval(u(N/2+1:N)) - 1, max(0._R_P, -minval(u(N/2+1:N))), abs(sum(u(1:N)) - mass0) < 1.e-12_R_P*N
enddo
contains
   pure function bump_average(a, b) result(v)
   !< The average over [a, b] of exp(-((x - 1/4)/0.06)^2), a 5-point Gauss-Legendre rule (exact to 1e-14 here).
   real(R_P), intent(in) :: a, b
   real(R_P)             :: v
   real(R_P), parameter  :: XG(5) = [-0.9061798459386640_R_P, -0.5384693101056831_R_P, 0._R_P, &
                                      0.5384693101056831_R_P, 0.9061798459386640_R_P]
   real(R_P), parameter  :: WG(5) = [0.2369268850561891_R_P, 0.4786286704993665_R_P, 0.5688888888888889_R_P, &
                                      0.4786286704993665_R_P, 0.2369268850561891_R_P]
   v = sum(WG*exp(-(((a + b)/2 + (b - a)/2*XG - 0.25_R_P)/0.06_R_P)**2))/2
   endfunction bump_average

   subroutine update(scheme, u, dt, unew)
   !< One forward Euler step of the flux difference: unew = u - dt/h (f(i+1/2) - f(i-1/2)).
   integer(I4P), intent(in)    :: scheme
   real(R_P),    intent(inout) :: u(1-G:N+G)
   real(R_P),    intent(in)    :: dt
   real(R_P),    intent(out)   :: unew(N)
   real(R_P)                   :: f(0:N), c, d, s, d1, d3, e2, e4, nu(0:N+1)
   integer(I4P)                :: i
   u(1-G:0) = u(N-G+1:N) ; u(N+1:N+G) = u(1:G)   ! periodic ghosts
   do i=0, N+1                                    ! a pressure-like switch: the normalised second difference
      nu(i) = abs(u(i+1) - 2*u(i) + u(i-1))/(abs(u(i+1)) + 2*abs(u(i)) + abs(u(i-1)) + 1.e-12_R_P)
   enddo
   do i=0, N                                      ! the face i+1/2, between cells i and i+1
      select case(scheme)
      case(1)
         call gideon_fv_face_c6(u(i-2:i+3), f(i))
      case(2, 3)
         call gideon_fv_face_c6(u(i-2:i+3), c)
         call gideon_fv_face_diss_bl5(u(i-2:i+3), d) ! bl5 = c6 + d, exactly in rationals
         s = 1._R_P
         if (scheme == 3) s = min(1._R_P, 4*max(nu(i), nu(i+1)))
         f(i) = c + s*d
      case(4)
         call gideon_fv_face_c2(u(i:i+1), c)
         call gideon_fv_face_diss_c2(u(i:i+1), d1)     ! delta^1 = u(i+1) - u(i)
         call gideon_fv_face_diss_c4(u(i-1:i+2), d3)   ! delta^3
         e2 = 0.5_R_P*max(nu(i), nu(i+1))              ! Jameson's k2 = 1/2, k4 = 1/32, the wave speed 1
         e4 = max(0._R_P, 1._R_P/32 - e2)
         f(i) = c - (e2*d1 - e4*d3)
      endselect
   enddo
   do i=1, N
      unew(i) = u(i) - dt/h*(f(i) - f(i-1))
   enddo
   endsubroutine update
endprogram cookbook_dissipation
