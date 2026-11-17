!< GIDEON verification (A3): linear advection of a wave packet with centred and upwind first derivatives.
program gideon_vv_advection
!< GIDEON verification (A3): linear advection of a wave packet with centred and upwind first derivatives.
!<
!< u_t + u_x = 0 on the periodic [0, 1), N = 128 cells, cell averages advanced by RK4 at CFL 0.25 for two periods,
!< the derivative du/dx averaged over each cell given by a GIDEON kernel: gideon_fv_d1_c2 and gideon_fv_d1_c4
!< (centred), gideon_fv_d1_ul2 (upwind for this flow, the stencil on the left side). The packet
!< u0 = exp(-((x - 1/2)/0.07)^2) sin(16 pi x) has its carrier at xi = k h = pi/8: there the centred kernels disperse
!< (the packet travels at their group velocity, below 1) and do not dissipate, the upwind one dissipates (its
!< modified wavenumber, A3). The exact solution is u0(x - t), averaged over the cells by 5-point Gauss-Legendre
!< quadrature, as the initial data.
!<
!< Frozen bounds (written before the run): the cell averages' sum (the mass) is kept within 1e3 eps sum|u0| by every
!< kernel (flux form); the discrete energy sum u^2 never grows (RK4 at CFL 0.25 inside the stability region of all
!< three); after two periods the upwind kernel has lost at least half of the energy and the fourth-order centred one
!< has a smaller error than the second-order one.
!<
!< Report mode (scripts/run_report.sh): with GIDEON_REPORT_DIR set, the frames of the animation are written to
!< advection-frames-<kernel>.dat there: the first line the cell centres, then per frame t, u(1:N), exact(1:N). The
!< animation of the documentation (docs/scripts/make_animations.py) is drawn from them: GIDEON's compiled kernels,
!< not a re-implementation.
use penf, only : I4P, R8P, R_P
use gideon
implicit none
integer(I4P), parameter :: N = 128, G = 5, NK = 3, PERIODS = 2, EVERY = 8
real(R_P),    parameter :: CFL = 0.25_R_P, H = 1._R_P/N, DT = CFL*H
real(R_P),    parameter :: X0 = 0.5_R_P, SIGMA = 0.07_R_P
real(R_P),    parameter :: PI = 3.14159265358979323846264338327950288_R_P
character(*), parameter :: NAMES(NK) = ['fv_d1_c2 ', 'fv_d1_c4 ', 'fv_d1_ul2']
real(R_P)    :: u(N), u0(N), k1(N), k2(N), k3(N), k4(N), ex(N), xc(N)
real(R_P)    :: mass0, energy0, energy, energy_prev, err(NK), energy_end(NK), t
integer(I4P) :: kk, step, nsteps, i, fails, uf
character(4096) :: report_dir
integer      :: rlen
logical      :: report

call get_environment_variable('GIDEON_REPORT_DIR', report_dir, length=rlen)
report = rlen > 0
fails = 0
do i=1, N
   xc(i) = (i - 0.5_R_P)*H
enddo
call averages(0._R_P, u0)
mass0 = sum(u0) ; energy0 = sum(u0**2)
nsteps = nint(PERIODS/DT)
do kk=1, NK
   u = u0
   energy_prev = energy0
   t = 0._R_P
   if (report) then
      open(newunit=uf, file=trim(report_dir)//'/advection-frames-'//trim(NAMES(kk))//'.dat', action='write', &
           status='replace')
      write(uf, '(*(1X,ES24.16E3))') real(xc, R8P)
      write(uf, '(*(1X,ES24.16E3))') real(t, R8P), real(u, R8P), real(u0, R8P)
   endif
   do step=1, nsteps
      call rhs(kk, u, k1)
      call rhs(kk, u + 0.5_R_P*DT*k1, k2)
      call rhs(kk, u + 0.5_R_P*DT*k2, k3)
      call rhs(kk, u + DT*k3, k4)
      u = u + DT*(k1 + 2*k2 + 2*k3 + k4)/6
      t = step*DT
      energy = sum(u**2)
      if (energy > energy_prev*(1 + 1.e3_R_P*epsilon(1._R_P))) then
         fails = fails + 1
         print '(A,A,A,I0)', 'ENERGY GROWS ', trim(NAMES(kk)), ' at step ', step
         exit
      endif
      energy_prev = energy
      if (abs(sum(u) - mass0) > 1.e3_R_P*epsilon(1._R_P)*sum(abs(u0))) then
         fails = fails + 1
         print '(A,A,A,I0)', 'MASS NOT KEPT ', trim(NAMES(kk)), ' at step ', step
         exit
      endif
      if (report .and. mod(step, EVERY) == 0) then
         call averages(t, ex)
         write(uf, '(*(1X,ES24.16E3))') real(t, R8P), real(u, R8P), real(ex, R8P)
      endif
   enddo
   if (report) close(uf)
   call averages(t, ex)
   err(kk) = maxval(abs(u - ex))
   energy_end(kk) = sum(u**2)/energy0
   print '(A,A,A,F6.3,A,ES10.3,A,F7.4)', 'gideon_vv_advection: ', NAMES(kk), ' t = ', t, ', max error ', err(kk), &
                                         ', energy kept ', energy_end(kk)
enddo
if (energy_end(3) > 0.5_R_P) then
   fails = fails + 1
   print '(A)', 'UPWIND DOES NOT DISSIPATE'
endif
if (err(2) >= err(1)) then
   fails = fails + 1
   print '(A)', 'FOURTH ORDER NOT BETTER THAN SECOND'
endif
if (fails > 0) error stop 'gideon_vv_advection: failed'
print '(A)', 'gideon_vv_advection: ok'
contains
   subroutine rhs(kk, v, dvdt)
   !< -du/dx averaged over every cell, periodic, by kernel kk.
   integer(I4P), intent(in)  :: kk
   real(R_P),    intent(in)  :: v(N)
   real(R_P),    intent(out) :: dvdt(N)
   real(R_P)                 :: w(1-G:N+G), d
   integer(I4P)              :: j
   w(1:N) = v
   w(1-G:0) = v(N-G+1:N)
   w(N+1:N+G) = v(1:G)
   do j=1, N
      select case(kk)
      case(1) ; call gideon_fv_d1_c2(w(j-1:j+1), H, d)
      case(2) ; call gideon_fv_d1_c4(w(j-2:j+2), H, d)
      case(3) ; call gideon_fv_d1_ul2(w(j-2:j), H, d)
      endselect
      dvdt(j) = -d
   enddo
   endsubroutine rhs

   pure function u_exact(x) result(v)
   !< The packet, periodic.
   real(R_P), intent(in) :: x
   real(R_P)             :: v, y
   y = x - floor(x)
   v = exp(-((y - X0)/SIGMA)**2)*sin(16*PI*y)
   endfunction u_exact

   subroutine averages(tt, v)
   !< Cell averages of u0(x - tt), 5-point Gauss-Legendre per cell.
   real(R_P), intent(in)  :: tt
   real(R_P), intent(out) :: v(N)
   real(R_P), parameter   :: GX(5) = [-0.906179845938663992797626878299392965_R_P, &
                                      -0.538469310105683091036314420700208805_R_P, 0._R_P, &
                                       0.538469310105683091036314420700208805_R_P, &
                                       0.906179845938663992797626878299392965_R_P]
   real(R_P), parameter   :: GW(5) = [0.236926885056189087514264040719917363_R_P, &
                                      0.478628670499366468041291514835638192_R_P, &
                                      0.568888888888888888888888888888888889_R_P, &
                                      0.478628670499366468041291514835638192_R_P, &
                                      0.236926885056189087514264040719917363_R_P]
   integer(I4P)           :: j, m
   do j=1, N
      v(j) = 0._R_P
      do m=1, 5
         v(j) = v(j) + GW(m)*u_exact(xc(j) + 0.5_R_P*H*GX(m) - tt)
      enddo
      v(j) = 0.5_R_P*v(j)
   enddo
   endsubroutine averages
endprogram gideon_vv_advection
