!< GIDEON test (#2, L3 and L5): the line solvers against a dense solve, on the host and on the device.
program gideon_test_solvers
!< GIDEON test (#2, L3 and L5): the line solvers against a dense solve, on the host and on the device.
!<
!< Random strictly diagonally dominant bands (off-diagonals in [-1, 1], the diagonal of either sign, its dominance
!< margin in [0.25, 1.25]) and random right sides, on lines of 3 (tridiagonal) or 5 (pentadiagonal) to 64 unknowns,
!< bounded and periodic, constant band and variable coefficients. Frozen bounds, written before the solvers were run:
!< - forward error: max|x - x_ref| <= 64 eps cond max|x_ref|, x_ref a dense solve with partial pivoting in RD (real128
!<   where the compiler has it, else real64), cond <= ||A|| / min_i (|b_i| - sum_j |a_ij|) (Varah's bound on the
!<   inverse of a strictly diagonally dominant matrix);
!< - residual, computed in RD: max|A x - d| <= 64 eps (||A|| max|x| + max|d|) (backward stability);
!< - the coefficients outside a bounded line are never read: the variable variant gets NaN there;
!< - the constant-band variant equals the variable one on the same band, bitwise (one template, the same operations);
!< - the masked variant (#2 4.1: the constant band, the identity on the rows of a random mask, ~30 %) equals the
!<   variable one on that band, bitwise, is within the bounds above, and returns on every masked row exactly the value
!<   given there;
!< - device = host: NL lines solved in one `!$acc parallel loop`, one line per thread, against the host, within the
!<   forward bound (the device contracts multiply-adds differently); bitwise in a build without OpenACC.
use, intrinsic :: iso_fortran_env, only : real64, real128
use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
use penf, only : I4P, R_P
use gideon_line_solvers
implicit none
integer,      parameter :: RD = merge(real128, real64, real128 > 0) !< Kind of the dense reference.
integer(I4P), parameter :: NMAX = 64_I4P      !< Longest line.
integer(I4P), parameter :: NSIZES = 6_I4P     !< Line lengths per band.
integer(I4P), parameter :: NTRIAL = 10_I4P    !< Systems per length.
integer(I4P), parameter :: NL = 4096_I4P      !< Lines of the device check.
integer(I4P), parameter :: NDEV = 64_I4P      !< Their length.
integer(I4P), parameter :: SIZES_TRI(NSIZES) = [3, 4, 5, 8, 17, 64]
integer(I4P), parameter :: SIZES_PENTA(NSIZES) = [5, 6, 7, 8, 17, 64]
real(R_P)    :: bd(-2:2,NMAX), d(NMAX), x(NMAX), xc(NMAX), w1(NMAX), w2(NMAX), w3(NMAX), w4(NMAX), nan
real(R_P)    :: bdm(-2:2,NMAX), um(NMAX)
logical      :: mask(NMAX)
integer(I4P) :: fails, checks, half, is, n, trial, periodic_case, i
logical      :: periodic
integer, allocatable :: seed(:)
integer      :: ns

call random_seed(size=ns)
allocate(seed(ns)) ; seed = 20261118 ; call random_seed(put=seed)
nan = ieee_value(1._R_P, ieee_quiet_nan)
fails = 0 ; checks = 0
do half=1, 2
   do is=1, NSIZES
      do trial=1, NTRIAL
         do periodic_case=0, 1
            periodic = periodic_case == 1
            n = merge(SIZES_TRI(is), SIZES_PENTA(is), half == 1)
            call random_band(half, n, bd)
            call random_number(d(1:n)) ; d(1:n) = 2*d(1:n) - 1
            ! variable coefficients, the corners outside a bounded line NaN
            call solve_var(half, n, periodic, bd, d, x)
            call check_dense(half, n, periodic, bd, d, x)
            ! constant band: the first row's band on every row, against the variable variant on that band
            call constant_band(half, n, bd)
            call solve_var(half, n, periodic, bd, d, x)
            call solve_const(half, n, periodic, bd(:,1), d, xc)
            call check_dense(half, n, periodic, bd, d, xc)
            checks = checks + 1
            if (any(xc(1:n) /= x(1:n))) then
               fails = fails + 1
               print '(A,I0,A,I0,A,L1)', 'CONST /= VAR: half-band ', half, ', n ', n, ', periodic ', periodic
            endif
            ! masked: the constant band, the identity on the masked rows
            call random_number(um(1:n)) ; mask(1:n) = um(1:n) < 0.3_R_P
            bdm = bd
            do i=1, n
               if (mask(i)) then
                  bdm(:,i) = 0._R_P ; bdm(0,i) = 1._R_P
               endif
            enddo
            call solve_var(half, n, periodic, bdm, d, x)
            call solve_masked(half, n, periodic, bd(:,1), mask, d, xc)
            call check_dense(half, n, periodic, bdm, d, xc)
            checks = checks + 2
            if (any(xc(1:n) /= x(1:n))) then
               fails = fails + 1
               print '(A,I0,A,I0,A,L1)', 'MASKED /= VAR: half-band ', half, ', n ', n, ', periodic ', periodic
            endif
            if (any(mask(1:n) .and. xc(1:n) /= d(1:n))) then
               fails = fails + 1
               print '(A,I0,A,I0,A,L1)', 'MASKED ROW NOT THE VALUE GIVEN: half-band ', half, ', n ', n, &
                                         ', periodic ', periodic
            endif
         enddo
      enddo
   enddo
enddo
call check_device
print '(A,I0,A,I0)', 'gideon_test_solvers: ', checks, ' checks, failures ', fails
if (fails > 0) error stop 'gideon_test_solvers: failed'
print '(A)', 'gideon_test_solvers: ok'
contains
   subroutine random_band(half, n, bd)
   !< A strictly diagonally dominant band: bd(o,i) multiplies x(i+o).
   integer(I4P), intent(in)  :: half, n
   real(R_P),    intent(out) :: bd(-2:2,NMAX)
   real(R_P)                 :: u(-2:2), margin
   integer(I4P)              :: i
   bd = 0._R_P
   do i=1, n
      call random_number(u)
      bd(-half:half,i) = 2*u(-half:half) - 1
      bd(0,i) = 0._R_P
      call random_number(margin)
      bd(0,i) = sum(abs(bd(:,i))) + 0.25_R_P + margin
      if (u(0) < 0.5_R_P) bd(0,i) = -bd(0,i)
   enddo
   endsubroutine random_band

   subroutine constant_band(half, n, bd)
   !< The first row's band on every row.
   integer(I4P), intent(in)    :: half, n
   real(R_P),    intent(inout) :: bd(-2:2,NMAX)
   integer(I4P)              :: i
   do i=2, n
      bd(-half:half,i) = bd(-half:half,1)
   enddo
   endsubroutine constant_band

   subroutine solve_var(half, n, periodic, bd, d, x)
   !< The variable-coefficient solvers; on a bounded line the coefficients outside it are NaN (never read).
   integer(I4P), intent(in)  :: half, n
   logical,      intent(in)  :: periodic
   real(R_P),    intent(in)  :: bd(-2:2,NMAX), d(NMAX)
   real(R_P),    intent(out) :: x(NMAX)
   real(R_P)                 :: b(-2:2,NMAX)
   integer(I4P)              :: i, o
   b = bd
   if (.not. periodic) then
      do i=1, n
         do o=-half, half
            if (i + o < 1 .or. i + o > n) b(o,i) = nan
         enddo
      enddo
   endif
   x = d
   if (half == 1) then
      if (periodic) then
         call gideon_tri_periodic_var(n, b(-1,1:n), b(0,1:n), b(1,1:n), x, w1, w2)
      else
         call gideon_tri_bounded_var(n, b(-1,1:n), b(0,1:n), b(1,1:n), x, w1)
      endif
   else
      if (periodic) then
         call gideon_penta_periodic_var(n, b(-2,1:n), b(-1,1:n), b(0,1:n), b(1,1:n), b(2,1:n), x, w1, w2, w3, w4)
      else
         call gideon_penta_bounded_var(n, b(-2,1:n), b(-1,1:n), b(0,1:n), b(1,1:n), b(2,1:n), x, w1, w2)
      endif
   endif
   endsubroutine solve_var

   subroutine solve_const(half, n, periodic, band, d, x)
   !< The constant-band solvers.
   integer(I4P), intent(in)  :: half, n
   logical,      intent(in)  :: periodic
   real(R_P),    intent(in)  :: band(-2:2), d(NMAX)
   real(R_P),    intent(out) :: x(NMAX)
   x = d
   if (half == 1) then
      if (periodic) then
         call gideon_tri_periodic_const(n, band(-1), band(0), band(1), x, w1, w2)
      else
         call gideon_tri_bounded_const(n, band(-1), band(0), band(1), x, w1)
      endif
   else
      if (periodic) then
         call gideon_penta_periodic_const(n, band(-2), band(-1), band(0), band(1), band(2), x, w1, w2, w3, w4)
      else
         call gideon_penta_bounded_const(n, band(-2), band(-1), band(0), band(1), band(2), x, w1, w2)
      endif
   endif
   endsubroutine solve_const

   subroutine solve_masked(half, n, periodic, band, mask, d, x)
   !< The masked solvers: the constant band, the identity on the masked rows.
   integer(I4P), intent(in)  :: half, n
   logical,      intent(in)  :: periodic, mask(NMAX)
   real(R_P),    intent(in)  :: band(-2:2), d(NMAX)
   real(R_P),    intent(out) :: x(NMAX)
   x = d
   if (half == 1) then
      if (periodic) then
         call gideon_tri_periodic_masked(n, band(-1), band(0), band(1), mask, x, w1, w2)
      else
         call gideon_tri_bounded_masked(n, band(-1), band(0), band(1), mask, x, w1)
      endif
   else
      if (periodic) then
         call gideon_penta_periodic_masked(n, band(-2), band(-1), band(0), band(1), band(2), mask, x, w1, w2, w3, w4)
      else
         call gideon_penta_bounded_masked(n, band(-2), band(-1), band(0), band(1), band(2), mask, x, w1, w2)
      endif
   endif
   endsubroutine solve_masked

   subroutine dense(half, n, periodic, bd, a)
   !< The dense matrix of the band (indices modulo n on a periodic line), in RD.
   integer(I4P), intent(in)  :: half, n
   logical,      intent(in)  :: periodic
   real(R_P),    intent(in)  :: bd(-2:2,NMAX)
   real(RD),     intent(out) :: a(NMAX,NMAX)
   integer(I4P)              :: i, o, j
   a = 0._RD
   do i=1, n
      do o=-half, half
         j = i + o
         if (periodic) then
            j = modulo(j - 1, n) + 1
         elseif (j < 1 .or. j > n) then
            cycle
         endif
         a(i,j) = a(i,j) + real(bd(o,i), RD)
      enddo
   enddo
   endsubroutine dense

   subroutine check_dense(half, n, periodic, bd, d, x)
   !< Forward error against the dense solve, and the residual, within the frozen bounds.
   integer(I4P), intent(in) :: half, n
   logical,      intent(in) :: periodic
   real(R_P),    intent(in) :: bd(-2:2,NMAX), d(NMAX), x(NMAX)
   real(RD)                 :: a(NMAX,NMAX), lu(NMAX,NMAX), xr(NMAX), r(NMAX), anorm, dom, cond, t
   integer(I4P)             :: i, j, k, p
   call dense(half, n, periodic, bd, a)
   anorm = 0._RD ; dom = huge(1._RD)
   do i=1, n
      anorm = max(anorm, sum(abs(a(i,1:n))))
      dom = min(dom, 2*abs(a(i,i)) - sum(abs(a(i,1:n))))
   enddo
   cond = anorm/dom
   ! Gaussian elimination with partial pivoting
   lu = a ; xr(1:n) = real(d(1:n), RD)
   do k=1, n - 1
      p = k - 1 + maxloc(abs(lu(k:n,k)), dim=1)
      if (p /= k) then
         r(1:n) = lu(k,1:n) ; lu(k,1:n) = lu(p,1:n) ; lu(p,1:n) = r(1:n)
         t = xr(k) ; xr(k) = xr(p) ; xr(p) = t
      endif
      do i=k + 1, n
         t = lu(i,k)/lu(k,k)
         lu(i,k:n) = lu(i,k:n) - t*lu(k,k:n)
         xr(i) = xr(i) - t*xr(k)
      enddo
   enddo
   do i=n, 1, -1
      xr(i) = (xr(i) - sum(lu(i,i+1:n)*xr(i+1:n)))/lu(i,i)
   enddo
   checks = checks + 2
   if (maxval(abs(real(x(1:n), RD) - xr(1:n))) > 64*epsilon(1._R_P)*cond*maxval(abs(xr(1:n)))) then
      fails = fails + 1
      print '(A,I0,A,I0,A,L1,A,ES10.2,A,ES10.2)', 'FORWARD ERROR half-band ', half, ', n ', n, ', periodic ', &
         periodic, ': ', real(maxval(abs(real(x(1:n), RD) - xr(1:n)))), ', cond ', real(cond)
   endif
   do i=1, n
      r(i) = real(d(i), RD)
      do j=1, n
         r(i) = r(i) - a(i,j)*real(x(j), RD)
      enddo
   enddo
   if (maxval(abs(r(1:n))) > 64*epsilon(1._R_P)*(anorm*maxval(abs(real(x(1:n), RD))) + maxval(abs(real(d(1:n), RD))))) &
      then
      fails = fails + 1
      print '(A,I0,A,I0,A,L1,A,ES10.2)', 'RESIDUAL half-band ', half, ', n ', n, ', periodic ', periodic, ': ', &
         real(maxval(abs(r(1:n))))
   endif
   endsubroutine check_dense

   subroutine check_device
   !< NL periodic and bounded lines, tridiagonal and pentadiagonal, variable coefficients: device against host.
   real(R_P), allocatable :: bdg(:,:,:), dd(:,:), xd(:,:), xh(:,:), s1(:,:), s2(:,:), s3(:,:), s4(:,:)
   real(RD)               :: a(NMAX,NMAX), anorm, dom, worst
   integer(I4P)           :: l, i, hb, pc, nbit
   logical                :: per
   allocate(bdg(NDEV,NL,-2:2), dd(NDEV,NL), xd(NDEV,NL), xh(NDEV,NL))
   allocate(s1(NDEV,NL), s2(NDEV,NL), s3(NDEV,NL), s4(NDEV,NL))
   do hb=1, 2
      do pc=0, 1
         per = pc == 1
         do l=1, NL
            call random_band(hb, NDEV, bd)
            do i=-2, 2
               bdg(:,l,i) = bd(i,1:NDEV)
            enddo
            call random_number(dd(:,l)) ; dd(:,l) = 2*dd(:,l) - 1
         enddo
         xd = dd ; xh = dd
         !$acc parallel loop copyin(bdg) copy(xd) create(s1, s2, s3, s4)
         do l=1, NL
            if (hb == 1 .and. per) then
               call gideon_tri_periodic_var(NDEV, bdg(:,l,-1), bdg(:,l,0), bdg(:,l,1), xd(:,l), s1(:,l), s2(:,l))
            elseif (hb == 1) then
               call gideon_tri_bounded_var(NDEV, bdg(:,l,-1), bdg(:,l,0), bdg(:,l,1), xd(:,l), s1(:,l))
            elseif (per) then
               call gideon_penta_periodic_var(NDEV, bdg(:,l,-2), bdg(:,l,-1), bdg(:,l,0), bdg(:,l,1), &
                                              bdg(:,l,2), xd(:,l), s1(:,l), s2(:,l), s3(:,l), s4(:,l))
            else
               call gideon_penta_bounded_var(NDEV, bdg(:,l,-2), bdg(:,l,-1), bdg(:,l,0), bdg(:,l,1), &
                                             bdg(:,l,2), xd(:,l), s1(:,l), s2(:,l))
            endif
         enddo
         nbit = 0 ; worst = 0._RD
         do l=1, NL
            do i=-2, 2
               bd(i,1:NDEV) = bdg(:,l,i)
            enddo
            call solve_var(hb, NDEV, per, bd, dd(:,l), x)
            xh(:,l) = x(1:NDEV)
            call dense(hb, NDEV, per, bd, a)
            anorm = 0._RD ; dom = huge(1._RD)
            do i=1, NDEV
               anorm = max(anorm, sum(abs(a(i,1:NDEV))))
               dom = min(dom, 2*abs(a(i,i)) - sum(abs(a(i,1:NDEV))))
            enddo
            checks = checks + 1
            if (all(xd(:,l) == xh(:,l))) nbit = nbit + 1
            worst = max(worst, maxval(abs(real(xd(:,l) - xh(:,l), RD)))/(anorm/dom*maxval(abs(real(xh(:,l), RD)))))
            if (maxval(abs(xd(:,l) - xh(:,l))) > 64*epsilon(1._R_P)*anorm/dom*maxval(abs(xh(:,l)))) then
               fails = fails + 1
               print '(A,I0,A,L1,A,I0)', 'DEVICE half-band ', hb, ', periodic ', per, ', line ', l
            endif
         enddo
#if defined _OPENACC
         print '(A,I0,A,L1,A,I0,A,I0,A,ES9.2)', 'gideon_test_solvers: device, half-band ', hb, ', periodic ', per, &
            ': lines bitwise ', nbit, ' of ', NL, ', worst error / (cond max|x|) ', real(worst)
#else
         checks = checks + 1
         if (nbit /= NL) then
            fails = fails + 1
            print '(A)', 'HOST LOOP NOT BITWISE'
         endif
#endif
      enddo
   enddo
   endsubroutine check_device
endprogram gideon_test_solvers
