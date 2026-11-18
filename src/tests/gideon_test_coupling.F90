!< GIDEON test (#2, 4.1 and 4.2): a compact-WENO coupling on one line, the compact rows where smooth, given values at a step.
program gideon_test_coupling
!< GIDEON test (#2, 4.1 and 4.2): a compact-WENO coupling on one line, the compact rows where smooth, given values at a step.
!<
!< u(x) = sin(2 pi x) + 0.5 cos(6 pi x + 0.3) + H(x), H = 1 on [1/2, 1), 0 on [0, 1/2): a smooth wave and a step, two
!< jumps (x = 1/2 up, x = 0 down) on faces of a periodic line of N cells of [0, 1). The face values: the upwind-biased
!< compact operator (gideon_fv_face_pbl5, pbl7, its band FV_FACE_PBL<o>_LEFT) on every face farther than 3 cells from a
!< jump; on the others (flagged: where a sensor would switch) the value given, here the first-order upwind one, the cell
!< average on the left (a stand-in for the nonlinear WENO/TENO value: GIDEON provides the operator and the solve, the
!< switch and the nonlinear value are the caller's, D7). One solve: gideon_tri_periodic_masked.
!<
!< Frozen bounds, written before the test was run:
!< - every flagged face returns the value given, bitwise;
!< - conservation: the flux differences f_j - f_(j-1), summed in RD (real128 where the compiler has it), vanish within
!<   64 eps(RD) sum|f|. A flux form conserves whatever its face values are: the check guards that every face has one
!<   value, the one both of its cells use. (The first version summed the update in R_P against 64 eps: a sum of N terms
!<   rounds up to N eps, and gfortran's debug build exceeded it at N = 131072.)
!< - away from the step (faces at least 0.2 from both jumps) the error to the exact face value has the design order:
!<   N = 2^(k+5), a plateau of three resolved sizes (error > 10 eps F sum|w| (max|u| + 2 max|P|/h), F the band's
!<   amplification, P the primitive the cell averages are differences of, over h) whose two observed orders differ by at
!<   most 0.25, both >= order - 0.25; measured for both operators in R16P (required). (The first version of the floor
!<   omitted the data's own round-off, 2 max|P|/h, which grows like N: the finest sizes counted as resolved and showed
!<   order -1. test_order has carried that term from the start.)
!< - the device: NLD coupled lines (random data and masks), one per thread, against the host within 64 eps F sum|w|
!<   max|u| (bitwise without OpenACC).
use, intrinsic :: iso_fortran_env, only : real64, real128
use penf, only : I4P, R_P
use gideon
implicit none
integer,      parameter :: RD = merge(real128, real64, real128 > 0) !< Kind of the conservation sum.
real(R_P),    parameter :: PI = 3.14159265358979323846264338327950288_R_P
integer(I4P), parameter :: KMAX = 12_I4P
integer(I4P), parameter :: FLAG_CELLS = 3_I4P
real(R_P),    parameter :: PMAX = 1.5_R_P !< A bound of |prim| on [0, 1].
integer(I4P)            :: fails, checks, iop

fails = 0 ; checks = 0
do iop=5, 7, 2
   call check_line(iop)
enddo
call check_device
print '(A,I0,A,I0)', 'gideon_test_coupling: 2 operators, ', checks, ' checks, failures ', fails
if (fails > 0) error stop 'gideon_test_coupling: failed'
print '(A)', 'gideon_test_coupling: ok'
contains
   pure function prim(x) result(v)
   !< A primitive of the smooth part and of the step.
   real(R_P), intent(in) :: x
   real(R_P)             :: v
   v = -cos(2*PI*x)/(2*PI) + 0.5_R_P*sin(6*PI*x + 0.3_R_P)/(6*PI) + max(x - 0.5_R_P, 0._R_P)
   endfunction prim

   pure function exact(x) result(v)
   !< u away from the jumps.
   real(R_P), intent(in) :: x
   real(R_P)             :: v
   v = sin(2*PI*x) + 0.5_R_P*cos(6*PI*x + 0.3_R_P) + merge(1._R_P, 0._R_P, x >= 0.5_R_P)
   endfunction exact

   subroutine band_of(order, am, ap, sw, ff)
   !< The band of pbl<order>, sum|w| of its right side, the amplification F of the round-off through the solve.
   integer(I4P), intent(in)  :: order
   real(R_P),    intent(out) :: am, ap, sw, ff
   real(R_P)                 :: s
   if (order == 5) then
      am = FV_FACE_PBL5_LEFT(-1) ; ap = FV_FACE_PBL5_LEFT(1) ; sw = sum(abs(FV_FACE_PBR5_RIGHT))
   else
      am = FV_FACE_PBL7_LEFT(-1) ; ap = FV_FACE_PBL7_LEFT(1) ; sw = sum(abs(FV_FACE_PBR7_RIGHT))
   endif
   s = 1._R_P - abs(am) - abs(ap) ! the band's diagonal dominance bounds its smallest symbol from below
   ff = (1._R_P + (1._R_P + abs(am) + abs(ap))/s)/s
   endsubroutine band_of

   pure subroutine rhs(order, n, u, j, r)
   !< The right side of pbl<order> at face j (between cells j and j+1) of a periodic line.
   integer(I4P), intent(in)  :: order, n, j
   real(R_P),    intent(in)  :: u(n)
   real(R_P),    intent(out) :: r
   real(R_P)                 :: q(-2:2)
   integer(I4P)              :: o
   do o=-2, 2
      q(o) = u(modulo(j + o - 1, n) + 1)
   enddo
   if (order == 5) then
      call gideon_fv_face_pbl5(q(-1:1), r)
   else
      call gideon_fv_face_pbl7(q(-2:2), r)
   endif
   endsubroutine rhs

   subroutine coupled(order, n, u, mask, f, w, z)
   !< The face values of the coupled line: compact rows where not masked, the upwind average where masked.
   integer(I4P), intent(in)  :: order, n
   real(R_P),    intent(in)  :: u(n)
   logical,      intent(in)  :: mask(n)
   real(R_P),    intent(out) :: f(n), w(n), z(n)
   real(R_P)                 :: am, ap, sw, ff
   integer(I4P)              :: j
   call band_of(order, am, ap, sw, ff)
   do j=1, n
      if (mask(j)) then
         f(j) = u(j)
      else
         call rhs(order, n, u, j, f(j))
      endif
   enddo
   call gideon_tri_periodic_masked(n, am, 1._R_P, ap, mask, f, w, z)
   endsubroutine coupled

   subroutine check_line(order)
   !< Values given where flagged, conservation, the order away from the step.
   integer(I4P), intent(in) :: order
   real(R_P), allocatable   :: u(:), f(:), w(:), z(:)
   real(RD)                 :: cons
   logical, allocatable     :: mask(:)
   real(R_P)    :: h, x, d, e(KMAX), noise(KMAX), p(KMAX), am, ap, sw, ff
   logical      :: ok(KMAX)
   integer(I4P) :: k, n, j, kbest
   call band_of(order, am, ap, sw, ff)
   ok = .false. ; e = 0._R_P
   do k=1, KMAX
      n = 2**(k + 5)
      h = 1._R_P/n
      allocate(u(n), f(n), w(n), z(n), mask(n))
      do j=1, n ! cell j = [(j-1) h, j h], face j at x = j h
         u(j) = (prim(j*h) - prim((j - 1)*h))/h
         d = min(abs(j*h - 0.5_R_P), abs(j*h - 1._R_P), j*h)
         mask(j) = d <= FLAG_CELLS*h*(1 + 1.e-3_R_P)
      enddo
      call coupled(order, n, u, mask, f, w, z)
      ! the values given
      checks = checks + 1
      if (any(mask .and. f /= u)) then
         fails = fails + 1
         print '(A,I0,A,I0)', 'FLAGGED FACE NOT THE VALUE GIVEN: pbl', order, ', N = ', n
      endif
      ! conservation: the flux differences, summed in RD
      cons = 0._RD
      do j=1, n
         cons = cons + (real(f(j), RD) - real(f(modulo(j - 2, n) + 1), RD))
      enddo
      checks = checks + 1
      if (abs(cons) > 64*epsilon(1._RD)*sum(abs(real(f, RD)))) then
         fails = fails + 1
         print '(A,I0,A,I0,A,ES10.2)', 'NOT CONSERVATIVE: pbl', order, ', N = ', n, ': ', real(abs(cons))
      endif
      ! the error away from the step
      do j=1, n
         x = j*h
         d = min(abs(x - 0.5_R_P), abs(x - 1._R_P), x)
         if (d >= 0.2_R_P) e(k) = max(e(k), abs(f(j) - exact(x)))
      enddo
      noise(k) = epsilon(1._R_P)*ff*sw*(maxval(abs(u)) + 2*PMAX/h)
      ok(k) = e(k) > 10*noise(k)
      deallocate(u, f, w, z, mask)
   enddo
   p = 0._R_P
   do k=1, KMAX - 1
      if (ok(k) .and. ok(k+1)) p(k) = log(e(k)/e(k+1))/log(2._R_P)
   enddo
   kbest = 0
   do k=1, KMAX - 2
      if (ok(k) .and. ok(k+1) .and. ok(k+2)) then
         if (abs(p(k) - p(k+1)) <= 0.25_R_P) kbest = k
      endif
   enddo
   checks = checks + 1
   if (kbest == 0) then
      print '(A,I0,A)', 'gideon_test_coupling: pbl', order, ' away from the step: unresolved in this precision'
      if (precision(1._R_P) >= 30) then
         fails = fails + 1
         print '(A,I0)', 'UNRESOLVED IN QUADRUPLE PRECISION: pbl', order
      endif
   elseif (min(p(kbest), p(kbest+1)) < order - 0.25_R_P) then
      fails = fails + 1
      print '(A,I0,A,2F7.2)', 'ORDER AWAY FROM THE STEP: pbl', order, ': observed ', p(kbest), p(kbest+1)
   else
      print '(A,I0,A,2F7.2,A,I0)', 'gideon_test_coupling: pbl', order, ' away from the step, observed orders ', &
         p(kbest), p(kbest+1), ' at N = ', 2**(kbest + 5)
   endif
   endsubroutine check_line

   subroutine check_device
   !< NLD coupled lines of random data and masks on the device against the host.
   integer(I4P), parameter :: NLD = 1024_I4P, NN = 128_I4P
   real(R_P), allocatable  :: u(:,:), fd(:,:), fh(:,:), w(:,:), z(:,:), um(:,:)
   logical, allocatable    :: mask(:,:)
   real(R_P)               :: am, ap, sw, ff, q(-2:2)
   integer(I4P)            :: l, j, o, order, nbit
   integer, allocatable    :: seed(:)
   integer                 :: ns
   call random_seed(size=ns)
   allocate(seed(ns)) ; seed = 41 ; call random_seed(put=seed)
   allocate(u(NN,NLD), fd(NN,NLD), fh(NN,NLD), w(NN,NLD), z(NN,NLD), um(NN,NLD), mask(NN,NLD))
   call random_number(u) ; u = 2*u - 1
   call random_number(um) ; mask = um < 0.2_R_P
   do order=5, 7, 2
      call band_of(order, am, ap, sw, ff)
      !$acc parallel loop copyin(u, mask) copyout(fd) create(w, z) private(q)
      do l=1, NLD
         !$acc loop seq
         do j=1, NN
            if (mask(j,l)) then
               fd(j,l) = u(j,l)
            else
               !$acc loop seq
               do o=-2, 2
                  q(o) = u(modulo(j + o - 1, NN) + 1,l)
               enddo
               if (order == 5) then
                  call gideon_fv_face_pbl5(q(-1:1), fd(j,l))
               else
                  call gideon_fv_face_pbl7(q(-2:2), fd(j,l))
               endif
            endif
         enddo
         call gideon_tri_periodic_masked(NN, am, 1._R_P, ap, mask(:,l), fd(:,l), w(:,l), z(:,l))
      enddo
      nbit = 0
      do l=1, NLD
         call coupled(order, NN, u(:,l), mask(:,l), fh(:,l), w(:,l), z(:,l))
         checks = checks + 1
         if (all(fd(:,l) == fh(:,l))) nbit = nbit + 1
         if (maxval(abs(fd(:,l) - fh(:,l))) > 64*epsilon(1._R_P)*ff*sw*maxval(abs(u(:,l)))) then
            fails = fails + 1
            if (fails <= 10) print '(A,I0,A,I0)', 'DEVICE COUPLED LINE pbl', order, ', line ', l
         endif
      enddo
#if defined _OPENACC
      print '(A,I0,A,I0,A,I0)', 'gideon_test_coupling: device, pbl', order, ': lines bitwise ', nbit, ' of ', NLD
#else
      checks = checks + 1
      if (nbit /= NLD) then
         fails = fails + 1
         print '(A)', 'HOST LOOP NOT BITWISE'
      endif
#endif
   enddo
   endsubroutine check_device
endprogram gideon_test_coupling
