!< GIDEON test (#2, L1, L2, L4, L5): the compact operators, row by row, on periodic lines, spectrally, on the device.
program gideon_test_lines
!< GIDEON test (#2, L1, L2, L4, L5): the compact operators, row by row, on periodic lines, spectrally, on the device.
!<
!< Frozen bounds, written before the operators were run:
!< - L1 (rows): every compiled row is exact on the monomials up to its degree K = order + n - 1: |sum_o alpha_o
!<   T(x^j)(o) - r| <= tol, and on x^(K+1) its residual is the generator's exact one (lerror) within tol; tol = 64 eps
!<   (sum|w| max|q| + sum|alpha| max|T|), h = 1.
!< - L2 (order on a periodic line): g(x) = sin(x + 0.4) + 0.5 cos(2x - 0.2) + 0.3 sin(3x + 1.1) on N = 8 .. 2^16
!<   points of [0, 2 pi), the right side by the kernel, the band solved by gideon_line_solvers; the error max|u -
!<   exact| has a plateau (three resolved steps, two observed orders within 0.25) with both orders >= design - 0.25.
!<   Resolved: error > 10 times the round-off floor eps (sum|w| max|q| / h^n) F + eps max|exact|, and 3 width h <= 1.
!<   Unresolved operators are listed in R4P and R8P; in R16P every operator must be measured.
!< - L4 (symbol): on N = 64 points, the solve of the Fourier modes k = 1 .. N/2 - 1 equals the right side's symbol (the
!<   kernel on the mode) over the band's symbol, within 64 eps sum|w| / h^n F; centred d1 and d2: the modified
!<   wavenumber is real (no dissipation) within the same tolerance.
!< - L5 (device): the right sides of NP random windows (64 eps sum|w| max|q| / h^n), and NLD periodic lines solved one
!<   per thread (that times F), on the device against the host (bitwise without OpenACC).
!< Compact filters (alpha_f at run time, band [alpha_f, 1, alpha_f]): L1 at alpha_f = 0, 0.3, 0.45 (the leading
!< residual lerror + alpha_f lerror_alpha), L2, L4, L5 at alpha_f = AF = 0.45; L4 also: the transfer function is real
!< and in [0, 1] within the tolerance, and the odd-even mode (-1)^i is removed within it.
!< F = (1 + sum|alpha| / s) / s, s = min over xi of |band symbol|: the round-off of the right side (relative eps) and
!< the solve's backward error (relative eps to the band, on a solution of size sum|w| max|q| / h^n / s), both through
!< the inverse band, whose norm is 1/s. (The first version of this test, before any run, had F = 1/s: it omitted the
!< solve's backward error, and the L4 check failed by up to 1.3 times on the fourth-order operators.)
use penf, only : I4P, R_P
use gideon_test_lines_registry
implicit none
real(R_P),    parameter :: PI = 3.14159265358979323846264338327950288_R_P
real(R_P),    parameter :: X0 = 0.3_R_P       !< The first unknown of a periodic line.
integer(I4P), parameter :: KMAX = 14_I4P      !< N = 2^(k+2), k = 1..KMAX.
integer(I4P), parameter :: NSYM = 64_I4P      !< Line of the symbol check.
real(R_P),    parameter :: AF = 0.45_R_P      !< alpha_f of the filters in L2, L4, L5.
integer(I4P)            :: fails, checks

fails = 0 ; checks = 0
call check_rows
call check_orders
call check_symbols
call check_device
print '(A,I0,A,I0,A,I0)', 'gideon_test_lines: ', NL, ' compact operators, ', checks, ' checks, failures ', fails
if (fails > 0) error stop 'gideon_test_lines: failed'
print '(A)', 'gideon_test_lines: ok'
contains
   pure function gp(k, x) result(v)
   !< k-th derivative (k >= -1, -1 a primitive) of the periodic g.
   integer(I4P), intent(in) :: k
   real(R_P),    intent(in) :: x
   real(R_P)                :: v
   v = sin(x + 0.4_R_P + k*PI/2) + 0.5_R_P*2._R_P**k*cos(2*x - 0.2_R_P + k*PI/2) + &
       0.3_R_P*3._R_P**k*sin(3*x + 1.1_R_P + k*PI/2)
   endfunction gp

   pure function dmono(n, j, x) result(v)
   !< n-th derivative of x^j at x.
   integer(I4P), intent(in) :: n, j
   real(R_P),    intent(in) :: x
   real(R_P)                :: v
   integer(I4P)             :: i
   v = 0._R_P
   if (j < n) return
   v = 1._R_P
   do i=j - n + 1, j
      v = v*i
   enddo
   if (j > n) v = v*x**(j - n)
   endfunction dmono

   pure function band(id, af, o) result(a)
   !< The band of operator id at offset o (a filter's: [af, 1, af]).
   integer(I4P), intent(in) :: id, o
   real(R_P),    intent(in) :: af
   real(R_P)                :: a
   if (lalpha_free(id)) then
      a = merge(1._R_P, merge(af, 0._R_P, abs(o) == 1), o == 0)
   else
      a = lleft(o,id)
   endif
   endfunction band

   pure function sumabs(id, af) result(s)
   !< sum|w| of the right side (a filter's at af).
   integer(I4P), intent(in) :: id
   real(R_P),    intent(in) :: af
   real(R_P)                :: s
   s = lsumabs(id) + abs(af)*lsumabs_alpha(id)
   endfunction sumabs

   pure function floor_factor(id, af) result(f)
   !< F = (1 + sum|alpha|/s)/s, s the band's smallest symbol: the amplification of the round-off through the solve.
   integer(I4P), intent(in) :: id
   real(R_P),    intent(in) :: af
   real(R_P)                :: f, sa
   integer(I4P)             :: o
   sa = 0._R_P
   do o=-2, 2
      sa = sa + abs(band(id, af, o))
   enddo
   f = (1._R_P + sa/band_min(id, af))/band_min(id, af)
   endfunction floor_factor

   pure function band_min(id, af) result(s)
   !< min over xi of |sum_o alpha_o e^(i o xi)|: the inverse of the band's amplification of the right side's errors.
   integer(I4P), intent(in) :: id
   real(R_P),    intent(in) :: af
   real(R_P)                :: s, re, im, xi
   integer(I4P)             :: j, o
   s = huge(1._R_P)
   do j=0, 512
      xi = j*PI/512
      re = 0._R_P ; im = 0._R_P
      do o=-2, 2
         re = re + band(id, af, o)*cos(o*xi) ; im = im + band(id, af, o)*sin(o*xi)
      enddo
      s = min(s, sqrt(re**2 + im**2))
   enddo
   endfunction band_min

   subroutine check_rows
   !< L1: every row exact on the monomials up to its degree, its residual beyond it the derived one.
   real(R_P)    :: q(-LWMAX:LWMAX), r, left, tmax, tol, a, b, af, sa
   real(R_P), parameter :: AFS(3) = [0._R_P, 0.3_R_P, 0.45_R_P]
   integer(I4P) :: id, j, o, kdeg, ia
   do id=1, NL
   do ia=1, merge(3, 1, lalpha_free(id))
      af = AFS(ia)
      kdeg = lorder(id) + ln(id) - 1
      do j=0, kdeg + 1
         q = 0._R_P
         do o=llo(id), lhi(id)
            select case(lgeom(id))
            case(1) ; q(o) = real(o, R_P)**j
            case(2) ; q(o) = (o - 0.5_R_P)**j
            case(3) ; a = o - 1._R_P ; b = o ; q(o) = (b**(j + 1) - a**(j + 1))/(j + 1)
            endselect
         enddo
         left = 0._R_P ; tmax = 0._R_P ; sa = 0._R_P
         do o=-2, 2
            left = left + band(id, af, o)*dmono(ln(id), j, real(o, R_P))
            tmax = max(tmax, abs(dmono(ln(id), j, real(o, R_P))))
            sa = sa + abs(band(id, af, o))
         enddo
         call line_rhs(id, q, 1._R_P, af, r)
         tol = 64*epsilon(1._R_P)*(sumabs(id, af)*max(maxval(abs(q)), 1._R_P) + sa*tmax)
         checks = checks + 1
         if (j <= kdeg) then
            if (abs(left - r) > tol) then
               fails = fails + 1
               print '(A,A,A,I0,A,ES10.2)', 'ROW NOT EXACT ', trim(lname(id)), ' on x^', j, ': ', abs(left - r)
            endif
         elseif (abs((left - r) - (lerror(id) + af*lerror_alpha(id))) > tol) then
            fails = fails + 1
            print '(A,A,A,I0,2ES12.4)', 'ROW LEADING ERROR ', trim(lname(id)), ' on x^', j, left - r, lerror(id)
         endif
      enddo
   enddo
   enddo
   endsubroutine check_rows

   subroutine periodic_data(id, n, h, q)
   !< The data of g on a periodic line of n unknowns, the unknown i at X0 + (i-1) h.
   integer(I4P), intent(in)  :: id, n
   real(R_P),    intent(in)  :: h
   real(R_P),    intent(out) :: q(n)
   real(R_P)                 :: x
   integer(I4P)              :: m
   do m=1, n
      x = X0 + (m - 1)*h
      select case(lgeom(id))
      case(1) ; q(m) = gp(0, x)
      case(2) ; q(m) = gp(0, x - h/2)
      case(3) ; q(m) = (gp(-1, x) - gp(-1, x - h))/h
      endselect
   enddo
   endsubroutine periodic_data

   subroutine apply_line(id, n, h, q, u, w)
   !< The operator on a periodic line: the right side point by point, then the band solved.
   integer(I4P), intent(in)    :: id, n
   real(R_P),    intent(in)    :: h, q(n)
   real(R_P),    intent(out)   :: u(n)
   real(R_P),    intent(inout) :: w(n,4)
   real(R_P)                   :: win(-LWMAX:LWMAX)
   integer(I4P)                :: i, o
   do i=1, n
      win = 0._R_P
      do o=llo(id), lhi(id)
         win(o) = q(modulo(i + o - 1, n) + 1)
      enddo
      call line_rhs(id, win, h, AF, u(i))
   enddo
   call line_solve(id, n, AF, u, w(:,1), w(:,2), w(:,3), w(:,4))
   endsubroutine apply_line

   subroutine check_orders
   !< L2: the order of every operator on a periodic line.
   real(R_P), allocatable :: q(:), u(:), w(:,:)
   real(R_P)    :: e(KMAX), noise(KMAX), p(KMAX), h, ex, emax, ff
   logical      :: ok(KMAX)
   integer(I4P) :: id, k, n, i, kbest, unresolved
   unresolved = 0
   do id=1, NL
      ff = floor_factor(id, AF)
      ok = .false. ; e = 0._R_P
      do k=1, KMAX
         n = 2**(k + 2)
         h = 2*PI/n
         allocate(q(n), u(n), w(n,4))
         call periodic_data(id, n, h, q)
         call apply_line(id, n, h, q, u, w)
         emax = 0._R_P ; ex = 0._R_P
         do i=1, n
            emax = max(emax, abs(gp(ln(id), X0 + (i - 1)*h)))
            ex = max(ex, abs(u(i) - gp(ln(id), X0 + (i - 1)*h)))
         enddo
         e(k) = ex
         noise(k) = epsilon(1._R_P)*(sumabs(id, AF)*(maxval(abs(q)) + merge(2/h, 0._R_P, lgeom(id) == 3))/h**ln(id)*ff &
                                     + emax)
         ok(k) = e(k) > 10*noise(k) .and. 3*(lhi(id) - llo(id) + 1)*h <= 1._R_P
         deallocate(q, u, w)
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
         unresolved = unresolved + 1
         print '(A,A,A,I0,A)', 'unresolved in this precision: ', trim(lname(id)), ' (order ', lorder(id), ')'
         if (precision(1._R_P) >= 30) then
            fails = fails + 1
            print '(A,A)', 'UNRESOLVED IN QUADRUPLE PRECISION ', trim(lname(id))
         endif
         cycle
      endif
      if (min(p(kbest), p(kbest+1)) < lorder(id) - 0.25_R_P) then
         fails = fails + 1
         print '(A,A,A,2F7.2,A,I0)', 'LINE ORDER ', trim(lname(id)), ': observed ', p(kbest), p(kbest+1), &
                                     ', design ', lorder(id)
      endif
   enddo
   print '(A,I0,A,I0,A)', 'gideon_test_lines: orders on periodic lines, ', NL - unresolved, ' measured, ', unresolved, &
                          ' unresolved in this precision'
   endsubroutine check_orders

   subroutine check_symbols
   !< L4: the line solve of a Fourier mode is the right symbol over the left one.
   real(R_P)    :: qc(NSYM), qs(NSYM), uc(NSYM), us(NSYM), w(NSYM,4), win(-LWMAX:LWMAX), h, xi, rc, rs, lre, lim
   real(R_P)    :: pre, pim, tol, x, err, mwre, mwim, sfac
   integer(I4P) :: id, kk, i, o
   h = 2*PI/NSYM
   do id=1, NL
      tol = 64*epsilon(1._R_P)*sumabs(id, AF)/h**ln(id)*floor_factor(id, AF)
      do kk=1, NSYM/2 - 1
         xi = kk*h
         do i=1, NSYM ! the mode cos(k x), sin(k x) through the geometry of the data, the unknown 1 at x = 0
            x = (i - 1)*h
            select case(lgeom(id))
            case(1) ; qc(i) = cos(kk*x) ; qs(i) = sin(kk*x)
            case(2) ; qc(i) = cos(kk*(x - h/2)) ; qs(i) = sin(kk*(x - h/2))
            case(3) ; qc(i) = (sin(kk*x) - sin(kk*(x - h)))/(kk*h) ; qs(i) = (cos(kk*(x - h)) - cos(kk*x))/(kk*h)
            endselect
         enddo
         ! the right symbol: the kernel on the mode at the unknown 1 (x = 0)
         win = 0._R_P
         do o=llo(id), lhi(id)
            win(o) = qc(modulo(o, NSYM) + 1)
         enddo
         call line_rhs(id, win, h, AF, rc)
         do o=llo(id), lhi(id)
            win(o) = qs(modulo(o, NSYM) + 1)
         enddo
         call line_rhs(id, win, h, AF, rs)
         lre = 0._R_P ; lim = 0._R_P
         do o=-2, 2
            lre = lre + band(id, AF, o)*cos(o*xi) ; lim = lim + band(id, AF, o)*sin(o*xi)
         enddo
         ! the response to e^(ikx) is (rc + i rs)/(lre + i lim) e^(ikx)
         pre = (rc*lre + rs*lim)/(lre**2 + lim**2)
         pim = (rs*lre - rc*lim)/(lre**2 + lim**2)
         call apply_line_mode(id, h, qc, uc, w)
         call apply_line_mode(id, h, qs, us, w)
         err = 0._R_P
         do i=1, NSYM
            x = (i - 1)*h
            err = max(err, abs(uc(i) - (pre*cos(kk*x) - pim*sin(kk*x))), abs(us(i) - (pim*cos(kk*x) + pre*sin(kk*x))))
         enddo
         checks = checks + 1
         if (err > tol) then
            fails = fails + 1
            print '(A,A,A,I0,A,ES10.2)', 'SYMBOL ', trim(lname(id)), ' k = ', kk, ': ', err
         endif
         if (lalpha_free(id)) then ! a filter: its transfer function real, in [0, 1]
            checks = checks + 1
            if (abs(pim) > tol .or. pre < -tol .or. pre > 1 + tol) then
               fails = fails + 1
               print '(A,A,A,I0,2ES12.4)', 'TRANSFER ', trim(lname(id)), ' k = ', kk, pre, pim
            endif
         elseif (lgeom(id) == 1 .and. (ln(id) == 1 .or. ln(id) == 2)) then
            ! centred: d1's modified wavenumber -i h (pre + i pim) real, d2's -h^2 (pre + i pim) real
            sfac = h**ln(id)
            if (ln(id) == 1) then
               mwre = sfac*pim ; mwim = -sfac*pre
            else
               mwre = -sfac*pre ; mwim = -sfac*pim
            endif
            checks = checks + 1
            if (abs(mwim) > tol*sfac) then
               fails = fails + 1
               print '(A,A,A,I0,A,ES10.2)', 'DISSIPATIVE ', trim(lname(id)), ' k = ', kk, ': ', mwim
            endif
         endif
      enddo
      if (lalpha_free(id)) then ! the odd-even mode removed
         do i=1, NSYM
            qc(i) = real(1 - 2*modulo(i, 2), R_P)
         enddo
         call apply_line_mode(id, h, qc, uc, w)
         checks = checks + 1
         if (maxval(abs(uc)) > tol) then
            fails = fails + 1
            print '(A,A,ES12.4)', 'ODD-EVEN MODE NOT REMOVED ', trim(lname(id)), maxval(abs(uc))
         endif
      endif
   enddo
   endsubroutine check_symbols

   subroutine apply_line_mode(id, h, q, u, w)
   !< apply_line on the symbol's line.
   integer(I4P), intent(in)    :: id
   real(R_P),    intent(in)    :: h, q(NSYM)
   real(R_P),    intent(out)   :: u(NSYM)
   real(R_P),    intent(inout) :: w(NSYM,4)
   call apply_line(id, NSYM, h, q, u, w)
   endsubroutine apply_line_mode

   subroutine check_device
   !< L5: the right sides and the periodic line solves on the device against the host.
   integer(I4P), parameter :: NP = 16384_I4P, NLD = 512_I4P, NN = 128_I4P
   real(R_P), allocatable  :: q(:,:), rd(:), rh(:), ql(:,:), ud(:,:), uh(:,:), w(:,:,:), win(:,:)
   real(R_P)               :: ds, tol, worst
   integer(I4P)            :: id, i, l, o, nbit, nbitl
   integer, allocatable    :: seed(:)
   integer                 :: ns
   call random_seed(size=ns)
   allocate(seed(ns)) ; seed = 1118 ; call random_seed(put=seed)
   allocate(q(-LWMAX:LWMAX,NP), rd(NP), rh(NP), ql(NN,NLD), ud(NN,NLD), uh(NN,NLD), w(NN,4,NLD), win(-LWMAX:LWMAX,NLD))
   call random_number(q) ; q = 2*q - 1
   call random_number(ql) ; ql = 2*ql - 1
   ds = 0.37_R_P
   nbit = 0 ; nbitl = 0 ; worst = 0._R_P
   do id=1, NL
      !$acc parallel loop copyin(q) copyout(rd)
      do i=1, NP
         call line_rhs_dev(id, q(:,i), ds, AF, rd(i))
      enddo
      do i=1, NP
         call line_rhs(id, q(:,i), ds, AF, rh(i))
      enddo
      tol = 64*epsilon(1._R_P)*sumabs(id, AF)/ds**ln(id)
      do i=1, NP
         checks = checks + 1
         if (rd(i) == rh(i)) nbit = nbit + 1
         if (abs(rd(i) - rh(i)) > tol*maxval(abs(q(:,i)))) then
            fails = fails + 1
            if (fails <= 10) print '(A,A,1X,I0,2ES25.16)', 'DEVICE RHS ', trim(lname(id)), i, rd(i), rh(i)
         endif
      enddo
      ! whole lines: one line per thread, the right side and the solve
      !$acc parallel loop copyin(ql) copyout(ud) create(w, win)
      do l=1, NLD
         !$acc loop seq
         do i=1, NN
            !$acc loop seq
            do o=llo(id), lhi(id)
               win(o,l) = ql(modulo(i + o - 1, NN) + 1,l)
            enddo
            call line_rhs_dev(id, win(:,l), ds, AF, ud(i,l))
         enddo
         call line_solve(id, NN, AF, ud(:,l), w(:,1,l), w(:,2,l), w(:,3,l), w(:,4,l))
      enddo
      do l=1, NLD
         call apply_line(id, NN, ds, ql(:,l), uh(:,l), w(:,:,l))
      enddo
      tol = tol*floor_factor(id, AF)
      do l=1, NLD
         checks = checks + 1
         if (all(ud(:,l) == uh(:,l))) nbitl = nbitl + 1
         worst = max(worst, maxval(abs(ud(:,l) - uh(:,l)))/(sumabs(id, AF)/ds**ln(id)*floor_factor(id, AF)))
         if (maxval(abs(ud(:,l) - uh(:,l))) > tol) then
            fails = fails + 1
            if (fails <= 10) print '(A,A,1X,I0)', 'DEVICE LINE ', trim(lname(id)), l
         endif
      enddo
   enddo
#if defined _OPENACC
   print '(A,I0,A,I0,A,I0,A,I0,A,ES9.2)', 'gideon_test_lines: device, right sides bitwise ', nbit, ' of ', NP*NL, &
      ', lines bitwise ', nbitl, ' of ', NLD*NL, ', worst relative ', worst
#else
   checks = checks + 1
   if (nbit /= NP*NL .or. nbitl /= NLD*NL) then
      fails = fails + 1
      print '(A)', 'HOST LOOPS NOT BITWISE'
   endif
#endif
   endsubroutine check_device
endprogram gideon_test_lines
