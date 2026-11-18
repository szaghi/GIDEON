!< GIDEON test (G7 of the compact operators): regression against stored golden values of every right side and line solve.
program gideon_test_lines_golden
!< GIDEON test (G7 of the compact operators): regression against stored golden values of every right side and line solve.
!<
!< As gideon_test_golden, for the compact operators (gideon_test_lines_registry): every operator's right side on NST
!< deterministic dyadic windows, and one deterministic periodic line of NLINE unknowns solved with its band (a filter's
!< alpha_f = AF), ds = 0.75. A value passes when it differs from the golden one by at most TOL of its magnitude,
!< sum|w| max|q| / ds^n, times F = (1 + sum|alpha|/s)/s for a solved value (s the band's smallest symbol, the
!< amplification of the round-off through the solve): other compilers contract multiply-adds differently.
!<
!<   gideon_test_lines_golden                 compare with the golden file of the build's kind (from the repository root)
!<   gideon_test_lines_golden write <file>    write a new golden file (refuses to overwrite an existing one)
!<
!< A golden file is never edited: a deliberate change of results gets a new file and this test's default path.
!< v1 (2026-11-18, #2): the 29 compact operators, written by the optimised gfortran builds.
use, intrinsic :: iso_fortran_env, only : compiler_version, compiler_options
use penf, only : I4P, I8P, R_P
use gideon_test_lines_registry
implicit none
integer(I4P), parameter :: NST = 8_I4P, NLINE = 16_I4P
character(*), parameter :: VERSION = 'v1'
real(R_P),    parameter :: DS = 0.75_R_P, AF = 0.45_R_P
real(R_P),    parameter :: PI = 3.14159265358979323846264338327950288_R_P
#if defined _R_P_IS_R4P
character(*), parameter :: GOLDEN = 'src/tests/golden/gideon_golden_lines_r4_v1.txt'
real(R_P),    parameter :: TOL = 1.e-5_R_P
character(*), parameter :: FMT = '(A24,1X,I3,1X,ES16.8E3)'
#elif defined _R_P_IS_R16P
character(*), parameter :: GOLDEN = 'src/tests/golden/gideon_golden_lines_r16_v1.txt'
real(R_P),    parameter :: TOL = 1.e-30_R_P
character(*), parameter :: FMT = '(A24,1X,I3,1X,ES45.36E4)'
#else
character(*), parameter :: GOLDEN = 'src/tests/golden/gideon_golden_lines_v1.txt'
real(R_P),    parameter :: TOL = 1.e-13_R_P
character(*), parameter :: FMT = '(A24,1X,I3,1X,ES25.17E3)'
#endif
real(R_P)      :: q(-LWMAX:LWMAX,NST), ql(NLINE), vals(NST+NLINE), mags(NST+NLINE), g, worst
integer(I4P)   :: id, k, kk, u, ios, n_values, n_bitwise, n_fail
integer(I8P)   :: lcg
character(256) :: arg, path, line
character(24)  :: name

lcg = 20261118_I8P
do k=1, NST
   do kk=-LWMAX, LWMAX
      q(kk,k) = next()
   enddo
enddo
do k=1, NLINE
   ql(k) = next()
enddo
call get_command_argument(1, arg)
if (trim(arg) == 'write') then
   call get_command_argument(2, path)
   open(newunit=u, file=trim(path), status='new', action='write', iostat=ios)
   if (ios /= 0) error stop 'gideon_test_lines_golden: cannot create the file (it exists?): never overwrite a golden file'
   write(u, '(A)') '# GIDEON golden values of the compact operators '//VERSION//': operator, index, value. Written by '// &
                   compiler_version()
   write(u, '(A)') '# options: '//compiler_options()
   do id=1, NL
      call values(id, vals, mags)
      do k=1, NST + NLINE
         write(u, FMT) lname(id), k, vals(k)
      enddo
   enddo
   close(u)
   print '(A)', 'gideon_test_lines_golden: written '//trim(path)
   stop
endif
open(newunit=u, file=GOLDEN, status='old', action='read', iostat=ios)
if (ios /= 0) error stop 'gideon_test_lines_golden: golden file not found (run from the repository root)'
n_values = 0 ; n_bitwise = 0 ; n_fail = 0 ; worst = 0._R_P
read(u, '(A)') line ; read(u, '(A)') line
do id=1, NL
   call values(id, vals, mags)
   do k=1, NST + NLINE
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) error stop 'gideon_test_lines_golden: golden file too short'
      read(line, FMT) name, kk, g
      if (trim(name) /= trim(lname(id)) .or. kk /= k) error stop 'gideon_test_lines_golden: golden file out of order'
      n_values = n_values + 1
      if (vals(k) == g) n_bitwise = n_bitwise + 1
      worst = max(worst, abs(vals(k) - g)/mags(k))
      if (abs(vals(k) - g) > TOL*mags(k)) then
         n_fail = n_fail + 1
         print '(A,A,1X,I0,2ES25.16)', 'GOLDEN ', trim(lname(id)), k, vals(k), g
      endif
   enddo
enddo
close(u)
print '(A,I0,A,I0,A,ES9.2,A,I0)', 'gideon_test_lines_golden: ', n_values, ' values, bitwise ', n_bitwise, &
                                  ', worst relative difference ', worst, ', failures ', n_fail
if (n_fail > 0) error stop 'gideon_test_lines_golden: failed'
print '(A)', 'gideon_test_lines_golden: ok'
contains
   function next() result(v)
   !< The next dyadic datum: an integer LCG, multiples of 2^-10 in [-32, 32], exact in every kind.
   real(R_P) :: v
   lcg = modulo(1103515245_I8P*lcg + 12345_I8P, 2147483648_I8P)
   v = real(modulo(lcg/65536_I8P, 65536_I8P) - 32768_I8P, R_P)/1024._R_P
   endfunction next

   subroutine values(id, vals, mags)
   !< The right sides on the windows, then the solved line; their magnitudes.
   integer(I4P), intent(in)  :: id
   real(R_P),    intent(out) :: vals(NST+NLINE), mags(NST+NLINE)
   real(R_P)                 :: win(-LWMAX:LWMAX), w(NLINE,4), sw, sa, s, re, im, xi
   integer(I4P)              :: k, i, o, j
   sw = lsumabs(id) + merge(AF*lsumabs_alpha(id), 0._R_P, lalpha_free(id))
   do k=1, NST
      call line_rhs(id, q(:,k), DS, AF, vals(k))
      mags(k) = sw*maxval(abs(q(llo(id):lhi(id),k)))/DS**ln(id)
   enddo
   do i=1, NLINE
      win = 0._R_P
      do o=llo(id), lhi(id)
         win(o) = ql(modulo(i + o - 1, NLINE) + 1)
      enddo
      call line_rhs(id, win, DS, AF, vals(NST+i))
   enddo
   call line_solve(id, NLINE, AF, vals(NST+1:NST+NLINE), w(:,1), w(:,2), w(:,3), w(:,4))
   ! F: the band's amplification of the round-off
   sa = 0._R_P ; s = huge(1._R_P)
   do o=-2, 2
      sa = sa + abs(band(id, o))
   enddo
   do j=0, 256
      xi = j*PI/256
      re = 0._R_P ; im = 0._R_P
      do o=-2, 2
         re = re + band(id, o)*cos(o*xi) ; im = im + band(id, o)*sin(o*xi)
      enddo
      s = min(s, sqrt(re**2 + im**2))
   enddo
   mags(NST+1:) = sw*maxval(abs(ql))/DS**ln(id)*(1 + sa/s)/s
   endsubroutine values

   pure function band(id, o) result(a)
   !< The band at offset o (a filter's: [AF, 1, AF]).
   integer(I4P), intent(in) :: id, o
   real(R_P)                :: a
   if (lalpha_free(id)) then
      a = merge(1._R_P, merge(AF, 0._R_P, abs(o) == 1), o == 0)
   else
      a = lleft(o,id)
   endif
   endfunction band
endprogram gideon_test_lines_golden
