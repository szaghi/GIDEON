!< GIDEON test (G7): regression against stored golden values of every kernel on deterministic data.
program gideon_test_golden
!< GIDEON test (G7): regression against stored golden values of every kernel on deterministic data.
!<
!< The other tests assert exactness, orders and properties: a change of a few ulp, or a wrong but still exact row of
!< a kernel outside the tests' reach, passes them. This test pins the values. Every kernel of the registry on NST
!< deterministic windows: an integer LCG makes dyadic data (multiples of 2^-10 in [-32, 32]), exact in every kind and
!< compiler, ds = 0.75. A value passes when it differs from the golden one by at most TOL of its magnitude
!< (sum_k |c_k| max|q| / ds^n): other compilers and flags contract multiply-adds differently from the build that
!< wrote the file (its compiler and options are in the file's header). The counts of bitwise-equal values and the
!< worst difference are printed.
!<
!<   gideon_test_golden                  compare with the golden file of the build's kind (run from the repository root)
!<   gideon_test_golden write <file>     write a new golden file (refuses to overwrite an existing one)
!<
!< A golden file is never edited: a deliberate change of results gets a new file and this test's default path.
!< v1 (2026-11-17): the 243 kernels of the first release, written by the optimised gfortran build; r16_v1 the same in
!< quadruple precision (the host-only verification kind).
use, intrinsic :: iso_fortran_env, only : compiler_version, compiler_options
use penf, only : I4P, I8P, R4P, R8P, R_P
use gideon_test_registry
implicit none
integer(I4P), parameter :: NST = 8_I4P
real(R_P),    parameter :: DS = 0.75_R_P
#if defined _R_P_IS_R4P
character(*), parameter :: GOLDEN = 'src/tests/golden/gideon_golden_r4_v1.txt'
real(R_P),    parameter :: TOL = 1.e-5_R_P
character(*), parameter :: FMT = '(A24,1X,I2,1X,ES16.8E3)'
#elif defined _R_P_IS_R16P
character(*), parameter :: GOLDEN = 'src/tests/golden/gideon_golden_r16_v1.txt'
real(R_P),    parameter :: TOL = 1.e-30_R_P
character(*), parameter :: FMT = '(A24,1X,I2,1X,ES45.36E4)'
#else
character(*), parameter :: GOLDEN = 'src/tests/golden/gideon_golden_v1.txt'
real(R_P),    parameter :: TOL = 1.e-13_R_P
character(*), parameter :: FMT = '(A24,1X,I2,1X,ES25.17E3)'
#endif
real(R_P)      :: q(-WMAX:WMAX,NST), r, g, mag, worst
integer(I4P)   :: id, k, u, ios, n_values, n_bitwise, n_fail, kk
integer(I8P)   :: lcg
character(256) :: arg, path, line
character(24)  :: name

! deterministic dyadic data
lcg = 20261117_I8P
do k=1, NST
   do kk=-WMAX, WMAX
      lcg = modulo(1103515245_I8P*lcg + 12345_I8P, 2147483648_I8P)
      q(kk,k) = real(modulo(lcg/65536_I8P, 65536_I8P) - 32768_I8P, R_P)/1024._R_P
   enddo
enddo
call get_command_argument(1, arg)
if (trim(arg) == 'write') then
   call get_command_argument(2, path)
   open(newunit=u, file=trim(path), status='new', action='write', iostat=ios)
   if (ios /= 0) error stop 'gideon_test_golden: cannot create the file (it exists?): never overwrite a golden file'
   write(u, '(A)') '# GIDEON golden values v1: kernel, window, value. Written by '//compiler_version()
   write(u, '(A)') '# options: '//compiler_options()
   do id=1, NK
      do k=1, NST
         call apply(id, q(:,k), DS, r)
         write(u, FMT) kname(id), k, r
      enddo
   enddo
   close(u)
   print '(A)', 'gideon_test_golden: written '//trim(path)
   stop
endif
open(newunit=u, file=GOLDEN, status='old', action='read', iostat=ios)
if (ios /= 0) error stop 'gideon_test_golden: golden file not found (run from the repository root)'
n_values = 0 ; n_bitwise = 0 ; n_fail = 0 ; worst = 0._R_P
read(u, '(A)') line ; read(u, '(A)') line
do id=1, NK
   do k=1, NST
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) error stop 'gideon_test_golden: golden file too short'
      read(line, FMT) name, kk, g
      if (trim(name) /= trim(kname(id)) .or. kk /= k) error stop 'gideon_test_golden: golden file out of order'
      call apply(id, q(:,k), DS, r)
      mag = ksumabs(id)*maxval(abs(q(klo(id):khi(id),k)))/DS**kn(id)
      n_values = n_values + 1
      if (r == g) n_bitwise = n_bitwise + 1
      worst = max(worst, abs(r - g)/mag)
      if (abs(r - g) > TOL*mag) then
         n_fail = n_fail + 1
         print '(A,A,1X,I0,2ES25.16)', 'GOLDEN ', trim(kname(id)), k, r, g
      endif
   enddo
enddo
close(u)
print '(A,I0,A,I0,A,ES9.2,A,I0)', 'gideon_test_golden: ', n_values, ' values, bitwise ', n_bitwise, &
                                  ', worst relative difference ', worst, ', failures ', n_fail
if (n_fail > 0) error stop 'gideon_test_golden: failed'
print '(A)', 'gideon_test_golden: ok'
endprogram gideon_test_golden
