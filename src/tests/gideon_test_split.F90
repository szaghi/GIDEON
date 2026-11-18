!< GIDEON test (#2, 1.2): the upwind-central split of every upwind-biased face value, on the compiled kernels.
program gideon_test_split
!< GIDEON test (#2, 1.2): the upwind-central split of every upwind-biased face value, on the compiled kernels.
!<
!< The generator proves, in exact rationals, gideon_<data>_face_b<s><p> = gideon_<data>_face_c<p+1> +
!< gideon_<data>_face_diss_b<s><p> for p = 3, 5, 7, 9, both sides, point values and cell averages. This test checks
!< the compiled kernels on random windows: a hybrid that computes c + sigma*d is the centred value at sigma = 0
!< (bitwise: a sensor switching the dissipation off leaves the central scheme untouched) and the biased value at
!< sigma = 1 within the round-off of the two sums, 64 eps (sum|b| + sum|c|) max|q| (d = b - c: sum|d| is no larger).
!< Frozen bound, written before the kernels were run.
use penf, only : I4P, R_P
use gideon
use gideon_test_registry, only : WMAX, NK, kname, ksumabs
implicit none
real(R_P)    :: q(-WMAX:WMAX), b, c, d
integer(I4P) :: trial, fails, checks

fails = 0 ; checks = 0
call random_seed_fixed
do trial=1, 50
   call random_number(q) ; q = 2*q - 1
   ! point values
   call gideon_fd_face_bl3(q(-1:1), b) ; call gideon_fd_face_c4(q(-1:2), c) ; call gideon_fd_face_diss_bl3(q(-1:2), d)
   call check('fd_face_bl3', 'fd_face_c4')
   call gideon_fd_face_br3(q(0:2), b) ; call gideon_fd_face_diss_br3(q(-1:2), d)
   call check('fd_face_br3', 'fd_face_c4')
   call gideon_fd_face_bl5(q(-2:2), b) ; call gideon_fd_face_c6(q(-2:3), c) ; call gideon_fd_face_diss_bl5(q(-2:3), d)
   call check('fd_face_bl5', 'fd_face_c6')
   call gideon_fd_face_br5(q(-1:3), b) ; call gideon_fd_face_diss_br5(q(-2:3), d)
   call check('fd_face_br5', 'fd_face_c6')
   call gideon_fd_face_bl7(q(-3:3), b) ; call gideon_fd_face_c8(q(-3:4), c) ; call gideon_fd_face_diss_bl7(q(-3:4), d)
   call check('fd_face_bl7', 'fd_face_c8')
   call gideon_fd_face_br7(q(-2:4), b) ; call gideon_fd_face_diss_br7(q(-3:4), d)
   call check('fd_face_br7', 'fd_face_c8')
   call gideon_fd_face_bl9(q(-4:4), b) ; call gideon_fd_face_c10(q(-4:5), c) ; call gideon_fd_face_diss_bl9(q(-4:5), d)
   call check('fd_face_bl9', 'fd_face_c10')
   call gideon_fd_face_br9(q(-3:5), b) ; call gideon_fd_face_diss_br9(q(-4:5), d)
   call check('fd_face_br9', 'fd_face_c10')
   ! cell averages
   call gideon_fv_face_bl3(q(-1:1), b) ; call gideon_fv_face_c4(q(-1:2), c) ; call gideon_fv_face_diss_bl3(q(-1:2), d)
   call check('fv_face_bl3', 'fv_face_c4')
   call gideon_fv_face_br3(q(0:2), b) ; call gideon_fv_face_diss_br3(q(-1:2), d)
   call check('fv_face_br3', 'fv_face_c4')
   call gideon_fv_face_bl5(q(-2:2), b) ; call gideon_fv_face_c6(q(-2:3), c) ; call gideon_fv_face_diss_bl5(q(-2:3), d)
   call check('fv_face_bl5', 'fv_face_c6')
   call gideon_fv_face_br5(q(-1:3), b) ; call gideon_fv_face_diss_br5(q(-2:3), d)
   call check('fv_face_br5', 'fv_face_c6')
   call gideon_fv_face_bl7(q(-3:3), b) ; call gideon_fv_face_c8(q(-3:4), c) ; call gideon_fv_face_diss_bl7(q(-3:4), d)
   call check('fv_face_bl7', 'fv_face_c8')
   call gideon_fv_face_br7(q(-2:4), b) ; call gideon_fv_face_diss_br7(q(-3:4), d)
   call check('fv_face_br7', 'fv_face_c8')
   call gideon_fv_face_bl9(q(-4:4), b) ; call gideon_fv_face_c10(q(-4:5), c) ; call gideon_fv_face_diss_bl9(q(-4:5), d)
   call check('fv_face_bl9', 'fv_face_c10')
   call gideon_fv_face_br9(q(-3:5), b) ; call gideon_fv_face_diss_br9(q(-4:5), d)
   call check('fv_face_br9', 'fv_face_c10')
enddo
print '(A,I0,A,I0)', 'gideon_test_split: 16 biased kernels, ', checks, ' checks, failures ', fails
if (fails > 0) error stop 'gideon_test_split: failed'
print '(A)', 'gideon_test_split: ok'
contains
   subroutine check(bname, cname)
   !< b = c + d within round-off, c + 0*d = c bitwise (b, c, d and q of the host).
   character(*), intent(in) :: bname, cname !< The biased and the centred kernel (registry names).
   real(R_P)                :: tol, sigma
   sigma = 0._R_P
   checks = checks + 2
   if (c + sigma*d /= c) then
      fails = fails + 1
      print '(A,A)', 'SIGMA 0 NOT CENTRED ', bname
   endif
   tol = 64*epsilon(1._R_P)*(ksumabs(id_of(bname)) + ksumabs(id_of(cname)))*maxval(abs(q))
   if (abs(c + d - b) > tol) then
      fails = fails + 1
      print '(A,A,A,ES12.4,A,ES12.4)', 'SPLIT ', bname, ': |c + d - b| = ', abs(c + d - b), ' > ', tol
   endif
   endsubroutine check

   function id_of(name) result(id)
   !< The registry id of a kernel by its name (without the gideon_ prefix).
   character(*), intent(in) :: name
   integer(I4P)             :: id
   do id=1, NK
      if (trim(kname(id)) == name) return
   enddo
   error stop 'gideon_test_split: no such kernel'
   endfunction id_of

   subroutine random_seed_fixed
   !< A fixed seed: the same data on every run.
   integer, allocatable :: seed(:)
   integer              :: n
   call random_seed(size=n)
   allocate(seed(n))
   seed = 20261118
   call random_seed(put=seed)
   endsubroutine random_seed_fixed
endprogram gideon_test_split
