!< GIDEON test (G3): the compiled kernels reach their design orders on a generic smooth function, where the
!< working precision can show it.
program gideon_test_order
!< GIDEON test (G3): the compiled kernels reach their design orders on a generic smooth function, where the
!< working precision can show it.
!<
!< g(x) = sin(1.3 x + 0.4) + 0.5 cos(0.7 x - 0.2) at x0 = 0.37 (gideon_test_functions), h = 2^-1 ... 2^-20. The
!< round-off floor of a kernel at step h is bounded from above: eps sum_k |c_k| (max|q| + 2/h for cell averages,
!< built as a difference of primitives over h) / h^n, plus eps |target| and, for a cell-average target, the 2 eps/h
!< of the difference that builds it. A step whose error exceeds 10 floors, and whose window spans less than a twelfth
!< of the shorter wavelength (1.3 width h <= 0.5), is resolved. A plateau is three consecutive resolved steps whose two
!< observed orders log2(e(h)/e(h/2)) differ by at most 0.25: the asymptotic range. Frozen
!< bound (written before the kernels): on the finest plateau, both observed orders >= design - 0.25.
!<
!< A kernel with no plateau in the working precision is listed as unresolved, not failed: in double precision the
!< widest upwind stencils (7 to 9 cells) leave their pre-asymptotic range only where round-off already dominates,
!< whatever the function (truncation over round-off depends on the wavenumber times h alone): R8P resolves 119 of the
!< 288 kernels, R4P 17. The quadruple-precision build (tests-gnu-r16, host only) resolves all 288: that is where the
!< orders of the compiled kernels are observed, every one. The generator measures the rows themselves at 100 digits.
!<
!< Report mode (scripts/run_report.sh, the documentation's data): with the environment variable GIDEON_REPORT_DIR set,
!< the test also writes order-rNN.csv (every kernel, every step: h, e(h), the round-off floor, resolved or not) and
!< order-summary-rNN.csv (every kernel: the plateau and its two observed orders) into that directory, NN the bits of
!< the working kind. The test itself is unchanged.
use penf, only : I4P, R8P, R_P
use gideon_test_registry
use gideon_test_functions, only : data_generic, target_generic
implicit none
integer(I4P), parameter :: KMAX = 20
real(R_P)    :: q(-WMAX:WMAX), r, t, h, e(KMAX), noise(KMAX), p(KMAX)
logical      :: ok(KMAX)
integer(I4P) :: id, k, kbest, fails, unresolved
integer(I4P) :: ksteps
character(4096) :: report_dir
integer      :: ulen, ust, usum, ios
logical      :: report

call get_environment_variable('GIDEON_REPORT_DIR', report_dir, length=ulen)
report = ulen > 0
if (report) then
   open(newunit=ust, file=trim(report_dir)//'/order-r'//kind_tag()//'.csv', action='write', status='replace', &
        iostat=ios)
   if (ios /= 0) error stop 'gideon_test_order: cannot write into GIDEON_REPORT_DIR'
   write(ust, '(A)') 'kernel,design,k,h,e,floor,resolved'
   open(newunit=usum, file=trim(report_dir)//'/order-summary-r'//kind_tag()//'.csv', action='write', &
        status='replace')
   write(usum, '(A)') 'kernel,design,status,kbest,p1,p2'
endif
fails = 0 ; unresolved = 0
do id=1, NK
   ok = .false.
   e = 0._R_P
   ksteps = 0
   do k=1, KMAX
      h = 0.5_R_P**k
      ! the step's scales 1/h^(n+1) must fit the working kind (R4P: 2^127), or a trapping build stops on overflow
      if (k*(kn(id) + 1) > 0.8_R_P*maxexponent(1._R_P)) exit
      ksteps = k
      call data_generic(id, h, q)
      call apply(id, q, h, r)
      t = target_generic(id, h)
      e(k) = abs(r - t)
      noise(k) = epsilon(1._R_P)*(ksumabs(id)*(maxval(abs(q)) + merge(2/h, 0._R_P, kgeom(id) >= 3))/h**kn(id) + &
                                  abs(t) + merge(2/h, 0._R_P, ktarget(id) == 2))
      ok(k) = e(k) > 10*noise(k) .and. 1.3_R_P*(khi(id) - klo(id) + 1)*h <= 0.5_R_P
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
   if (report) then
      do k=1, ksteps
         write(ust, '(A,",",I0,",",I0,",",ES16.8E4,",",ES16.8E4,",",ES16.8E4,",",L1)') trim(kname(id)), korder(id), &
            k, real(0.5_R_P**k, R8P), real(e(k), R8P), real(noise(k), R8P), ok(k)
      enddo
      if (kbest == 0) then
         write(usum, '(A,",",I0,",unresolved,0,,")') trim(kname(id)), korder(id)
      else
         write(usum, '(A,",",I0,",measured,",I0,",",F0.4,",",F0.4)') trim(kname(id)), korder(id), kbest, &
            real(p(kbest), R8P), real(p(kbest+1), R8P)
      endif
   endif
   if (kbest == 0) then
      unresolved = unresolved + 1
      print '(A,A,A,I0,A)', 'unresolved in this precision: ', trim(kname(id)), ' (order ', korder(id), ')'
      cycle
   endif
   if (min(p(kbest), p(kbest+1)) < korder(id) - 0.25_R_P) then
      fails = fails + 1
      print '(A,A,A,2F7.2,A,I0,A,ES9.2)', 'ORDER ', trim(kname(id)), ': observed ', p(kbest), p(kbest+1), &
                                          ', design ', korder(id), ' at h = ', 0.5_R_P**kbest
   endif
enddo
if (report) then
   close(ust) ; close(usum)
endif
print '(A,I0,A,I0,A,I0,A,I0)', 'gideon_test_order: ', NK, ' kernels, ', NK - unresolved, ' measured, ', unresolved, &
                               ' unresolved in this precision, failures ', fails
if (fails > 0) error stop 'gideon_test_order: failed'
print '(A)', 'gideon_test_order: ok'
contains
   function kind_tag() result(tag)
   !< The bits of the working kind, the suffix of the report files: 32, 64 or 128.
   character(:), allocatable :: tag
   character(8)              :: buf
   write(buf, '(I0)') storage_size(1._R_P)
   tag = trim(buf)
   endfunction kind_tag
endprogram gideon_test_order
