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
!< A kernel with no plateau in the working precision is listed as unresolved, not failed: the widest upwind stencils
!< (7 to 9 cells) leave their pre-asymptotic range in double precision only where round-off already dominates,
!< whatever the function (truncation over round-off depends on the wavenumber times h alone). Their orders, and
!< every other row's, are measured round-off free by the generator (scripts/gen_gideon.py, 100-digit arithmetic,
!< same bound), and the compiled rows are proven exact to their degree in every kind (gideon_test_tables).
use penf, only : I4P, R_P
use gideon_test_registry
use gideon_test_functions, only : data_generic, target_generic
implicit none
integer(I4P), parameter :: KMAX = 20
real(R_P)    :: q(-WMAX:WMAX), r, t, h, e(KMAX), noise(KMAX), p(KMAX)
logical      :: ok(KMAX)
integer(I4P) :: id, k, kbest, fails, unresolved

fails = 0 ; unresolved = 0
do id=1, NK
   ok = .false.
   e = 0._R_P
   do k=1, KMAX
      h = 0.5_R_P**k
      ! the step's scales 1/h^(n+1) must fit the working kind (R4P: 2^127), or a trapping build stops on overflow
      if (k*(kn(id) + 1) > 0.8_R_P*maxexponent(1._R_P)) exit
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
print '(A,I0,A,I0,A,I0,A,I0)', 'gideon_test_order: ', NK, ' kernels, ', NK - unresolved, ' measured, ', unresolved, &
                               ' unresolved in this precision, failures ', fails
if (fails > 0) error stop 'gideon_test_order: failed'
print '(A)', 'gideon_test_order: ok'
endprogram gideon_test_order
