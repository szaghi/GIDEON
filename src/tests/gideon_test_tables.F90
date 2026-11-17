!< GIDEON test (G2): every compiled kernel is exact on the monomials up to its degree, and not beyond.
program gideon_test_tables
!< GIDEON test (G2): every compiled kernel is exact on the monomials up to its degree, and not beyond.
!<
!< For each kernel and each x^j, j <= K (its exact degree, K = order + n - 1): |kernel - target| <= tol, tol =
!< 64 eps sum_k |c_k| max|q|, the round-off of the dot product. For x^(K+1): the error equals the exact leading error
!< the generator derived (kerror, non-zero) within the same tol, so the compiled row is neither less nor more exact
!< than derived, in every kind. Frozen bound (written before the kernels): 64 eps.
use penf, only : I4P, R_P
use gideon_test_registry
use gideon_test_functions, only : data_monomial, target_monomial
implicit none
real(R_P)    :: q(-WMAX:WMAX), r, t, tol
integer(I4P) :: id, j, fails, checks

fails = 0 ; checks = 0
do id=1, NK
   do j=0, kdegree(id) + 1
      call data_monomial(id, j, q)
      call apply(id, q, 1._R_P, r)
      t = target_monomial(id, j)
      tol = 64*epsilon(1._R_P)*ksumabs(id)*max(maxval(abs(q)), 1._R_P)
      checks = checks + 1
      if (j <= kdegree(id)) then
         if (abs(r - t) > tol) then
            fails = fails + 1
            print '(A,A,A,I0,A,ES10.2,A,ES10.2)', 'NOT EXACT ', trim(kname(id)), ' on x^', j, ': error ', abs(r - t), &
                                                 ' > tol ', tol
         endif
      else
         if (abs((r - t) - kerror(id)) > tol) then
            fails = fails + 1
            print '(A,A,A,I0,A,ES12.4,A,ES12.4)', 'LEADING ERROR ', trim(kname(id)), ' on x^', j, ': ', r - t, &
                                                  ' derived ', kerror(id)
         endif
      endif
   enddo
enddo
print '(A,I0,A,I0,A,I0)', 'gideon_test_tables: ', NK, ' kernels, ', checks, ' checks, failures ', fails
if (fails > 0) error stop 'gideon_test_tables: failed'
print '(A)', 'gideon_test_tables: ok'
endprogram gideon_test_tables
