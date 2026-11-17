!< GIDEON test: every run-time order dispatcher calls the kernel of its order, and returns NaN for an order without one.
program gideon_test_dispatch
!< GIDEON test: every run-time order dispatcher calls the kernel of its order, and returns NaN for an order without one.
!<
!< For each kernel: dispatcher(order) on the window of that order equals the kernel on the same data, bitwise; and
!< the dispatcher at an order with no kernel (gideon_test_registry's apply_dispatch_bad) returns GIDEON_NAN, never a
!< silent number (ADAM returned 0 for an unsupported order: audit item C).
use penf, only : I4P, R_P
use gideon_test_registry
implicit none
real(R_P)    :: q(-WMAX:WMAX), r, rd
integer(I4P) :: id, fails
integer, allocatable :: seed(:)
integer      :: n

call random_seed(size=n)
allocate(seed(n)) ; seed = 1117 ; call random_seed(put=seed)
call random_number(q)
fails = 0
do id=1, NK
   call apply(id, q, 0.37_R_P, r)
   call apply_dispatch(id, q, 0.37_R_P, rd)
   if (rd /= r) then
      fails = fails + 1
      print '(A,A,2ES24.16)', 'DISPATCH ', trim(kname(id)), rd, r
   endif
   call apply_dispatch_bad(id, q, 0.37_R_P, rd)
   if (rd == rd) then
      fails = fails + 1
      print '(A,A,ES24.16)', 'NOT NAN ', trim(kname(id)), rd
   endif
enddo
print '(A,I0,A,I0)', 'gideon_test_dispatch: ', NK, ' kernels, failures ', fails
if (fails > 0) error stop 'gideon_test_dispatch: failed'
print '(A)', 'gideon_test_dispatch: ok'
endprogram gideon_test_dispatch
