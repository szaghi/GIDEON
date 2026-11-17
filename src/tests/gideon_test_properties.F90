!< GIDEON test (G4): the mirror identity, the (anti)symmetries and the exact scalings, bitwise.
program gideon_test_properties
!< GIDEON test (G4): the mirror identity, the (anti)symmetries and the exact scalings, bitwise.
!<
!< - mirror: every left-side upwind kernel on mirrored data equals (-1)^n its right-side twin on the data, bitwise
!<   (the left rows are the right rows read mirrored, never stored: audit item A cannot recur);
!< - symmetry: a centred kernel whose result is odd in the data's symmetry returns exactly 0 (pair form);
!< - constant data: the pair and flux-difference forms return exactly 0 for a derivative; the others the exact
!<   target within 64 eps sum|c|;
!< - scaling: kernel(8 q) = 8 kernel(q) and kernel(q, 2 ds) = kernel(q, ds)/2^n, bitwise (powers of two).
use penf, only : I4P, R_P
use gideon_test_registry
implicit none
real(R_P)    :: q(-WMAX:WMAX), qm(-WMAX:WMAX), r, rm, t
integer(I4P) :: id, jr, o, m, fails, trial
real(R_P)    :: u(-WMAX:WMAX)

fails = 0
call random_seed_fixed
do trial=1, 20
   call random_number(u)
   q = u - 0.5_R_P
   ! mirror identity
   do id=1, NK
      if (kmirror(id) == 0) cycle
      jr = kmirror(id)
      qm = 0._R_P
      do o=-WMAX + 1, WMAX - 1
         if (kgeom(id) == 2 .or. kgeom(id) == 3) then ! face at 0 (points at o - 1/2, or cells [o-1, o])
            qm(o) = q(1 - o)
         else
            qm(o) = q(-o)
         endif
      enddo
      call apply(jr, q, 1._R_P, r)
      call apply(id, qm, 1._R_P, rm)
      if (rm /= (-1)**kn(id)*r) then
         fails = fails + 1
         print '(A,A,A,A,2ES24.16)', 'MIRROR ', trim(kname(id)), ' vs ', trim(kname(jr)), rm, r
      endif
   enddo
   ! symmetry of centred kernels: odd results vanish exactly
   do id=1, NK
      if (kbias(id) /= 'c ') cycle
      qm = 0._R_P
      do m=0, WMAX - 1
         if (kgeom(id) == 2 .or. kgeom(id) == 3) then ! face at 0: pairs (m+1, -m)
            if (mod(kn(id), 2) == 1) then
               qm(m+1) = q(m+1) ; qm(-m) = q(m+1)
            else
               qm(m+1) = q(m+1) ; qm(-m) = -q(m+1)
            endif
         else ! node or cell 0: pairs (m, -m)
            if (mod(kn(id), 2) == 1) then
               qm(m) = q(m) ; qm(-m) = q(m)
            else
               qm(m) = q(m) ; qm(-m) = -q(m)
               if (m == 0) qm(0) = 0._R_P
            endif
         endif
      enddo
      call apply(id, qm, 1._R_P, r)
      if (r /= 0._R_P) then
         fails = fails + 1
         print '(A,A,ES12.4)', 'SYMMETRY ', trim(kname(id)), r
      endif
   enddo
   ! scalings
   do id=1, NK
      call apply(id, q, 0.75_R_P, r)
      call apply(id, 8*q, 0.75_R_P, rm)
      if (rm /= 8*r) then
         fails = fails + 1
         print '(A,A)', 'SCALE DATA ', trim(kname(id))
      endif
      call apply(id, q, 1.5_R_P, rm)
      if (rm /= r/2._R_P**kn(id)) then
         fails = fails + 1
         print '(A,A)', 'SCALE STEP ', trim(kname(id))
      endif
   enddo
enddo
! constant data
do id=1, NK
   q = 0.3_R_P
   call apply(id, q, 1._R_P, r)
   t = merge(0.3_R_P, 0._R_P, kn(id) == 0)
   if (kzero_const(id) .and. kn(id) >= 1) then
      if (r /= 0._R_P) then
         fails = fails + 1
         print '(A,A,ES12.4)', 'CONSTANT (exact) ', trim(kname(id)), r
      endif
   elseif (abs(r - t) > 64*epsilon(1._R_P)*ksumabs(id)*0.3_R_P) then
      fails = fails + 1
      print '(A,A,ES12.4)', 'CONSTANT ', trim(kname(id)), r - t
   endif
enddo
print '(A,I0,A,I0)', 'gideon_test_properties: ', NK, ' kernels, failures ', fails
if (fails > 0) error stop 'gideon_test_properties: failed'
print '(A)', 'gideon_test_properties: ok'
contains
   subroutine random_seed_fixed
   !< A fixed seed: the same data on every run.
   integer, allocatable :: seed(:)
   integer              :: n
   call random_seed(size=n)
   allocate(seed(n))
   seed = 20261117
   call random_seed(put=seed)
   endsubroutine random_seed_fixed
endprogram gideon_test_properties
