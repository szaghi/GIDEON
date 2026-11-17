!< GIDEON test (G8): every kernel on the device equals the host within the multiply-add contraction tolerance.
program gideon_test_device
!< GIDEON test (G8): every kernel on the device equals the host within the multiply-add contraction tolerance.
!<
!< NP windows of random data, each kernel applied to all of them inside one `!$acc parallel loop` (apply_dev, the
!< registry's device-callable dispatcher: one select case per call, as a caller with a run-time choice would do),
!< against the host's apply. Tolerance 64 eps sum_k |c_k| max|q| / ds^n: the device and the host contract multiply-
!< adds differently. The number of bitwise-equal values is printed. Built without OpenACC the loop runs on the host,
!< and the test checks apply_dev against apply (bitwise).
use penf, only : I4P, R_P
use gideon_test_registry
implicit none
integer(I4P), parameter :: NP = 65536_I4P
real(R_P),    parameter :: DS = 0.37_R_P
real(R_P), allocatable  :: q(:,:), rd(:), rh(:)
real(R_P)               :: mag, worst
integer(I4P)            :: id, i, n_bitwise, n_values, n_fail
integer, allocatable    :: seed(:)
integer                 :: ns

call random_seed(size=ns)
allocate(seed(ns)) ; seed = 4096 ; call random_seed(put=seed)
allocate(q(-WMAX:WMAX,NP), rd(NP), rh(NP))
call random_number(q) ; q = 2*q - 1
n_bitwise = 0 ; n_values = 0 ; n_fail = 0 ; worst = 0._R_P
!$acc data copyin(q) create(rd)
do id=1, NK
   !$acc parallel loop present(q, rd)
   do i=1, NP
      call apply_dev(id, q(:,i), DS, rd(i))
   enddo
   !$acc update host(rd)
   do i=1, NP
      call apply(id, q(:,i), DS, rh(i))
   enddo
   do i=1, NP
      mag = ksumabs(id)*maxval(abs(q(klo(id):khi(id),i)))/DS**kn(id)
      n_values = n_values + 1
      if (rd(i) == rh(i)) n_bitwise = n_bitwise + 1
      worst = max(worst, abs(rd(i) - rh(i))/mag)
      if (abs(rd(i) - rh(i)) > 64*epsilon(1._R_P)*mag) then
         n_fail = n_fail + 1
         if (n_fail <= 10) print '(A,A,1X,I0,2ES25.16)', 'DEVICE ', trim(kname(id)), i, rd(i), rh(i)
      endif
   enddo
enddo
!$acc end data
#if defined _OPENACC
print '(A)', 'gideon_test_device: kernels on the device (OpenACC)'
#else
print '(A)', 'gideon_test_device: kernels on the host (no OpenACC in this build)'
#endif
print '(A,I0,A,I0,A,I0,A,ES9.2,A,I0)', 'gideon_test_device: ', NK, ' kernels, ', n_values, ' values, bitwise ', &
                                       n_bitwise, ', worst relative difference ', worst, ', failures ', n_fail
if (n_fail > 0) error stop 'gideon_test_device: failed'
print '(A)', 'gideon_test_device: ok'
endprogram gideon_test_device
