!< GIDEON benchmark: GIDEON's kernels against ADAM's FDV operators (compiled verbatim), ns per point.
program gideon_bench_adam
!< GIDEON benchmark: GIDEON's kernels against ADAM's FDV operators (compiled verbatim), ns per point.
!<
!< SETH's method: N = 2^24 points of a 1-D array, every kernel applied at every point inside one loop (`!$acc
!< parallel loop` in the OpenACC build, a plain loop on one CPU core otherwise), 2 warm-up passes, then the median of
!< NPASS timed passes (system_clock around the loop and `!$acc wait`). ADAM's operators take a run-time half-stencil
!< s, tables mirrored to the device and assumed-shape windows; GIDEON's a compile-time row and an explicit-shape
!< window. Frozen budget (PERF1, roadmap #1): GIDEON <= 1.0 x ADAM, same device and build. The run-time dispatcher
!< (gideon_fd_d1_c with order 6) is timed against its kernel: budget <= 1.05.
use penf, only : I4P, I8P, R8P
use gideon
use adam_fdv_operators_library
implicit none
integer(I4P), parameter :: N = 2**24, NPASS = 20, G = 8
real(R8P),    parameter :: DS = 0.37_R8P
real(R8P), allocatable  :: q(:), r(:)
real(R8P)               :: t(NPASS), tg, ta
integer(I4P)            :: c, s, k

allocate(q(1-G:N+G), r(N))
call random_number(q)
!$acc enter data copyin(q) create(r)
#if defined _OPENACC
print '(A)', '# GIDEON vs ADAM, OpenACC (GPU), ns per point, median of 20 passes, N = 2^24'
#else
print '(A)', '# GIDEON vs ADAM, one CPU core, ns per point, median of 20 passes, N = 2^24'
#endif
print '(A)', '# case                   GIDEON     ADAM      ratio'
do c=1, 9
   tg = bench(c, .true.) ; ta = bench(c, .false.)
   print '(A24,2F10.4,F10.3)', case_name(c), tg, ta, tg/ta
enddo
tg = bench(10, .true.) ; ta = bench(1 + 2, .true.) ! dispatcher at order 6 against the kernel fd_d1_c6
print '(A24,2F10.4,F10.3,A)', 'dispatch fd_d1_c(6)', tg, ta, tg/ta, '   (against gideon_fd_d1_c6)'
!$acc exit data delete(q, r)
contains
   function case_name(c) result(name)
   integer(I4P), intent(in) :: c
   character(24)            :: name
   character(24), parameter :: NAMES(9) = [character(24) :: 'fd d1 c2', 'fd d1 c4', 'fd d1 c6', 'fd d1 c10', &
                                           'fd d2 c6', 'fv face c6', 'fv d1 c6', 'fv d2 c6', 'fv face ur3']
   name = NAMES(c)
   endfunction case_name

   function bench(c, gid) result(tns)
   !< Median ns per point of case c, GIDEON (gid) or ADAM.
   integer(I4P), intent(in) :: c
   logical,      intent(in) :: gid
   real(R8P)                :: tns
   integer(I8P)             :: c0, c1, rate
   integer(I4P)             :: p, i
   real(R8P)                :: v
   do p=-1, NPASS
      call system_clock(c0, rate)
      select case(c)
      case(1)
         if (gid) then
            !$acc parallel loop present(q, r) private(v)
            do i=1, N ; call gideon_fd_d1_c2(q(i-1:i+1), DS, v) ; r(i) = v ; enddo
         else
            !$acc parallel loop present(q, r) private(v)
            do i=1, N ; call compute_derivative1_fd_centered(s=1, ds=DS, q=q(i-1:i+1), dq_ds=v) ; r(i) = v ; enddo
         endif
      case(2)
         if (gid) then
            !$acc parallel loop present(q, r) private(v)
            do i=1, N ; call gideon_fd_d1_c4(q(i-2:i+2), DS, v) ; r(i) = v ; enddo
         else
            !$acc parallel loop present(q, r) private(v)
            do i=1, N ; call compute_derivative1_fd_centered(s=2, ds=DS, q=q(i-2:i+2), dq_ds=v) ; r(i) = v ; enddo
         endif
      case(3)
         if (gid) then
            !$acc parallel loop present(q, r) private(v)
            do i=1, N ; call gideon_fd_d1_c6(q(i-3:i+3), DS, v) ; r(i) = v ; enddo
         else
            !$acc parallel loop present(q, r) private(v)
            do i=1, N ; call compute_derivative1_fd_centered(s=3, ds=DS, q=q(i-3:i+3), dq_ds=v) ; r(i) = v ; enddo
         endif
      case(4)
         if (gid) then
            !$acc parallel loop present(q, r) private(v)
            do i=1, N ; call gideon_fd_d1_c10(q(i-5:i+5), DS, v) ; r(i) = v ; enddo
         else
            !$acc parallel loop present(q, r) private(v)
            do i=1, N ; call compute_derivative1_fd_centered(s=5, ds=DS, q=q(i-5:i+5), dq_ds=v) ; r(i) = v ; enddo
         endif
      case(5)
         if (gid) then
            !$acc parallel loop present(q, r) private(v)
            do i=1, N ; call gideon_fd_d2_c6(q(i-3:i+3), DS, v) ; r(i) = v ; enddo
         else
            !$acc parallel loop present(q, r) private(v)
            do i=1, N ; call compute_derivative2_fd_centered(s=3, ds=DS, q=q(i-3:i+3), d2q_ds2=v) ; r(i) = v ; enddo
         endif
      case(6)
         if (gid) then
            !$acc parallel loop present(q, r) private(v)
            do i=1, N ; call gideon_fv_face_c6(q(i-2:i+3), v) ; r(i) = v ; enddo
         else
            !$acc parallel loop present(q, r) private(v)
            do i=1, N ; call compute_reconstruction_r_fv_centered(s=3, q=q(i-2:i+3), qr=v) ; r(i) = v ; enddo
         endif
      case(7)
         if (gid) then
            !$acc parallel loop present(q, r) private(v)
            do i=1, N ; call gideon_fv_d1_c6(q(i-3:i+3), DS, v) ; r(i) = v ; enddo
         else
            !$acc parallel loop present(q, r) private(v)
            do i=1, N ; call compute_derivative1_fv_centered(s=3, ds=DS, q=q(i-3:i+3), dq_ds=v) ; r(i) = v ; enddo
         endif
      case(8)
         if (gid) then
            !$acc parallel loop present(q, r) private(v)
            do i=1, N ; call gideon_fv_d2_c6(q(i-3:i+3), DS, v) ; r(i) = v ; enddo
         else
            !$acc parallel loop present(q, r) private(v)
            do i=1, N ; call compute_derivative2_fv_centered(s=3, ds=DS, q=q(i-3:i+3), d2q_ds2=v) ; r(i) = v ; enddo
         endif
      case(9)
         if (gid) then
            !$acc parallel loop present(q, r) private(v)
            do i=1, N ; call gideon_fv_face_ur3(q(i+1:i+3), v) ; r(i) = v ; enddo
         else
            !$acc parallel loop present(q, r) private(v)
            do i=1, N ; call compute_reconstruction_r_fv_rupwind(s=3, q=q(i:i+3), qr=v) ; r(i) = v ; enddo
         endif
      case(10)
         s = 6
         !$acc parallel loop present(q, r) private(v)
         do i=1, N ; call gideon_fd_d1_c(s, q(i-3:i+3), DS, v) ; r(i) = v ; enddo
      endselect
      !$acc wait
      call system_clock(c1)
      if (p >= 1) t(p) = real(c1 - c0, R8P)/real(rate, R8P)*1.e9_R8P/N
   enddo
   do k=1, NPASS - 1 ! insertion sort, then the median
      do i=k + 1, NPASS
         if (t(i) < t(k)) then
            v = t(k) ; t(k) = t(i) ; t(i) = v
         endif
      enddo
   enddo
   tns = 0.5_R8P*(t(NPASS/2) + t(NPASS/2 + 1))
   endfunction bench
endprogram gideon_bench_adam
