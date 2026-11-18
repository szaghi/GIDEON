!< GIDEON benchmark (#2, L6): a compact operator (right side + line solve) against the explicit one of equal order.
program gideon_bench_compact
!< GIDEON benchmark (#2, L6): a compact operator (right side + line solve) against the explicit one of equal order.
!<
!< 2^24 points: NL = 65536 periodic lines of NN = 256 points, q(i,l) (a line contiguous: the layout the line solvers
!< take, one column per line). The explicit operator and a compact operator's right side: one kernel call per point,
!< `!$acc parallel loop collapse(2)` (one thread per point) on the GPU. The solve: one line per call, one line per
!< thread on the GPU (`!$acc parallel loop`, the scratch a column of the caller's arrays). Method (SETH's, as
!< gideon_bench_adam): 2 warm-up passes, the median of NPASS timed passes (system_clock, `!$acc wait`), ns per point.
!< Frozen budget (L6, written before this program): gideon_fd_d1_p6 (right side + gideon_tri_periodic_const) <= 4 x
!< gideon_fd_d1_c6 on one CPU core, <= 8 x on the GPU. The other cases are measured without a budget.
use penf, only : I4P, I8P, R8P
use gideon
implicit none
integer(I4P), parameter :: NN = 256, NL = 65536, NPASS = 20, G = 4
real(R8P),    parameter :: DS = 0.37_R8P, ALPHA = 0.45_R8P
real(R8P), allocatable  :: q(:,:), r(:,:), w1(:,:), w2(:,:), w3(:,:), w4(:,:)
real(R8P)               :: t(8)
integer(I4P)            :: c, l, i
character(32), parameter :: NAMES(8) = [character(32) :: 'explicit fd_d1_c6', 'compact fd_d1_p6, right side', &
   'compact fd_d1_p6, solve', 'compact fd_d1_p6, both', 'explicit fd_d1_c8', 'compact fd_d1_p8, both', &
   'explicit fd_filter_c8', 'compact fd_filter_p8, both']

allocate(q(1-G:NN+G,NL), r(NN,NL), w1(NN,NL), w2(NN,NL), w3(NN,NL), w4(NN,NL))
call random_number(q)
do l=1, NL ! periodic ghosts
   q(1-G:0,l) = q(NN-G+1:NN,l) ; q(NN+1:NN+G,l) = q(1:G,l)
enddo
!$acc enter data copyin(q) create(r, w1, w2, w3, w4)
#if defined _OPENACC
print '(A)', '# compact against explicit, OpenACC (GPU), ns per point, median of 20 passes, 65536 lines of 256'
#else
print '(A)', '# compact against explicit, one CPU core, ns per point, median of 20 passes, 65536 lines of 256'
#endif
do c=1, 8
   t(c) = bench(c)
enddo
print '(A)', '# case                              ns/point   / explicit of equal order'
do c=1, 8
   select case(c)
   case(2:4) ; print '(A32,F12.4,F12.3)', NAMES(c), t(c), t(c)/t(1)
   case(6)   ; print '(A32,F12.4,F12.3)', NAMES(c), t(c), t(c)/t(5)
   case(8)   ; print '(A32,F12.4,F12.3)', NAMES(c), t(c), t(c)/t(7)
   case default ; print '(A32,F12.4)', NAMES(c), t(c)
   endselect
enddo
!$acc exit data delete(q, r, w1, w2, w3, w4)
contains
   function bench(c) result(tns)
   !< Median ns per point of case c.
   integer(I4P), intent(in) :: c
   real(R8P)                :: tns, tp(NPASS)
   integer(I8P)             :: c0, c1, rate
   integer(I4P)             :: p, k
   do p=-1, NPASS
      call system_clock(c0, rate)
      select case(c)
      case(1)
         !$acc parallel loop collapse(2) present(q, r)
         do l=1, NL
            do i=1, NN
               call gideon_fd_d1_c6(q(i-3:i+3,l), DS, r(i,l))
            enddo
         enddo
      case(2)
         call rhs_p6
      case(3)
         call solve_tri(FD_D1_P6_LEFT(1), FD_D1_P6_LEFT(0))
      case(4)
         call rhs_p6
         call solve_tri(FD_D1_P6_LEFT(1), FD_D1_P6_LEFT(0))
      case(5)
         !$acc parallel loop collapse(2) present(q, r)
         do l=1, NL
            do i=1, NN
               call gideon_fd_d1_c8(q(i-4:i+4,l), DS, r(i,l))
            enddo
         enddo
      case(6)
         !$acc parallel loop collapse(2) present(q, r)
         do l=1, NL
            do i=1, NN
               call gideon_fd_d1_p8(q(i-2:i+2,l), DS, r(i,l))
            enddo
         enddo
         !$acc parallel loop present(r, w1, w2, w3, w4)
         do l=1, NL
            call gideon_penta_periodic_const(NN, FD_D1_P8_LEFT(2), FD_D1_P8_LEFT(1), FD_D1_P8_LEFT(0), &
                                             FD_D1_P8_LEFT(1), FD_D1_P8_LEFT(2), r(:,l), w1(:,l), w2(:,l), w3(:,l), &
                                             w4(:,l))
         enddo
      case(7)
         !$acc parallel loop collapse(2) present(q, r)
         do l=1, NL
            do i=1, NN
               call gideon_fd_filter_c8(q(i-4:i+4,l), r(i,l))
            enddo
         enddo
      case(8)
         !$acc parallel loop collapse(2) present(q, r)
         do l=1, NL
            do i=1, NN
               call gideon_fd_filter_p8(q(i-4:i+4,l), ALPHA, r(i,l))
            enddo
         enddo
         call solve_tri(ALPHA, 1._R8P)
      endselect
      !$acc wait
      call system_clock(c1)
      if (p >= 1) tp(p) = real(c1 - c0, R8P)/rate*1.e9_R8P/(real(NN, R8P)*NL)
   enddo
   do p=2, NPASS ! insertion sort, then the median
      tns = tp(p) ; k = p - 1
      do while (k >= 1)
         if (tp(k) <= tns) exit
         tp(k+1) = tp(k) ; k = k - 1
      enddo
      tp(k+1) = tns
   enddo
   tns = 0.5_R8P*(tp(NPASS/2) + tp(NPASS/2 + 1))
   endfunction bench

   subroutine rhs_p6
   !< The right side of gideon_fd_d1_p6 at every point.
   !$acc parallel loop collapse(2) present(q, r)
   do l=1, NL
      do i=1, NN
         call gideon_fd_d1_p6(q(i-2:i+2,l), DS, r(i,l))
      enddo
   enddo
   endsubroutine rhs_p6

   subroutine solve_tri(a, b)
   !< The periodic tridiagonal solve [a, b, a] of every line, one line per thread.
   real(R8P), intent(in) :: a, b
   !$acc parallel loop present(r, w1, w2)
   do l=1, NL
      call gideon_tri_periodic_const(NN, a, b, a, r(:,l), w1(:,l), w2(:,l))
   enddo
   endsubroutine solve_tri
endprogram gideon_bench_compact
