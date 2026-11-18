program cookbook_filter
!< Cookbook: filtering a central scheme. Point values of a smooth wave plus an odd-even (Nyquist) mode on a periodic
!< line, filtered by the explicit filters gideon_fd_filter_c<2k> and by the compact ones gideon_fd_filter_p<2k> at
!< alpha_f = 0.45: the odd-even mode is removed, the smooth wave kept to the filter's order.
use gideon
implicit none
integer(I4P), parameter :: G = 5
real(R_P),    parameter :: PI = acos(-1._R_P), ALPHA = 0.45_R_P
real(R_P)               :: e(2,5,2), noise(2,5)
integer(I4P)            :: k, ig, n, kind
character(20)           :: label

do ig=1, 2
   n = 16*ig
   do kind=1, 2
      do k=1, 5
         call filtered(kind, k, n, e(kind,k,ig), noise(kind,k))
      enddo
   enddo
enddo
print '(A)', 'u = sin x + cos(2x)/2 + 0.1 (-1)^i on a periodic line; the filter error on the smooth wave'
print '(A)', 'filter                   n = 16      n = 32   observed   odd-even mode removed'
do kind=1, 2
   do k=1, 5
      write(label, '(A,I0,A,I0)') merge('gideon_fd_filter_c', 'gideon_fd_filter_p', kind == 1), 2*k
      print '(A20,2ES12.3,F11.2,L14)', label, e(kind,k,1), e(kind,k,2), log(e(kind,k,1)/e(kind,k,2))/log(2._R_P), &
         noise(kind,k) < 1.e-12_R_P
   enddo
enddo
contains
   subroutine filtered(kind, k, n, err, nyq)
   !< Filter of order 2k (explicit: kind 1; compact: kind 2) on n points: its error on the smooth part, and what is left
   !< of the odd-even mode (filtering is linear: each part filtered alone).
   integer(I4P), intent(in)  :: kind, k, n
   real(R_P),    intent(out) :: err, nyq
   real(R_P)                 :: x(n), s(1-G:n+G), o(1-G:n+G), fs(n), fo(n), w(n,2)
   integer(I4P)              :: i
   do i=1, n
      x(i) = 2*PI*(i - 1)/n
      s(i) = sin(x(i)) + 0.5_R_P*cos(2*x(i))
      o(i) = 0.1_R_P*(1 - 2*modulo(i, 2))
   enddo
   s(1-G:0) = s(n-G+1:n) ; s(n+1:n+G) = s(1:G)
   o(1-G:0) = o(n-G+1:n) ; o(n+1:n+G) = o(1:G)
   do i=1, n
      call apply(kind, k, s(i-k:i+k), fs(i))
      call apply(kind, k, o(i-k:i+k), fo(i))
   enddo
   if (kind == 2) then ! the compact filter's band [alpha, 1, alpha]
      call gideon_tri_periodic_const(n, ALPHA, 1._R_P, ALPHA, fs, w(:,1), w(:,2))
      call gideon_tri_periodic_const(n, ALPHA, 1._R_P, ALPHA, fo, w(:,1), w(:,2))
   endif
   err = maxval(abs(fs - s(1:n)))
   nyq = maxval(abs(fo))
   endsubroutine filtered

   subroutine apply(kind, k, q, v)
   !< The filter of order 2k at one point (the compact one: its right side).
   integer(I4P), intent(in)  :: kind, k
   real(R_P),    intent(in)  :: q(-k:k)
   real(R_P),    intent(out) :: v
   select case(10*kind + k)
   case(11) ; call gideon_fd_filter_c2(q, v)
   case(12) ; call gideon_fd_filter_c4(q, v)
   case(13) ; call gideon_fd_filter_c6(q, v)
   case(14) ; call gideon_fd_filter_c8(q, v)
   case(15) ; call gideon_fd_filter_c10(q, v)
   case(21) ; call gideon_fd_filter_p2(q, ALPHA, v)
   case(22) ; call gideon_fd_filter_p4(q, ALPHA, v)
   case(23) ; call gideon_fd_filter_p6(q, ALPHA, v)
   case(24) ; call gideon_fd_filter_p8(q, ALPHA, v)
   case(25) ; call gideon_fd_filter_p10(q, ALPHA, v)
   endselect
   endsubroutine apply
endprogram cookbook_filter
