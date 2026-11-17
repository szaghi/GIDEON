!< Cookbook, face derivatives and face values from point values: the viscous flux, upwind-biased interpolation and
!< one-sided extrapolation to a boundary face, each with its error halving the step.
program cookbook_face_from_points
!< Cookbook, face derivatives and face values from point values: the viscous flux, upwind-biased interpolation and
!< one-sided extrapolation to a boundary face, each with its error halving the step.
use gideon
implicit none
real(R_P), parameter :: PI = acos(-1._R_P)
integer              :: k
real(R_P)            :: e1, e2

print '(A)', 'the first derivative at the faces i+1/2 of a periodic grid (the viscous flux), max error on sin(2 pi x)'
do k=1, 3
   e1 = viscous(k, 32) ; e2 = viscous(k, 64)
   print '(A,I0,A,2ES11.3,A,F5.2)', '  gideon_fd_face_d1_c', 2*k, '  N = 32, 64:', e1, e2, '   observed order ', &
                                    log(e1/e2)/log(2._R_P)
enddo
print '(A)', 'the value at the faces, centred and upwind-biased, max error on sin(2 pi x)'
do k=1, 4
   e1 = interpolation(k, 32) ; e2 = interpolation(k, 64)
   print '(A,A,2ES11.3,A,F5.2)', '  ', kernel_i(k), e1, e2, '   observed order ', log(e1/e2)/log(2._R_P)
enddo
print '(A)', 'the value at the boundary face x = 0 from the points inside, exp(x) at x = h/2, 3h/2, ...'
do k=2, 4
   e1 = extrapolation(k, 32) ; e2 = extrapolation(k, 64)
   print '(A,I0,A,2ES11.3,A,F5.2)', '  gideon_fd_face_ur', k, '        N = 32, 64:', e1, e2, '   observed order ', &
                                    log(e1/e2)/log(2._R_P)
enddo
contains
   function viscous(k, n) result(e)
   !< Max error of the face first derivative of order 2k on n points, nodes x_i = i h.
   integer, intent(in) :: k, n
   real(R_P)           :: e, u(-4:n+5), h, d
   integer             :: i
   h = 1._R_P/n
   do i=-4, n+5
      u(i) = sin(2*PI*i*h)
   enddo
   e = 0._R_P
   do i=0, n-1 ! the face between nodes i and i+1: offset 0 is node i, offset 1 node i+1
      select case(k)
      case(1) ; call gideon_fd_face_d1_c2(u(i:i+1), h, d)
      case(2) ; call gideon_fd_face_d1_c4(u(i-1:i+2), h, d)     ! FLUME's (27(W1-W0) - (W2-W-1))/(24 dx)
      case(3) ; call gideon_fd_face_d1_c6(u(i-2:i+3), h, d)
      endselect
      e = max(e, abs(d - 2*PI*cos(2*PI*(i + 0.5_R_P)*h)))
   enddo
   endfunction viscous

   function kernel_i(k) result(name)
   !< The kernel of the k-th interpolation.
   integer, intent(in)     :: k
   character(30)           :: name
   character(30), parameter :: NAMES(4) = [character(30) :: 'gideon_fd_face_c4  N = 32, 64:', &
                                           'gideon_fd_face_bl3 N = 32, 64:', 'gideon_fd_face_bl5 N = 32, 64:', &
                                           'gideon_fd_face_c6  N = 32, 64:']
   name = NAMES(k)
   endfunction kernel_i

   function interpolation(k, n) result(e)
   !< Max error of the k-th face interpolation on n points.
   integer, intent(in) :: k, n
   real(R_P)           :: e, u(-4:n+5), h, f
   integer             :: i
   h = 1._R_P/n
   do i=-4, n+5
      u(i) = sin(2*PI*i*h)
   enddo
   e = 0._R_P
   do i=0, n-1
      select case(k)
      case(1) ; call gideon_fd_face_c4(u(i-1:i+2), f)  ! centred, four points
      case(2) ; call gideon_fd_face_bl3(u(i-1:i+1), f) ! two points left of the face, one right: upwind for a > 0
      case(3) ; call gideon_fd_face_bl5(u(i-2:i+2), f) ! three left, two right
      case(4) ; call gideon_fd_face_c6(u(i-2:i+3), f)  ! centred, six points
      endselect
      e = max(e, abs(f - sin(2*PI*(i + 0.5_R_P)*h)))
   enddo
   endfunction interpolation

   function extrapolation(k, n) result(e)
   !< Error of the one-sided value at x = 0 from the k points x = h/2, 3h/2, ... (offsets 1..k).
   integer, intent(in) :: k, n
   real(R_P)           :: e, u(1:n), h, f
   integer             :: j
   h = 1._R_P/n
   do j=1, n
      u(j) = exp((j - 0.5_R_P)*h)
   enddo
   select case(k)
   case(2) ; call gideon_fd_face_ur2(u(1:2), f)
   case(3) ; call gideon_fd_face_ur3(u(1:3), f)
   case(4) ; call gideon_fd_face_ur4(u(1:4), f)
   endselect
   e = abs(f - 1._R_P)
   endfunction extrapolation
endprogram cookbook_face_from_points
