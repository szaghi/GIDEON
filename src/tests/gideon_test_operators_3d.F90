!< GIDEON test: the 3-D operators, FD and FV, every order through the run-time dispatchers, exact on polynomials.
program gideon_test_operators_3d
!< GIDEON test: the 3-D operators, FD and FV, every order through the run-time dispatchers, exact on polynomials.
!<
!< Random polynomial fields of total degree p (the order): f for the gradient and the Laplacian, v = (f1, f2, f3) for
!< the divergence and the curl, around (0.3, -0.2, 0.1) with the unequal steps (0.5, 0.25, 0.75). FD data are point
!< values along each axis; FV data are averages along the axis of each stencil (the 1-D cells the FV operators
!< difference), and the FV targets are the averages over that cell of the derivatives. Both are exact for degree <= p
!< (d1 of order p, d2 of order p + 1): tolerance 1e4 eps (1 + max|data|) / min(h)^n, frozen. An order without a kernel
!< returns NaN.
use penf, only : I4P, R_P
use gideon
implicit none
integer(I4P), parameter :: PMAX = 10
real(R_P),    parameter :: X0(3) = [0.3_R_P, -0.2_R_P, 0.1_R_P], H(3) = [0.5_R_P, 0.25_R_P, 0.75_R_P]
real(R_P) :: cf(0:PMAX,0:PMAX,0:PMAX,3)
real(R_P) :: qx(-5:5), qy(-5:5), qz(-5:5), vx(-5:5,3), vy(-5:5,3), vz(-5:5,3)
real(R_P) :: g(3), ge(3), c(3), ce(3), div, dive, lap, lape, tol
integer(I4P) :: p, s, fails, fv, comp
character(2) :: data

fails = 0
do fv=0, 1
   data = merge('fv', 'fd', fv == 1)
   do p=2, PMAX, 2
      s = p/2
      call random_polynomials(p)
      ! stencils of f (component 1) and of every component of v along each axis
      call line(1, 1, qx(-s:s), s, fv) ; call line(1, 2, qy(-s:s), s, fv) ; call line(1, 3, qz(-s:s), s, fv)
      do comp=1, 3
         call line(comp, 1, vx(-s:s,comp), s, fv)
         call line(comp, 2, vy(-s:s,comp), s, fv)
         call line(comp, 3, vz(-s:s,comp), s, fv)
      enddo
      tol = 1.e4_R_P*epsilon(1._R_P)*(1 + max(maxval(abs(vx(-s:s,:))), maxval(abs(vy(-s:s,:))), &
                                              maxval(abs(vz(-s:s,:)))))
      ! exact targets
      ge = [dmean(1, 1, 1, fv), dmean(1, 2, 1, fv), dmean(1, 3, 1, fv)]
      dive = dmean(1, 1, 1, fv) + dmean(2, 2, 1, fv) + dmean(3, 3, 1, fv)
      ce = [dmean(3, 2, 1, fv) - dmean(2, 3, 1, fv), dmean(1, 3, 1, fv) - dmean(3, 1, 1, fv), &
            dmean(2, 1, 1, fv) - dmean(1, 2, 1, fv)]
      lape = dmean(1, 1, 2, fv) + dmean(1, 2, 2, fv) + dmean(1, 3, 2, fv)
      if (fv == 0) then
         call gideon_fd_gradient_c(p, qx(-s:s), qy(-s:s), qz(-s:s), H, g)
         call gideon_fd_divergence_c(p, vx(-s:s,1), vy(-s:s,2), vz(-s:s,3), H, div)
         call gideon_fd_curl_c(p, vx(-s:s,2), vx(-s:s,3), vy(-s:s,1), vy(-s:s,3), vz(-s:s,1), vz(-s:s,2), H, c)
         call gideon_fd_laplacian_c(p, qx(-s:s), qy(-s:s), qz(-s:s), H, lap)
      else
         call gideon_fv_gradient_c(p, qx(-s:s), qy(-s:s), qz(-s:s), H, g)
         call gideon_fv_divergence_c(p, vx(-s:s,1), vy(-s:s,2), vz(-s:s,3), H, div)
         call gideon_fv_curl_c(p, vx(-s:s,2), vx(-s:s,3), vy(-s:s,1), vy(-s:s,3), vz(-s:s,1), vz(-s:s,2), H, c)
         call gideon_fv_laplacian_c(p, qx(-s:s), qy(-s:s), qz(-s:s), H, lap)
      endif
      call check('gradient', maxval(abs(g - ge)), tol/minval(H))
      call check('divergence', abs(div - dive), tol/minval(H))
      call check('curl', maxval(abs(c - ce)), tol/minval(H))
      call check('laplacian', abs(lap - lape), tol/minval(H)**2)
   enddo
   ! an order without a kernel: NaN (order 3: window -1:1)
   if (fv == 0) then
      call gideon_fd_gradient_c(3, qx(-1:1), qy(-1:1), qz(-1:1), H, g)
      call gideon_fd_laplacian_c(3, qx(-1:1), qy(-1:1), qz(-1:1), H, lap)
   else
      call gideon_fv_gradient_c(3, qx(-1:1), qy(-1:1), qz(-1:1), H, g)
      call gideon_fv_laplacian_c(3, qx(-1:1), qy(-1:1), qz(-1:1), H, lap)
   endif
   if (all(g == g) .or. lap == lap) then
      fails = fails + 1
      print '(A,A)', 'NOT NAN ', data
   endif
enddo
print '(A,I0)', 'gideon_test_operators_3d: 4 operators x 2 data x 5 orders, failures ', fails
if (fails > 0) error stop 'gideon_test_operators_3d: failed'
print '(A)', 'gideon_test_operators_3d: ok'
contains
   subroutine check(what, err, t)
   character(*), intent(in) :: what
   real(R_P),    intent(in) :: err, t
   if (.not. err <= t) then
      fails = fails + 1
      print '(A,A,A,A,I0,A,ES10.2,A,ES10.2)', data, ' ', what, ' order ', p, ': error ', err, ' > tol ', t
   endif
   endsubroutine check

   subroutine random_polynomials(p)
   !< Three random polynomials of total degree p, coefficients in [-1, 1].
   integer(I4P), intent(in) :: p
   integer(I4P)             :: a, b, cc
   integer, allocatable     :: seed(:)
   integer                  :: n
   call random_seed(size=n)
   allocate(seed(n)) ; seed = 31 + p ; call random_seed(put=seed)
   call random_number(cf)
   cf = 2*cf - 1
   do a=0, PMAX ; do b=0, PMAX ; do cc=0, PMAX
      if (a + b + cc > p) cf(a,b,cc,:) = 0._R_P
   enddo ; enddo ; enddo
   endsubroutine random_polynomials

   pure function pmono(j, k, x) result(v)
   !< k-th derivative of x^j (k = -1: the primitive x^(j+1)/(j+1)).
   integer(I4P), intent(in) :: j, k
   real(R_P),    intent(in) :: x
   real(R_P)                :: v
   integer(I4P)             :: i
   if (k == -1) then
      v = x**(j + 1)/(j + 1)
      return
   endif
   v = 0._R_P
   if (j < k) return
   v = 1._R_P
   do i=j - k + 1, j
      v = v*i
   enddo
   v = v*x**(j - k)
   endfunction pmono

   pure function poly(comp, xyz, dir, k) result(v)
   !< Component comp at xyz, differentiated k times (-1: integrated once) along the axis dir.
   integer(I4P), intent(in) :: comp, dir, k
   real(R_P),    intent(in) :: xyz(3)
   real(R_P)                :: v
   integer(I4P)             :: e(3), a, b, cc, ks(3)
   v = 0._R_P
   ks = 0 ; ks(dir) = k
   do a=0, PMAX ; do b=0, PMAX ; do cc=0, PMAX
      if (cf(a,b,cc,comp) == 0._R_P) cycle
      e = [a, b, cc]
      v = v + cf(a,b,cc,comp)*pmono(e(1), ks(1), xyz(1))*pmono(e(2), ks(2), xyz(2))*pmono(e(3), ks(3), xyz(3))
   enddo ; enddo ; enddo
   endfunction poly

   subroutine line(comp, dir, q, s, fv)
   !< The stencil of component comp along the axis dir: point values (fv = 0) or averages over the 1-D cells (fv = 1).
   integer(I4P), intent(in)  :: comp, dir, s, fv
   real(R_P),    intent(out) :: q(-s:s)
   real(R_P)                 :: xa(3), xb(3)
   integer(I4P)              :: o
   do o=-s, s
      xa = X0 ; xb = X0
      if (fv == 0) then
         xa(dir) = X0(dir) + o*H(dir)
         q(o) = poly(comp, xa, dir, 0)
      else
         xa(dir) = X0(dir) + (o - 0.5_R_P)*H(dir) ; xb(dir) = X0(dir) + (o + 0.5_R_P)*H(dir)
         q(o) = (poly(comp, xb, dir, -1) - poly(comp, xa, dir, -1))/H(dir)
      endif
   enddo
   endsubroutine line

   function dmean(comp, dir, k, fv) result(v)
   !< The k-th derivative of component comp along dir: at X0 (fv = 0) or averaged over the 1-D cell of X0 (fv = 1).
   integer(I4P), intent(in) :: comp, dir, k, fv
   real(R_P)                :: v, xa(3), xb(3)
   if (fv == 0) then
      v = poly(comp, X0, dir, k)
   else
      xa = X0 ; xb = X0
      xa(dir) = X0(dir) - 0.5_R_P*H(dir) ; xb(dir) = X0(dir) + 0.5_R_P*H(dir)
      v = (poly(comp, xb, dir, k - 1) - poly(comp, xa, dir, k - 1))/H(dir)
   endif
   endfunction dmean
endprogram gideon_test_operators_3d
