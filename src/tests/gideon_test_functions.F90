!< GIDEON tests, the data and the exact targets: monomials (exactness) and a generic smooth function (orders).
module gideon_test_functions
!< GIDEON tests, the data and the exact targets: monomials (exactness) and a generic smooth function (orders).
!<
!< The generic function is g(x) = sin(A x + P) + B cos(C x + Q), two incommensurate waves at a generic phase: no
!< derivative vanishes at the test point, so no leading error term cancels by accident (a single sine at a symmetric
!< point, or exp, would hide an order). Its k-th derivative is closed-form for every k >= -1 (k = -1: a primitive).
use penf, only : I4P, R_P
use gideon_test_registry, only : WMAX, kgeom, ktarget, kn, klo, khi
implicit none
private
public :: GX0, data_monomial, target_monomial, data_generic, target_generic, gk

real(R_P), parameter :: A = 1.3_R_P, P = 0.4_R_P, B = 0.5_R_P, C = 0.7_R_P, Q = -0.2_R_P
real(R_P), parameter :: GX0 = 0.37_R_P !< The test point of the generic function.
real(R_P), parameter :: HALF_PI = 1.5707963267948966192313216916397514_R_P

contains
   pure function gk(k, x) result(v)
   !< k-th derivative of g at x (k = -1: the primitive -cos(A x + P)/A + B sin(C x + Q)/C).
   integer(I4P), intent(in) :: k !< Derivative order, >= -1.
   real(R_P),    intent(in) :: x !< Abscissa.
   real(R_P)                :: v !< Value.
   v = A**k*sin(A*x + P + k*HALF_PI) + B*C**k*cos(C*x + Q + k*HALF_PI)
   endfunction gk

   pure subroutine window_cell(id, o, x0, h, a, b, xp)
   !< The datum at offset o: the cell [a, b] (geometries 3, 4) or the point xp (geometries 1, 2).
   integer(I4P), intent(in)  :: id, o
   real(R_P),    intent(in)  :: x0, h
   real(R_P),    intent(out) :: a, b, xp
   a = 0._R_P ; b = 0._R_P ; xp = 0._R_P
   select case(kgeom(id))
   case(1) ; xp = x0 + o*h
   case(2) ; xp = x0 + (o - 0.5_R_P)*h
   case(3) ; a = x0 + (o - 1)*h ; b = x0 + o*h
   case(4) ; a = x0 + (o - 0.5_R_P)*h ; b = x0 + (o + 0.5_R_P)*h
   endselect
   endsubroutine window_cell

   pure subroutine data_monomial(id, j, q)
   !< Data of x^j (h = 1, x0 = 0) on the window of kernel id.
   integer(I4P), intent(in)  :: id, j
   real(R_P),    intent(out) :: q(-WMAX:WMAX)
   real(R_P)                 :: a, b, xp
   integer(I4P)              :: o
   q = 0._R_P
   do o=klo(id), khi(id)
      call window_cell(id, o, 0._R_P, 1._R_P, a, b, xp)
      if (kgeom(id) <= 2) then
         q(o) = xp**j
      else
         q(o) = (b**(j + 1) - a**(j + 1))/(j + 1)
      endif
   enddo
   endsubroutine data_monomial

   pure function target_monomial(id, j) result(t)
   !< Exact target of kernel id on x^j (h = 1, x0 = 0).
   integer(I4P), intent(in) :: id, j
   real(R_P)                :: t
   integer(I4P)             :: n, i
   real(R_P)                :: f
   n = kn(id)
   t = 0._R_P
   if (j < n) return
   f = 1._R_P
   do i=j - n + 1, j
      f = f*i
   enddo
   if (ktarget(id) == 1) then
      if (j == n) t = f
   else ! average over [-1/2, 1/2] of f x^(j-n)
      t = f*(0.5_R_P**(j - n + 1) - (-0.5_R_P)**(j - n + 1))/(j - n + 1)
   endif
   endfunction target_monomial

   pure subroutine data_generic(id, h, q)
   !< Data of g at step h around GX0 on the window of kernel id.
   integer(I4P), intent(in)  :: id
   real(R_P),    intent(in)  :: h
   real(R_P),    intent(out) :: q(-WMAX:WMAX)
   real(R_P)                 :: a, b, xp
   integer(I4P)              :: o
   q = 0._R_P
   do o=klo(id), khi(id)
      call window_cell(id, o, GX0, h, a, b, xp)
      if (kgeom(id) <= 2) then
         q(o) = gk(0, xp)
      else
         q(o) = (gk(-1, b) - gk(-1, a))/h
      endif
   enddo
   endsubroutine data_generic

   pure function target_generic(id, h) result(t)
   !< Exact target of kernel id on g at step h around GX0.
   integer(I4P), intent(in) :: id
   real(R_P),    intent(in) :: h
   real(R_P)                :: t
   if (ktarget(id) == 1) then
      t = gk(kn(id), GX0)
   else
      t = (gk(kn(id) - 1, GX0 + 0.5_R_P*h) - gk(kn(id) - 1, GX0 - 0.5_R_P*h))/h
   endif
   endfunction target_generic
endmodule gideon_test_functions
