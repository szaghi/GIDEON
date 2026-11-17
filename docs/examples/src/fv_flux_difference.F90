program fv_flux_difference
!< Tutorial T2: cell averages of u = sin x + cos(3x)/2 on a periodic grid; the sixth-order face values, the cell-average
!< first derivative as their difference, and the telescoping sum that makes a finite-volume update conservative.
use gideon
implicit none
integer(I4P), parameter :: G = 3                 ! ghost cells on each side: the half-width of the c6 cell window
integer(I4P), parameter :: NGRID = 4             ! number of grids: 16, 32, 64, 128 cells
real(R_P),    parameter :: PI = acos(-1._R_P)
real(R_P), allocatable  :: xf(:), q(:), f(:), d(:) ! faces, cell averages (with ghosts), face values, cell derivatives
real(R_P)               :: h, ef(NGRID), ed(NGRID), telescoped
integer(I4P)            :: n, i, k, nbitwise

print '(A)', 'cell averages of sin x + cos(3x)/2 on [0, 2 pi), order 6'
print '(A)', '     n  face error    d1 error   d1 == flux difference   telescoped to round-off'
do k=1, NGRID
   n = 8 * 2**k
   h = 2*PI / n
   allocate(xf(-1:n-1), q(-G:n-1+G), f(-1:n-1), d(0:n-1))
   xf = [(i*h + h, i=-1, n-1)]                    ! xf(i): the face i+1/2, right of cell i = [i h, (i+1) h]
   q(0:n-1) = (prim(xf(0:n-1)) - prim(xf(-1:n-2))) / h ! exact cell averages
   q(-G:-1) = q(n-G:n-1)                          ! periodic ghost cells
   q(n:n-1+G) = q(0:G-1)
   do i=0, n-1
      call gideon_fv_face_c6(q(i-2:i+3), f(i))   ! window q(-2:3): cells i-2..i+3, the face between i and i+1
      call gideon_fv_d1_c6(q(i-3:i+3), h, d(i))  ! window q(-3:3): the derivative averaged over cell i
   enddo
   f(-1) = f(n-1)                                 ! the face 0 - 1/2 is the face n - 1/2
   ef(k) = maxval(abs(f(0:n-1) - fun(xf(0:n-1))))
   ed(k) = maxval(abs(d - (fun(xf(0:n-1)) - fun(xf(-1:n-2))) / h))
   nbitwise = count(d == (f(0:n-1) - f(-1:n-2)) / h)
   telescoped = sum(f(0:n-1) - f(-1:n-2))        ! every face value enters twice, with opposite signs
   print '(I6,2ES12.3,6X,I4,A,I4,A,15X,L1)', n, ef(k), ed(k), nbitwise, ' of ', n, ' cells', &
         abs(telescoped) <= n * epsilon(1._R_P) * maxval(abs(f))
   deallocate(xf, q, f, d)
enddo
print '(A)', 'observed orders, log2(e(n/2) / e(n))'
do k=2, NGRID
   print '(I6,2F12.2)', 8 * 2**k, log(ef(k-1) / ef(k)) / log(2._R_P), log(ed(k-1) / ed(k)) / log(2._R_P)
enddo
contains
   elemental function fun(x)
   !< The function.
   real(R_P), intent(in) :: x
   real(R_P)             :: fun
   fun = sin(x) + cos(3*x) / 2
   endfunction fun

   elemental function prim(x)
   !< A primitive of the function: its cell averages are differences of prim.
   real(R_P), intent(in) :: x
   real(R_P)             :: prim
   prim = -cos(x) + sin(3*x) / 6
   endfunction prim
endprogram fv_flux_difference
