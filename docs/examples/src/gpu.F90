program gpu
!< Tutorial T6: the 3-D Laplacian and gradient of order 6 inside an OpenACC parallel loop. Built with nvfortran -acc the
!< loop runs on the GPU; built without OpenACC the directives are comments and it runs on the host. The output is the
!< same: the errors against the exact values, and the device result against a host loop.
use gideon
implicit none
integer(I4P), parameter :: S = 3                 ! half-width of the order-6 stencils, the ghost nodes
integer(I4P), parameter :: NS(2) = [32, 64]      ! nodes per direction
real(R_P),    parameter :: PI = acos(-1._R_P)
real(R_P), allocatable  :: x(:), phi(:,:,:)      ! node coordinates; sin x sin y sin z, with ghosts
real(R_P), allocatable  :: lap(:,:,:), grad(:,:,:,:)       ! computed in the parallel loop
real(R_P), allocatable  :: lap_h(:,:,:), grad_h(:,:,:,:)   ! the same, by a host loop
real(R_P)               :: qx(-S:S), qy(-S:S), qz(-S:S)    ! the gathered stencils, private to each iteration
real(R_P)               :: dxyz(3), h, elap(2), egrad(2), tol
logical                 :: agree(2)
integer(I4P)            :: l, n, i, j, k

do l=1, 2
   n = NS(l)
   h = 2*PI / n
   dxyz = h
   allocate(x(-S:n-1+S), phi(-S:n-1+S,-S:n-1+S,-S:n-1+S))
   allocate(lap(0:n-1,0:n-1,0:n-1), lap_h(0:n-1,0:n-1,0:n-1), grad(3,0:n-1,0:n-1,0:n-1), grad_h(3,0:n-1,0:n-1,0:n-1))
   x = [(i*h, i=-S, n-1+S)]
   do k=-S, n-1+S ; do j=-S, n-1+S ; do i=-S, n-1+S
      phi(i,j,k) = sin(x(i)) * sin(x(j)) * sin(x(k))
   enddo ; enddo ; enddo

   !$acc data copyin(phi, dxyz) copyout(lap, grad)
   !$acc parallel loop collapse(3) private(qx, qy, qz)
   do k=0, n-1 ; do j=0, n-1 ; do i=0, n-1
      qx = phi(i-S:i+S,j,k) ; qy = phi(i,j-S:j+S,k) ; qz = phi(i,j,k-S:k+S)
      call gideon_fd_laplacian_c6(qx, qy, qz, dxyz, lap(i,j,k))
      call gideon_fd_gradient_c6(qx, qy, qz, dxyz, grad(:,i,j,k))
   enddo ; enddo ; enddo
   !$acc end data

   do k=0, n-1 ; do j=0, n-1 ; do i=0, n-1         ! the host's values, the reference of the device's
      qx = phi(i-S:i+S,j,k) ; qy = phi(i,j-S:j+S,k) ; qz = phi(i,j,k-S:k+S)
      call gideon_fd_laplacian_c6(qx, qy, qz, dxyz, lap_h(i,j,k))
      call gideon_fd_gradient_c6(qx, qy, qz, dxyz, grad_h(:,i,j,k))
   enddo ; enddo ; enddo

   elap(l) = 0._R_P ; egrad(l) = 0._R_P
   do k=0, n-1 ; do j=0, n-1 ; do i=0, n-1
      elap(l) = max(elap(l), abs(lap(i,j,k) + 3*phi(i,j,k)))
      egrad(l) = max(egrad(l), maxval(abs(grad(:,i,j,k) - [cos(x(i))*sin(x(j))*sin(x(k)), &
                                                            sin(x(i))*cos(x(j))*sin(x(k)), &
                                                            sin(x(i))*sin(x(j))*cos(x(k))])))
   enddo ; enddo ; enddo
   ! device and host code are compiled apart, and a compiler may contract their multiply-adds differently: compare
   ! within a few ulps of the operator's scale, not bitwise
   tol = 64 * epsilon(1._R_P)
   agree(l) = maxval(abs(lap - lap_h)) <= tol * maxval(abs(lap_h)) / h**2 .and. &
              maxval(abs(grad - grad_h)) <= tol * maxval(abs(grad_h)) / h
   deallocate(x, phi, lap, lap_h, grad, grad_h)
enddo

print '(A)', 'order 6, sin x sin y sin z on a periodic grid of n^3 nodes, max error'
print '(A)', '     n   Laplacian    gradient   parallel loop == host loop (64 eps)'
do l=1, 2
   print '(I6,2ES12.3,L12)', NS(l), elap(l), egrad(l), agree(l)
enddo
print '(A,2F6.2)', 'observed orders, log2(e(32) / e(64)):', log(elap(1) / elap(2)) / log(2._R_P), &
                                                           log(egrad(1) / egrad(2)) / log(2._R_P)
endprogram gpu
