program operators_3d
!< Tutorial T5: gradient, Laplacian, divergence and curl of analytic fields on a periodic 3-D grid of point values,
!< from directional stencils gathered out of the 3-D arrays; orders 2, 4, 6 and their observed accuracy.
use gideon
implicit none
integer(I4P), parameter :: G = 3                 ! ghost nodes: the half-width of the order-6 stencils
integer(I4P), parameter :: ORDERS(3) = [2, 4, 6]
integer(I4P), parameter :: NS(3) = [16, 32, 64]  ! nodes per direction
real(R_P),    parameter :: PI = acos(-1._R_P)
real(R_P), allocatable  :: x(:)                  ! node coordinates, the same along x, y, z
real(R_P), allocatable  :: phi(:,:,:)            ! the scalar sin x sin y sin z
real(R_P), allocatable  :: vx(:,:,:), vy(:,:,:), vz(:,:,:) ! the vector (sin x cos y, sin y cos z, sin z cos x)
real(R_P)               :: qx(-G:G), qy(-G:G), qz(-G:G)     ! phi along x, y, z through the node
real(R_P)               :: vx_x(-G:G), vy_y(-G:G), vz_z(-G:G) ! each component along its own axis: the divergence
real(R_P)               :: vy_x(-G:G), vz_x(-G:G), vx_y(-G:G), vz_y(-G:G), vx_z(-G:G), vy_z(-G:G) ! the curl
real(R_P)               :: dxyz(3), grad(3), curl(3), div, lap, h, sx, sy, sz, cx, cy, cz
real(R_P)               :: err(4,3,3)            ! gradient, divergence, curl, Laplacian; order; grid
integer(I4P)            :: l, o, order, n, s, i, j, k

do l=1, 3
   n = NS(l)
   h = 2*PI / n
   dxyz = h
   allocate(x(-G:n-1+G))
   allocate(phi(-G:n-1+G,-G:n-1+G,-G:n-1+G))
   allocate(vx, vy, vz, mold=phi)
   x = [(i*h, i=-G, n-1+G)]
   do k=-G, n-1+G ; do j=-G, n-1+G ; do i=-G, n-1+G  ! the ghosts hold the fields too, as a halo exchange would
      phi(i,j,k) = sin(x(i)) * sin(x(j)) * sin(x(k))
      vx(i,j,k) = sin(x(i)) * cos(x(j))
      vy(i,j,k) = sin(x(j)) * cos(x(k))
      vz(i,j,k) = sin(x(k)) * cos(x(i))
   enddo ; enddo ; enddo
   do o=1, 3
      order = ORDERS(o)
      s = order / 2
      err(:,o,l) = 0._R_P
      do k=0, n-1 ; do j=0, n-1 ; do i=0, n-1
         ! gather: one contiguous stencil per direction (and component), the node (i,j,k) at offset 0
         qx(-s:s) = phi(i-s:i+s,j,k) ; qy(-s:s) = phi(i,j-s:j+s,k) ; qz(-s:s) = phi(i,j,k-s:k+s)
         vx_x(-s:s) = vx(i-s:i+s,j,k) ; vy_y(-s:s) = vy(i,j-s:j+s,k) ; vz_z(-s:s) = vz(i,j,k-s:k+s)
         vy_x(-s:s) = vy(i-s:i+s,j,k) ; vz_x(-s:s) = vz(i-s:i+s,j,k)
         vx_y(-s:s) = vx(i,j-s:j+s,k) ; vz_y(-s:s) = vz(i,j-s:j+s,k)
         vx_z(-s:s) = vx(i,j,k-s:k+s) ; vy_z(-s:s) = vy(i,j,k-s:k+s)
         call gideon_fd_gradient_c(order, qx(-s:s), qy(-s:s), qz(-s:s), dxyz, grad)
         call gideon_fd_laplacian_c(order, qx(-s:s), qy(-s:s), qz(-s:s), dxyz, lap)
         call gideon_fd_divergence_c(order, vx_x(-s:s), vy_y(-s:s), vz_z(-s:s), dxyz, div)
         call gideon_fd_curl_c(order, vy_x(-s:s), vz_x(-s:s), vx_y(-s:s), vz_y(-s:s), vx_z(-s:s), vy_z(-s:s), &
                               dxyz, curl)
         ! the exact values
         sx = sin(x(i)) ; sy = sin(x(j)) ; sz = sin(x(k)) ; cx = cos(x(i)) ; cy = cos(x(j)) ; cz = cos(x(k))
         err(1,o,l) = max(err(1,o,l), maxval(abs(grad - [cx*sy*sz, sx*cy*sz, sx*sy*cz])))
         err(2,o,l) = max(err(2,o,l), abs(div - (cx*cy + cy*cz + cz*cx)))
         err(3,o,l) = max(err(3,o,l), maxval(abs(curl - [sy*sz, sz*sx, sx*sy])))
         err(4,o,l) = max(err(4,o,l), abs(lap + 3*sx*sy*sz))
      enddo ; enddo ; enddo
   enddo
   deallocate(x, phi, vx, vy, vz)
enddo

print '(A)', 'periodic grid of n^3 nodes, point values, max error over the nodes'
print '(A)', 'order     n    gradient  divergence        curl   Laplacian'
do o=1, 3
   do l=1, 3
      print '(I5,I6,4ES12.3)', ORDERS(o), NS(l), err(:,o,l)
   enddo
enddo
print '(A)', 'observed orders, log2(e(n/2) / e(n))'
do o=1, 3
   do l=2, 3
      print '(I5,I6,4F12.2)', ORDERS(o), NS(l), log(err(:,o,l-1) / err(:,o,l)) / log(2._R_P)
   enddo
enddo
endprogram operators_3d
