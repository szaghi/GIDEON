!< GIDEON test (G6): equivalence with ADAM's FDV operators, compiled verbatim from ADAM ff879596 (R8P only).
program gideon_test_adam
!< GIDEON test (G6): equivalence with ADAM's FDV operators, compiled verbatim from ADAM ff879596 (R8P only).
!<
!< adam_fdv_operators_library.F90 and adam_fnl_fdv_operators_library.F90 in src/tests are byte copies of
!< notes/adam-reference/. Every ADAM operator, at every half-stencil s = 1..5, against its GIDEON kernel on random
!< data, ds = 0.37, dxyz = (0.37, 0.53, 0.29). Index maps: ADAM's node or cell sits at q(1) where GIDEON's sits at
!< offset 0 (faces: the same indices); ADAM's "rupwind" d1 reads the right-side cells and its d2-d5 the left-side
!< ones, its "lupwind" the reverse (audit F): GIDEON names the side the stencil lies on.
!<
!< - bitwise: FD and FV centred (face, face derivatives, cell derivatives), deconvolution, the 3-D operators (FNL
!<   directional and CPU cube signatures), the right-side face value, the right-side d1, the left-side d2-d5 and the
!<   left-side face derivatives;
!< - within 64 eps of the scale: the right-side d2-d5 and face derivatives, which ADAM sums farthest cell first;
!< - ADAM's defects, reproduced (not equivalence: GIDEON is exact where ADAM is not): the left-side face value and d1
!<   for s >= 2 miss a linear function, its d1 a quadratic one (audit A); a row ADAM has not got returns 0 (audit C).
use penf, only : I4P, R8P
use gideon
use adam_fdv_operators_library
use adam_fnl_fdv_operators_library
implicit none
real(R8P), parameter :: DS = 0.37_R8P, DXYZ(3) = [0.37_R8P, 0.53_R8P, 0.29_R8P]
real(R8P)    :: q(-12:12), lin(-12:12), cube(-5:5,-5:5,-5:5), vcube(3,-5:5,-5:5,-5:5)
real(R8P)    :: a, g, a3(3), g3(3), tol
integer(I4P) :: s, n, o, fails, bitwise, tolerant, defects, absent
integer, allocatable :: seed(:)
integer      :: ns

call random_seed(size=ns)
allocate(seed(ns)) ; seed = 271828 ; call random_seed(put=seed)
call random_number(q) ; q = 2*q - 1
call random_number(cube) ; cube = 2*cube - 1
call random_number(vcube) ; vcube = 2*vcube - 1
fails = 0 ; bitwise = 0 ; tolerant = 0 ; defects = 0 ; absent = 0

do s=1, 5
   ! FD centred: face interpolation
   call compute_reconstruction_r_fd_centered(s=s, q=q(1-s:s), qr=a)
   call gideon_fd_face_c(2*s, q(1-s:s), g)
   call same('fd face', a, g)
   ! FD centred: node derivatives (ADAM node at q(1))
   do n=1, 6
      a = adam_fd(n, s)
      if (2*s - 2*((n - 1)/2) < 2) then
         call zero('fd d', a)
      else
         call same('fd d', a, gideon_fd(n, 2*s - 2*((n - 1)/2)))
      endif
   enddo
   ! FV centred: face value, face derivatives (ADAM unscaled: GIDEON with ds = 1), cell derivatives
   call compute_reconstruction_r_fv_centered(s=s, q=q(1-s:s), qr=a)
   call gideon_fv_face_c(2*s, q(1-s:s), g)
   call same('fv face', a, g)
   do n=1, 5
      a = adam_fv_face(n, s)
      if (2*s - 2*(n/2) < 2) then
         call zero('fv face d', a)
      else
         call same('fv face d', a, gideon_fv_face(n, 2*s - 2*(n/2)))
      endif
   enddo
   do n=1, 6
      a = adam_fv(n, s)
      if (2*s - 2*((n - 1)/2) < 2) then
         call zero('fv d', a)
      else
         call same('fv d', a, gideon_fv(n, 2*s - 2*((n - 1)/2)))
      endif
   enddo
   ! deconvolution
   call compute_avg_to_point(s=s, q=q(1-s:1+s), qp=a)
   call gideon_fv_point_c(2*s + 2, q(1-s:1+s), g)
   call same('fv point', a, g)
   call compute_point_to_avg(s=s, q=q(1-s:1+s), qa=a)
   call gideon_fd_average_c(2*s + 2, q(1-s:1+s), g)
   call same('fd average', a, g)
   ! upwind face values: right side bitwise; left side: ADAM's double mirror (audit A)
   call compute_reconstruction_r_fv_rupwind(s=s, q=q(0:s), qr=a)
   call gideon_fv_face_ur(s, q(1:s), g)
   call same('fv face ur', a, g)
   do o=-12, 12
      lin(o) = o - 0.5_R8P ! averages of x over the cells [o-1, o]: the face value at 0 is 0
   enddo
   call compute_reconstruction_r_fv_lupwind(s=s, q=lin(1-s:0), qr=a)
   call gideon_fv_face_ul(s, lin(1-s:0), g)
   if (s == 1) then
      call same('fv face ul', a, g)
   else
      call defect('fv face ul (A)', a, g, 0._R8P)
   endif
   ! upwind d1: right side (ADAM window [0:1+s], cell at q(1), q(0) unread) bitwise; left side: audit A
   call compute_derivative1_fv_rupwind(s=s, ds=DS, q=q(0:1+s), dq_ds=a)
   call gideon_fv_d1_ur(s, q(1:1+s), DS, g)
   call same('fv d1 ur', a, g)
   do o=-12, 12
      lin(o) = o**2 + 1._R8P/12 ! averages of x^2 over the cells [o-1/2, o+1/2]: the average of 2x over cell 0 is 0
   enddo                         ! (a linear function hides the defect: both faces err by the same amount)
   call compute_derivative1_fv_lupwind(s=s, ds=1._R8P, q=lin(-s:0), dq_ds=a)
   call gideon_fv_d1_ul(s, lin(-s:0), 1._R8P, g)
   if (s == 1) then
      call same('fv d1 ul', a, g)
   else
      call defect('fv d1 ul (A)', a, g, 0._R8P)
   endif
   ! upwind d2-d5: ADAM's rupwind is GIDEON's left side (bitwise), its lupwind the right side (summation order)
   do n=2, 5
      a = adam_rupwind(n, s)
      if (s + 1 - n < 1) then
         call zero('fv d ul', a)
      else
         call same('fv d ul', a, gideon_fv_ul(n, s + 1 - n))
      endif
      a = adam_lupwind(n, s)
      if (s + 1 - n < 1) then
         call zero('fv d ur', a)
      else
         g = gideon_fv_ur(n, s + 1 - n)
         tol = 64*epsilon(1._R8P)*sum(abs(q(0:s)))*1.e3_R8P/DS**n
         call close('fv d ur', a, g, tol)
      endif
   enddo
   ! upwind face derivatives (ADAM unscaled: ds = 1): rupwind_d is the left side, lupwind_d the right side
   do n=1, 4
      if (s - n < 1) cycle
      a = adam_rupwind_face(n, s)
      call same('fv face d ul', a, gideon_fv_face_u(n, s - n, .false.))
      a = adam_lupwind_face(n, s)
      call close('fv face d ur', a, gideon_fv_face_u(n, s - n, .true.), 64*epsilon(1._R8P)*sum(abs(q(1:s)))*1.e3_R8P)
   enddo
   ! d6 upwind: ADAM returns 0 at every s (S_MAX = 5), GIDEON has it from 6 cells (audit C)
   call compute_derivative6_fv_rupwind(s=s, ds=DS, q=q(0:1+s), d6q_ds6=a)
   call zero('fv d6 rupwind', a)
   ! 3-D operators: FNL directional signature and CPU cube signature
   call three_d(s)
enddo
print '(A,I0,A,I0,A,I0,A,I0,A,I0)', 'gideon_test_adam: bitwise ', bitwise, ', within tolerance ', tolerant, &
                                    ', ADAM defects reproduced ', defects, ', rows ADAM has not (returns 0) ', absent, &
                                    ', failures ', fails
if (fails > 0) error stop 'gideon_test_adam: failed'
print '(A)', 'gideon_test_adam: ok'
contains
   subroutine same(what, a, g)
   character(*), intent(in) :: what
   real(R8P),    intent(in) :: a, g
   if (a /= g) then
      fails = fails + 1
      print '(A,A,A,I0,A,I0,2ES25.17)', 'NOT BITWISE ', what, ' s=', s, ' n=', n, a, g
   else
      bitwise = bitwise + 1
   endif
   endsubroutine same

   subroutine close(what, a, g, tol)
   character(*), intent(in) :: what
   real(R8P),    intent(in) :: a, g, tol
   if (abs(a - g) > tol) then
      fails = fails + 1
      print '(A,A,A,I0,A,I0,3ES12.4)', 'NOT CLOSE ', what, ' s=', s, ' n=', n, a, g, tol
   else
      tolerant = tolerant + 1
   endif
   endsubroutine close

   subroutine zero(what, a)
   character(*), intent(in) :: what
   real(R8P),    intent(in) :: a
   if (a /= 0._R8P) then
      fails = fails + 1
      print '(A,A,A,I0,A,I0,ES12.4)', 'ADAM NOT ZERO ', what, ' s=', s, ' n=', n, a
   else
      absent = absent + 1
   endif
   endsubroutine zero

   subroutine defect(what, a, g, exact)
   !< ADAM misses the exact value of a linear function, GIDEON hits it.
   character(*), intent(in) :: what
   real(R8P),    intent(in) :: a, g, exact
   if (abs(g - exact) > 64*epsilon(1._R8P)*10 .or. abs(a - exact) < 1.e-3_R8P) then
      fails = fails + 1
      print '(A,A,A,I0,3ES12.4)', 'DEFECT NOT AS DOCUMENTED ', what, ' s=', s, a, g, exact
   else
      defects = defects + 1
   endif
   endsubroutine defect

   function adam_fd(n, s) result(v)
   integer(I4P), intent(in) :: n, s
   real(R8P)                :: v
   select case(n)
   case(1) ; call compute_derivative1_fd_centered(s=s, ds=DS, q=q(1-s:1+s), dq_ds=v)
   case(2) ; call compute_derivative2_fd_centered(s=s, ds=DS, q=q(1-s:1+s), d2q_ds2=v)
   case(3) ; call compute_derivative3_fd_centered(s=s, ds=DS, q=q(1-s:1+s), d3q_ds3=v)
   case(4) ; call compute_derivative4_fd_centered(s=s, ds=DS, q=q(1-s:1+s), d4q_ds4=v)
   case(5) ; call compute_derivative5_fd_centered(s=s, ds=DS, q=q(1-s:1+s), d5q_ds5=v)
   case(6) ; call compute_derivative6_fd_centered(s=s, ds=DS, q=q(1-s:1+s), d6q_ds6=v)
   endselect
   endfunction adam_fd

   function gideon_fd(n, p) result(v)
   integer(I4P), intent(in) :: n, p
   real(R8P)                :: v
   integer(I4P)             :: h
   h = p/2 + (n - 1)/2
   select case(n)
   case(1) ; call gideon_fd_d1_c(p, q(1-h:1+h), DS, v)
   case(2) ; call gideon_fd_d2_c(p, q(1-h:1+h), DS, v)
   case(3) ; call gideon_fd_d3_c(p, q(1-h:1+h), DS, v)
   case(4) ; call gideon_fd_d4_c(p, q(1-h:1+h), DS, v)
   case(5) ; call gideon_fd_d5_c(p, q(1-h:1+h), DS, v)
   case(6) ; call gideon_fd_d6_c(p, q(1-h:1+h), DS, v)
   endselect
   endfunction gideon_fd

   function adam_fv_face(n, s) result(v)
   integer(I4P), intent(in) :: n, s
   real(R8P)                :: v
   select case(n)
   case(1) ; call compute_reconstruction_r_fv_centered_d1(s=s, q=q(1-s:s), qr=v)
   case(2) ; call compute_reconstruction_r_fv_centered_d2(s=s, q=q(1-s:s), qr=v)
   case(3) ; call compute_reconstruction_r_fv_centered_d3(s=s, q=q(1-s:s), qr=v)
   case(4) ; call compute_reconstruction_r_fv_centered_d4(s=s, q=q(1-s:s), qr=v)
   case(5) ; call compute_reconstruction_r_fv_centered_d5(s=s, q=q(1-s:s), qr=v)
   endselect
   endfunction adam_fv_face

   function gideon_fv_face(n, p) result(v)
   integer(I4P), intent(in) :: n, p
   real(R8P)                :: v
   integer(I4P)             :: h
   h = p/2 + n/2
   select case(n)
   case(1) ; call gideon_fv_face_d1_c(p, q(1-h:h), 1._R8P, v)
   case(2) ; call gideon_fv_face_d2_c(p, q(1-h:h), 1._R8P, v)
   case(3) ; call gideon_fv_face_d3_c(p, q(1-h:h), 1._R8P, v)
   case(4) ; call gideon_fv_face_d4_c(p, q(1-h:h), 1._R8P, v)
   case(5) ; call gideon_fv_face_d5_c(p, q(1-h:h), 1._R8P, v)
   endselect
   endfunction gideon_fv_face

   function adam_fv(n, s) result(v)
   integer(I4P), intent(in) :: n, s
   real(R8P)                :: v
   select case(n)
   case(1) ; call compute_derivative1_fv_centered(s=s, ds=DS, q=q(1-s:1+s), dq_ds=v)
   case(2) ; call compute_derivative2_fv_centered(s=s, ds=DS, q=q(1-s:1+s), d2q_ds2=v)
   case(3) ; call compute_derivative3_fv_centered(s=s, ds=DS, q=q(1-s:1+s), d3q_ds3=v)
   case(4) ; call compute_derivative4_fv_centered(s=s, ds=DS, q=q(1-s:1+s), d4q_ds4=v)
   case(5) ; call compute_derivative5_fv_centered(s=s, ds=DS, q=q(1-s:1+s), d5q_ds5=v)
   case(6) ; call compute_derivative6_fv_centered(s=s, ds=DS, q=q(1-s:1+s), d6q_ds6=v)
   endselect
   endfunction adam_fv

   function gideon_fv(n, p) result(v)
   integer(I4P), intent(in) :: n, p
   real(R8P)                :: v
   integer(I4P)             :: h
   h = p/2 + (n - 1)/2
   select case(n)
   case(1) ; call gideon_fv_d1_c(p, q(1-h:1+h), DS, v)
   case(2) ; call gideon_fv_d2_c(p, q(1-h:1+h), DS, v)
   case(3) ; call gideon_fv_d3_c(p, q(1-h:1+h), DS, v)
   case(4) ; call gideon_fv_d4_c(p, q(1-h:1+h), DS, v)
   case(5) ; call gideon_fv_d5_c(p, q(1-h:1+h), DS, v)
   case(6) ; call gideon_fv_d6_c(p, q(1-h:1+h), DS, v)
   endselect
   endfunction gideon_fv

   function adam_rupwind(n, s) result(v)
   !< ADAM's rupwind d2-d5: cell at q(1), the left-side cells q(1-s:1).
   integer(I4P), intent(in) :: n, s
   real(R8P)                :: v
   select case(n)
   case(2) ; call compute_derivative2_fv_rupwind(s=s, ds=DS, q=q(1-s:1), d2q_ds2=v)
   case(3) ; call compute_derivative3_fv_rupwind(s=s, ds=DS, q=q(1-s:1), d3q_ds3=v)
   case(4) ; call compute_derivative4_fv_rupwind(s=s, ds=DS, q=q(1-s:1), d4q_ds4=v)
   case(5) ; call compute_derivative5_fv_rupwind(s=s, ds=DS, q=q(1-s:1), d5q_ds5=v)
   endselect
   endfunction adam_rupwind

   function gideon_fv_ul(n, p) result(v)
   !< GIDEON's left-side cell derivative at the cell q(1): window q(1-S:1), S = p + n - 1.
   integer(I4P), intent(in) :: n, p
   real(R8P)                :: v
   integer(I4P)             :: h
   h = p + n - 1
   select case(n)
   case(2) ; call gideon_fv_d2_ul(p, q(1-h:1), DS, v)
   case(3) ; call gideon_fv_d3_ul(p, q(1-h:1), DS, v)
   case(4) ; call gideon_fv_d4_ul(p, q(1-h:1), DS, v)
   case(5) ; call gideon_fv_d5_ul(p, q(1-h:1), DS, v)
   endselect
   endfunction gideon_fv_ul

   function adam_lupwind(n, s) result(v)
   !< ADAM's lupwind d2-d5: cell at q(0), the right-side cells q(0:s).
   integer(I4P), intent(in) :: n, s
   real(R8P)                :: v
   select case(n)
   case(2) ; call compute_derivative2_fv_lupwind(s=s, ds=DS, q=q(0:s), d2q_ds2=v)
   case(3) ; call compute_derivative3_fv_lupwind(s=s, ds=DS, q=q(0:s), d3q_ds3=v)
   case(4) ; call compute_derivative4_fv_lupwind(s=s, ds=DS, q=q(0:s), d4q_ds4=v)
   case(5) ; call compute_derivative5_fv_lupwind(s=s, ds=DS, q=q(0:s), d5q_ds5=v)
   endselect
   endfunction adam_lupwind

   function gideon_fv_ur(n, p) result(v)
   integer(I4P), intent(in) :: n, p
   real(R8P)                :: v
   integer(I4P)             :: h
   h = p + n - 1
   select case(n)
   case(2) ; call gideon_fv_d2_ur(p, q(0:h), DS, v)
   case(3) ; call gideon_fv_d3_ur(p, q(0:h), DS, v)
   case(4) ; call gideon_fv_d4_ur(p, q(0:h), DS, v)
   case(5) ; call gideon_fv_d5_ur(p, q(0:h), DS, v)
   endselect
   endfunction gideon_fv_ur

   function adam_rupwind_face(n, s) result(v)
   !< ADAM's rupwind_d (unscaled): the face right of cell 0 from the left-side cells q(1-s:0), passed ascending.
   integer(I4P), intent(in) :: n, s
   real(R8P)                :: v
   select case(n)
   case(1) ; call compute_reconstruction_r_fv_rupwind_d1(s=s, q=q(1-s:0), qr=v)
   case(2) ; call compute_reconstruction_r_fv_rupwind_d2(s=s, q=q(1-s:0), qr=v)
   case(3) ; call compute_reconstruction_r_fv_rupwind_d3(s=s, q=q(1-s:0), qr=v)
   case(4) ; call compute_reconstruction_r_fv_rupwind_d4(s=s, q=q(1-s:0), qr=v)
   endselect
   endfunction adam_rupwind_face

   function adam_lupwind_face(n, s) result(v)
   !< ADAM's lupwind_d (unscaled): the face left of cell 1 from the right-side cells q(1:s).
   integer(I4P), intent(in) :: n, s
   real(R8P)                :: v
   select case(n)
   case(1) ; call compute_reconstruction_r_fv_lupwind_d1(s=s, q=q(1:s), qr=v)
   case(2) ; call compute_reconstruction_r_fv_lupwind_d2(s=s, q=q(1:s), qr=v)
   case(3) ; call compute_reconstruction_r_fv_lupwind_d3(s=s, q=q(1:s), qr=v)
   case(4) ; call compute_reconstruction_r_fv_lupwind_d4(s=s, q=q(1:s), qr=v)
   endselect
   endfunction adam_lupwind_face

   function gideon_fv_face_u(n, p, right) result(v)
   integer(I4P), intent(in) :: n, p
   logical,      intent(in) :: right
   real(R8P)                :: v
   integer(I4P)             :: h
   h = p + n
   if (right) then
      select case(n)
      case(1) ; call gideon_fv_face_d1_ur(p, q(1:h), 1._R8P, v)
      case(2) ; call gideon_fv_face_d2_ur(p, q(1:h), 1._R8P, v)
      case(3) ; call gideon_fv_face_d3_ur(p, q(1:h), 1._R8P, v)
      case(4) ; call gideon_fv_face_d4_ur(p, q(1:h), 1._R8P, v)
      endselect
   else
      select case(n)
      case(1) ; call gideon_fv_face_d1_ul(p, q(1-h:0), 1._R8P, v)
      case(2) ; call gideon_fv_face_d2_ul(p, q(1-h:0), 1._R8P, v)
      case(3) ; call gideon_fv_face_d3_ul(p, q(1-h:0), 1._R8P, v)
      case(4) ; call gideon_fv_face_d4_ul(p, q(1-h:0), 1._R8P, v)
      endselect
   endif
   endfunction gideon_fv_face_u

   subroutine three_d(s)
   !< The 3-D operators: ADAM's FNL (directional) and CPU (cube) signatures against GIDEON's, bitwise.
   integer(I4P), intent(in) :: s
   real(R8P) :: qx(-s:s), qy(-s:s), qz(-s:s), vxx(-s:s), vyy(-s:s), vzz(-s:s)
   real(R8P) :: vyx(-s:s), vzx(-s:s), vxy(-s:s), vzy(-s:s), vxz(-s:s), vyz(-s:s)
   real(R8P) :: qc(1-s:1+s,1-s:1+s,1-s:1+s), vc(3,1-s:1+s,1-s:1+s,1-s:1+s), b, c
   qx = cube(-s:s,0,0) ; qy = cube(0,-s:s,0) ; qz = cube(0,0,-s:s)
   vxx = vcube(1,-s:s,0,0) ; vyy = vcube(2,0,-s:s,0) ; vzz = vcube(3,0,0,-s:s)
   vyx = vcube(2,-s:s,0,0) ; vzx = vcube(3,-s:s,0,0)
   vxy = vcube(1,0,-s:s,0) ; vzy = vcube(3,0,-s:s,0)
   vxz = vcube(1,0,0,-s:s) ; vyz = vcube(2,0,0,-s:s)
   qc = cube(-s:s,-s:s,-s:s) ; vc = vcube(:,-s:s,-s:s,-s:s) ! ADAM's cube: the point at index 1
   n = 1
   ! FD
   call compute_gradient_fd_centered_dev(s=s, dxyz=DXYZ, qsx=qx, qsy=qy, qsz=qz, gradient=a3)
   call gideon_fd_gradient_c(2*s, qx, qy, qz, DXYZ, g3)
   call same3('fd gradient (FNL)', a3, g3)
   call compute_gradient_fd_centered(s=s, dxyz=DXYZ, q=qc, gradient=a3)
   call same3('fd gradient (CPU)', a3, g3)
   call compute_divergence_fd_centered_dev(s=s, dxyz=DXYZ, qsx=vxx, qsy=vyy, qsz=vzz, divergence=b)
   call gideon_fd_divergence_c(2*s, vxx, vyy, vzz, DXYZ, c)
   call same('fd divergence (FNL)', b, c)
   call compute_divergence_fd_centered(s=s, dxyz=DXYZ, q=vc, divergence=b)
   call same('fd divergence (CPU)', b, c)
   call compute_curl_fd_centered_dev(s=s, dxyz=DXYZ, qsx_y=vyx, qsx_z=vzx, qsy_x=vxy, qsy_z=vzy, qsz_x=vxz, &
                                     qsz_y=vyz, curl=a3)
   call gideon_fd_curl_c(2*s, vyx, vzx, vxy, vzy, vxz, vyz, DXYZ, g3)
   call same3('fd curl (FNL)', a3, g3)
   call compute_curl_fd_centered(s=s, dxyz=DXYZ, q=vc, curl=a3)
   call same3('fd curl (CPU)', a3, g3)
   call compute_laplacian_fd_centered_dev(s=s, dxyz=DXYZ, q=qc, laplacian=b)
   call gideon_fd_laplacian_c(2*s, qx, qy, qz, DXYZ, c)
   call same('fd laplacian (FNL)', b, c)
   call compute_laplacian_fd_centered(s=s, dxyz=DXYZ, q=qc, laplacian=b)
   call same('fd laplacian (CPU)', b, c)
   ! FV
   call compute_gradient_fv_centered_dev(s=s, dxyz=DXYZ, qsx=qx, qsy=qy, qsz=qz, gradient=a3)
   call gideon_fv_gradient_c(2*s, qx, qy, qz, DXYZ, g3)
   call same3('fv gradient (FNL)', a3, g3)
   call compute_gradient_fv_centered(s=s, dxyz=DXYZ, q=qc, gradient=a3)
   call same3('fv gradient (CPU)', a3, g3)
   call compute_divergence_fv_centered_dev(s=s, dxyz=DXYZ, qsx=vxx, qsy=vyy, qsz=vzz, divergence=b)
   call gideon_fv_divergence_c(2*s, vxx, vyy, vzz, DXYZ, c)
   call same('fv divergence (FNL)', b, c)
   call compute_divergence_fv_centered(s=s, dxyz=DXYZ, q=vc, divergence=b)
   call same('fv divergence (CPU)', b, c)
   call compute_curl_fv_centered_dev(s=s, dxyz=DXYZ, qsx_y=vyx, qsx_z=vzx, qsy_x=vxy, qsy_z=vzy, qsz_x=vxz, &
                                     qsz_y=vyz, curl=a3)
   call gideon_fv_curl_c(2*s, vyx, vzx, vxy, vzy, vxz, vyz, DXYZ, g3)
   call same3('fv curl (FNL)', a3, g3)
   call compute_curl_fv_centered(s=s, dxyz=DXYZ, q=vc, curl=a3)
   call same3('fv curl (CPU)', a3, g3)
   call compute_laplacian_fv_centered_dev(s=s, dxyz=DXYZ, qsx=qx, qsy=qy, qsz=qz, laplacian=b)
   call gideon_fv_laplacian_c(2*s, qx, qy, qz, DXYZ, c)
   call same('fv laplacian (FNL)', b, c)
   call compute_laplacian_fv_centered(s=s, dxyz=DXYZ, q=qc, laplacian=b)
   call same('fv laplacian (CPU)', b, c)
   endsubroutine three_d

   subroutine same3(what, a, g)
   character(*), intent(in) :: what
   real(R8P),    intent(in) :: a(3), g(3)
   integer(I4P)             :: i
   do i=1, 3
      call same(what, a(i), g(i))
   enddo
   endsubroutine same3
endprogram gideon_test_adam
