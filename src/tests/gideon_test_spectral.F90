!< GIDEON test (A3): the modified wavenumber of the compiled first- and second-derivative kernels.
program gideon_test_spectral
!< GIDEON test (A3): the modified wavenumber of the compiled first- and second-derivative kernels.
!<
!< A derivative kernel applied to the Fourier mode e^(i o xi) of its window (h = 1) returns its symbol
!< sigma(xi) = sum_o w_o e^(i o xi), from two real calls, on cos(o xi) and sin(o xi): the kernels are linear. The
!< same holds for cell averages (the averaging factor of the data and of the target cancel), so node (fd) and cell
!< (fv) operators are read the same way. First derivative: the modified wavenumber xi~ = -i sigma, exact xi~ = xi;
!< Re xi~ is the dispersion, Im xi~ the dissipation (Im xi~ > 0 amplifies a mode of u_t + u_x = 0). Second
!< derivative: xi~^2 = -sigma, exact xi^2. xi = j pi/256, j = 1..256.
!<
!< Frozen bounds (written before the kernels were run), tolerance 64 eps sum_k |c_k|:
!< - centred kernels (d1, d2): Im = 0 within the tolerance at every xi (the pair form: no dissipation);
!< - the left-side upwind d1 of 1 and 2 cells: Im <= tolerance (dissipative for a flow to the right);
!< - every left-side kernel is the mirror of its right twin: sigma_ul = (-1)^n conj(sigma_ur), within the tolerance.
!< The left-side d1 that amplify some band are listed (a property of the exact rows, documented in A3, not a failure).
!<
!< Report mode (scripts/run_report.sh): with GIDEON_REPORT_DIR set, every symbol is written to spectral-rNN.csv
!< (kernel, xi, Re, Im of xi~ for d1, of xi~^2 for d2), NN the bits of the working kind.
use penf, only : I4P, R8P, R_P
use gideon_test_registry
implicit none
integer(I4P), parameter :: NXI = 256
real(R_P),    parameter :: PI = 3.14159265358979323846264338327950288_R_P
real(R_P)    :: qc(-WMAX:WMAX), qs(-WMAX:WMAX), xi, rc, rs, tol, re(NXI), im(NXI), rem(NXI), imm(NXI)
integer(I4P) :: id, j, o, fails, checked, amplifying
character(4096) :: report_dir
integer      :: rlen, urep
logical      :: report

call get_environment_variable('GIDEON_REPORT_DIR', report_dir, length=rlen)
report = rlen > 0
if (report) then
   open(newunit=urep, file=trim(report_dir)//'/spectral-r'//kind_tag()//'.csv', action='write', status='replace')
   write(urep, '(A)') 'kernel,xi,re,im'
endif
fails = 0 ; checked = 0 ; amplifying = 0
do id=1, NK
   if (.not. derivative(id)) cycle
   checked = checked + 1
   tol = 64*epsilon(1._R_P)*ksumabs(id)
   call symbol(id, re, im)
   if (report) then
      do j=1, NXI
         write(urep, '(A,",",ES24.16E3,",",ES24.16E3,",",ES24.16E3)') trim(kname(id)), real(j*PI/NXI, R8P), &
                                                                     real(re(j), R8P), real(im(j), R8P)
      enddo
   endif
   if (kbias(id) == 'c ') then
      if (any(abs(im) > tol)) then
         fails = fails + 1
         print '(A,A,ES10.2)', 'CENTRED DISSIPATES ', trim(kname(id)), maxval(abs(im))
      endif
   endif
   if (kbias(id) == 'ul' .and. kn(id) == 1) then
      if (maxval(im) > tol) then
         amplifying = amplifying + 1
         print '(A,A,A,ES10.3,A,F6.3,A,F6.3,A)', 'amplifying band: ', trim(kname(id)), ', max Im ', maxval(im), &
                                                 ' (xi/pi in ', band(im, tol, .true.), ', ', band(im, tol, .false.), ')'
         if (korder(id) <= 2) then
            fails = fails + 1
            print '(A,A)', 'NOT DISSIPATIVE ', trim(kname(id))
         endif
      endif
   endif
   if (kmirror(id) /= 0) then
      call symbol(kmirror(id), rem, imm)
      ! sigma_ul = (-1)^n conj(sigma_ur): for d1 (xi~ = -i sigma) Re equal, Im opposite; for d2 (-sigma) the same
      if (any(abs(re - rem) > tol) .or. any(abs(im + imm) > tol)) then
         fails = fails + 1
         print '(A,A)', 'MIRROR SYMBOL ', trim(kname(id))
      endif
   endif
enddo
if (report) close(urep)
print '(A,I0,A,I0,A,I0)', 'gideon_test_spectral: ', checked, &
                          ' d1 and d2 kernels, left-side d1 with an amplifying band ', amplifying, ', failures ', fails
if (fails > 0) error stop 'gideon_test_spectral: failed'
print '(A)', 'gideon_test_spectral: ok'
contains
   pure function derivative(id) result(yes)
   !< A first or second derivative at a node (fd) or averaged over a cell (fv): the operators of a modified wavenumber.
   integer(I4P), intent(in) :: id
   logical                  :: yes
   yes = (kgeom(id) == 1 .or. kgeom(id) == 4) .and. (kn(id) == 1 .or. kn(id) == 2)
   endfunction derivative

   subroutine symbol(id, re, im)
   !< xi~ (d1) or xi~^2 (d2) of kernel id at xi = j pi/NXI, from the compiled kernel.
   integer(I4P), intent(in)  :: id
   real(R_P),    intent(out) :: re(NXI), im(NXI)
   do j=1, NXI
      xi = j*PI/NXI
      do o=-WMAX, WMAX
         qc(o) = cos(o*xi) ; qs(o) = sin(o*xi)
      enddo
      call apply(id, qc, 1._R_P, rc)
      call apply(id, qs, 1._R_P, rs)
      if (kn(id) == 1) then ! xi~ = -i (rc + i rs)
         re(j) = rs ; im(j) = -rc
      else                  ! xi~^2 = -(rc + i rs)
         re(j) = -rc ; im(j) = -rs
      endif
   enddo
   endsubroutine symbol

   function band(im, tol, lower) result(b)
   !< The lower or upper end (xi/pi) of the wavenumbers where Im > tol.
   real(R_P),    intent(in) :: im(NXI), tol
   logical,      intent(in) :: lower
   real(R_P)                :: b
   integer(I4P)             :: jj
   b = 0._R_P
   do jj=1, NXI
      if (im(jj) > tol) then
         b = real(jj, R_P)/NXI
         if (lower) return
      endif
   enddo
   endfunction band

   function kind_tag() result(tag)
   !< The bits of the working kind, the suffix of the report file: 32, 64 or 128.
   character(:), allocatable :: tag
   character(8)              :: buf
   write(buf, '(I0)') storage_size(1._R_P)
   tag = trim(buf)
   endfunction kind_tag
endprogram gideon_test_spectral
