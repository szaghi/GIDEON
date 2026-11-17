program cookbook_precision
!< Cookbook, precision: the kind the library was built with (R_P), its NaN, and where its round-off stops a
!< convergence. The same source compiles against the R8P library (static-gnu) and the R4P one (static-gnu-r4).
use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
use gideon
implicit none
real(R_P)    :: q(-4:4), h, d, e
integer(I4P) :: k, o

print '(A,I0,A,I0,A,ES9.2)', 'R_P: ', storage_size(1._R_P), ' bits, ', precision(1._R_P), ' digits, epsilon ', &
                              epsilon(1._R_P)
print '(A,L1)', 'GIDEON_NAN is a NaN: ', ieee_is_nan(GIDEON_NAN)
print '(A)', 'gideon_fd_d1_c8 on sin(x) at x = 0.5: error, and its round-off scale eps sum|c| max|q| / h'
do k=1, 12, 2
   h = 0.5_R_P**k
   do o=-4, 4
      q(o) = sin(0.5_R_P + o*h)
   enddo
   call gideon_fd_d1_c8(q, h, d)
   e = abs(d - cos(0.5_R_P))
   ! sum|c| of the order-8 d1 row, both halves of FD_D1_S4: 2 (4/5 + 1/5 + 4/105 + 1/280) = 25/12
   print '(A,I0,T13,2ES10.2)', '  h = 2^-', k, e, epsilon(1._R_P)*(25._R_P/12)*maxval(abs(q))/h
enddo
endprogram cookbook_precision
