subroutine fixaray
   use mod_konst
   use mod_curve3d
   implicit none

   !*******************************************************************
   !    This routine stores toroidal and poloidal mode number arrays
   !
   !    MPOL = NUMBER OF POLOIDAL MODES USED IN CURVE-FIT
   !    NPHI = NUMBER OF EQUALLY-SPACED PHI PLANES PER FIELD
   !           PERIOD IN WHICH THE CURVE-FITTING IS PERFORMED
   !    IMPORTANT: NPHI MUST BE EVEN FOR NON-SYMMETRIC SYSTEMS
   !
   !    MPNT = NUMBER OF R AND Z MODES APPEARING IN FIT
   !    NTHETA=NUMBER OF THETA POINTS IN EACH PHI PLANE
   !    N2   = TOTAL NUMBER OF RESIDUALS IN FSQ (PER PHI PLANE)
   !    NFP  = NUMBER OF TOROIDAL FIELD PERIODS (IRRELEVANT)
   !*******************************************************************

   integer :: l, n, m, ntor, nn0

   l = 0
   ntor = max0(1,nphi-1)
   nn0 = 1 - (ntor+1)/2
   do n=1,ntor
      nn(n) = (nn0 + (n-1))*nfp
   end do
   do m = 1,mpol
      mm(m) = m-1
      do n = 1,ntor
         if( mm(m).eq.0 .and. nn(n).lt.0 )  cycle
         l = l+1
         m1(l) = mm(m)
         n1(l) = nn(n)
      end do
   end do
   mpnt = l
   dnorm = 2./ntheta


   do m = 1,mpol
      dm1(m) = real (m-1)
      faccon(m) = .125*dnorm/(1.+dm1(m))**pexp
      if( m > 1 ) then
         xmpq(m,1) = dm1(m)**pexp
         xmpq(m,2) = dm1(m)**qexp
         xmpq(m,4) = xmpq(m,1)
         if ( m.le.2 ) xmpq(m,4) = 0.
         xmpq(m,3) = sqrt(xmpq(m,4))
      end if
   end do

end subroutine fixaray
