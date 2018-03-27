subroutine fftrans (rmc,rms,zmc,zms,rbc,zbs,rbs,zbc,rmnaxis,zmnaxis)
   use mod_konst
   use mod_curve3d
   implicit none

   double precision, intent(in)  :: rmc(nphi2,mpol), rms(nphi2,mpol)
   double precision, intent(in)  :: zmc(nphi2,mpol), zms(nphi2,mpol)
   double precision, intent(out) :: rbc(0:mpol1,-nphi2:nphi2)
   double precision, intent(out) :: zbs(0:mpol1,-nphi2:nphi2)
   double precision, intent(out) :: rbs(0:mpol1,-nphi2:nphi2)
   double precision, intent(out) :: zbc(0:mpol1,-nphi2:nphi2)
   double precision, intent(out) :: rmnaxis(0:nphi2), zmnaxis(0:nphi2)

   integer          :: i, mn, m, mreal, nreal
   double precision :: intgrate(nphi2), argi(nphi2)
   double precision :: dn, tcosn, tsinn, delphi, arg

   !*******************************************************************
   !    PERFORM FOURIER TRANSFORM IN phi
   !*******************************************************************

   delphi = 2./ real (nphi)
   do i = 1,nphi2
      intgrate(i) = delphi
      argi(i) = twopi* real (i-1)/ real (nphi*nfp)
      if( (i.eq.1).or.(i.eq.nphi2) )intgrate(i) = .5*delphi
   end do

   rbc = 0.
   rbs = 0.
   zbc = 0.
   zbs = 0.

   do mn = 1,mpnt
      mreal = m1(mn)
      nreal = n1(mn)/nfp
      m = mreal + 1
      dn= real (n1(mn))
      rbc(mreal,nreal) = 0.
      zbs(mreal,nreal) = 0.
      rbs(mreal,nreal) = 0.
      zbc(mreal,nreal) = 0.
      do i = 1,nphi2
         arg = dn*argi(i)
         tcosn = cos(arg)
         tsinn = sin(arg)
         rbc(mreal,nreal) = rbc(mreal,nreal)  &
                          + intgrate(i)*(tcosn*rmc(i,m) + tsinn*rms(i,m))
         zbs(mreal,nreal) = zbs(mreal,nreal)  &
                         + intgrate(i)*(tcosn*zms(i,m) - tsinn*zmc(i,m))
         zbc(mreal,nreal) = zbc(mreal,nreal)  &
                          + intgrate(i)*(tcosn*zmc(i,m) + tsinn*zms(i,m))
         rbs(mreal,nreal) = rbs(mreal,nreal)  &
                          + intgrate(i)*(tcosn*rms(i,m) - tsinn*rmc(i,m))
      end do

      if( (rbc(mreal,nreal).eq.0.).and.(zbs(mreal,nreal).eq.0.) .and.  &
          (zbc(mreal,nreal).eq.0.).and.(rbs(mreal,nreal).eq.0.) )  &
      cycle

      if ( ( mreal.eq.0).and.(nreal.eq.0) ) then
         rmnaxis(0) = 0.
         rmnaxis(0) = dot_product (intgrate,raxis)
         zmnaxis(0) = 0.
      else if (mreal.eq.0.and.nreal.gt.0) then
         rbc(0,nreal) = 2.0*rbc(0,nreal)
         zbs(0,nreal) = 2.0*zbs(0,nreal)
         rmnaxis(nreal) = 0.
         zmnaxis(nreal) = 0.
         do i = 1,nphi2
            rmnaxis(nreal) = rmnaxis(nreal)  &
                           + 2.*intgrate(i)*raxis(i)*cos(dn*argi(i))
            zmnaxis(nreal) = zmnaxis(nreal)  &
                           - 2.*intgrate(i)*zaxis(i)*sin(dn*argi(i))
         end do
      end if
   end do

end subroutine fftrans
