subroutine getangle
   use mod_konst
   use mod_curve3d
   implicit none

   !*******************************************************************
   !     Compute angle offset consistent with constraint Z1n = Z1,-n
   !    Note: This is done iteratively, since elongation is unknown
   !*******************************************************************

   integer          :: i, j, iterate
   double precision :: rcos(nphi), rsin(nphi) 
   double precision :: zcos(nphi), zsin(nphi), phiangle(nphi)
   double precision :: xc, yc, dnum, denom, delangle

   do i = 1,nphi
      do j = 1,ntheta
         angle(j,i) = twopi* real (j-1)/ real (ntheta)
      end do
   end do

   do iterate = 1,max_iterate
      do i = 1,nphi
         rcos(i) = 0.
         rsin(i) = 0.
         zcos(i) = 0.
         zsin(i) = 0.
         do j = 1,ntheta
            xc = rin(j,i) - r0n(i)
            yc = zin(j,i) - z0n(i)
            rcos(i) = rcos(i) + cos(angle(j,i))*xc
            rsin(i) = rsin(i) + sin(angle(j,i))*xc
            zcos(i) = zcos(i) + cos(angle(j,i))*yc
            zsin(i) = zsin(i) + sin(angle(j,i))*yc
         end do
      end do

      !****************************************************************
      ! Compute new angles starting from offset phiangle(i)
      !****************************************************************
       dnum = 0.
       denom= 0.
       do i = 1,nphi
          dnum = dnum + zsin(i)
          denom= denom+ rcos(i)
       end do
       elongate = dnum/denom
       delangle = 0.
       do i = 1,nphi
          phiangle(i) = atan2(elongate*zcos(i)-rsin(i),elongate*zsin(i)+rcos(i))
          delangle = max(delangle,abs(phiangle(i)))
          do j = 1,ntheta
             angle(j,i) = angle(j,i) + phiangle(i)
          end do
       end do
       if (delangle.lt.0.02) then
          write(fdout,*) ' Average elongation = ',elongate
!          write(fdoc,*) ' Average elongation = ',elongate
!          write(fdlog,*) ' Average elongation = ',elongate
          exit
       end if
   end do

end subroutine getangle
