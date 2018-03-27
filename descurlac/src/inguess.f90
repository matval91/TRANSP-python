subroutine inguess(ok)
   use mod_konst
   use mod_curve3d
   implicit none

   !*******************************************************************
   !    This subroutine obtains initial guesses for the centroid at each
   !    toroidal plane. By default, the polar axis is set equal to this
   !    centroid.  It is imperative that the polar axis be enclosed by
   !    the surface (otherwise the computed theta angles will not span
   !    [0,2pi]). For certain complex cross-sections, this simple estimate
   !    of the polar axis may fail, and the user must supply values for
   !    raxis(i),zaxis(i).  In addition, for non-starlike domains, the
   !    points along the curve must be monitonically increasing as a
   !    function of the arclength from any fixed point on the curve. This
   !    ordering is attempted by the subroutine ORDER, but may fail in
   !    general if too few theta points are used.
   !*******************************************************************
   !    COMPUTE CENTROID
   !*******************************************************************

   integer :: ok, i, j

   ok = 0
   do i = 1,nphi
      do j = 1,ntheta
         r0n(i) = r0n(i)+rin(j,i)/ real (ntheta)
         z0n(i) = z0n(i)+zin(j,i)/ real (ntheta)
      end do
      raxis(i) = r0n(i)
      zaxis(i) = z0n(i)
   end do

   !*******************************************************************
   !    ORDER POINTS ON FLUX SURFACE AT EACH TOROIDAL ANGLES
   !*******************************************************************
   write(fdout,*) 'ORDERING SURFACE POINTS'
   write(fdlog,*) 'ORDERING SURFACE POINTS'

   do i = 1,nphi
      call order(ok,rin(1,i),zin(1,i),raxis(i),zaxis(i))
      if (ok /= 0) return
   end do

   !*******************************************************************
   !    COMPUTE OPTIMUM ANGLE OFFSETS FOR M=1 MODES
   !*******************************************************************
   call getangle

end subroutine inguess
