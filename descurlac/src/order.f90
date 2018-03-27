subroutine order(ok,rval,zval,xaxis,yaxis)
   use mod_konst
   implicit none

   integer                      :: ok
   double precision             :: rval(ntheta), zval(ntheta)
   double precision, intent(in) :: xaxis, yaxis

   !*******************************************************************
   !    Program ORDER : orders points on a magnetic surface at a
   !    fixed toroidal plane and assigns right-handed circulation
   !    around flux surface. XAXIS, YAXIS:  Polar-type axis (must lie
   !    inside curve to check sign of rotation)
   !*******************************************************************

   integer          :: next, i, ip1, i1, j
   double precision :: tempr(ntheta), tempz(ntheta)
   double precision :: olddist, shortest, newdist, saver, savez, residue
   double precision :: x, y, dx, dy

   ok = 0
   olddist = 1.e20

   i1 = 1
   ip1 = 2
   shortest = 1.e20
   do i = 1,ntheta-1

      do j = ip1,ntheta
         if (i1.gt.1) &
            olddist = (rval(i1-1)-rval(j))**2 + (zval(i1-1)-zval(j))**2
         newdist = (rval(i1)-rval(j))**2 + (zval(i1)-zval(j))**2
         if ((newdist.le.olddist) .and. (newdist.lt.shortest)) then
            next = j
            shortest = newdist
         endif
      end do

      !****************************************************************
      !  Swap nearest point (next) with current point (ip1)
      !****************************************************************
      if (shortest.ge.1.e10) then
         saver = rval(i-1)
         rval(i-1) = rval(i)
         rval(i)= saver
         savez = zval(i-1)
         zval(i-1) = zval(i)
         zval(i)= savez
         i1 = i1 - 1
         ip1 = ip1 - 1
      else
         saver = rval(ip1)
         rval(ip1) = rval(next)
         rval(next)= saver
         savez = zval(ip1)
         zval(ip1) = zval(next)
         zval(next)= savez
         ip1 = i + 1
         i1 = i
         shortest = 1.e20
      end if

   end do

   !*******************************************************************
   !    Check that xaxis,yaxis is inside surface and
   !    ascertain that the angle rotates counterclockwise
   !    using Cauchy theorem in "complex"-plane
   !*******************************************************************
   residue = 0.
   do i = 1,ntheta-1
      x = 0.5*(rval(i)+rval(i+1)) - xaxis
      y = 0.5*(zval(i)+zval(i+1)) - yaxis
      dx= rval(i+1) - rval(i)
      dy= zval(i+1) - zval(i)
      residue = residue + (x*dy - y*dx)/(x**2 + y**2 + 1.e-10)
   end do
   x = 0.5*(rval(1)+rval(ntheta)) - xaxis
   y = 0.5*(zval(1)+zval(ntheta)) - yaxis
   dx= rval(1) - rval(ntheta)
   dy= zval(1) - zval(ntheta)
   residue = residue + (x*dy - y*dx)/(x**2 + y**2 + 1.e-10)

   if ( residue .lt.(-.90*twopi)) then
      do i = 2,ntheta
         j = ntheta - i + 2
         tempr(i) = rval(j)
         tempz(i) = zval(j)
      end do
      do i = 2,ntheta
         rval(i) = tempr(i)
         zval(i) = tempz(i)
      end do
   else if ( abs(residue) .lt. (.90*twopi) ) then
      write(fdout,*) ' The magnetic axis is not enclosed by boundary '
      write(fdoc,*) ' The magnetic axis is not enclosed by boundary '
      write(fdlog,*) ' The magnetic axis is not enclosed by boundary '
      ok = -1
      return
   end if

end subroutine order
