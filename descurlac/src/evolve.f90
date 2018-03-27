subroutine evolve(g11)
   use mod_konst
   use mod_curve3d
   implicit none

   double precision :: g11, ftest, dtau, otav, b1, fac
   integer :: i
 
   ftest = gnorm/g11
   dtau = abs(1.-ftest)
   g11 = gnorm
   otav = dtau/delt
   dtau = delt*otav+.001
   if(dtau.gt.bmax)dtau = bmax
   b1 = 1.-.5*dtau
   fac = 1./(1.+.5*dtau)
   do  i  =  1,n2
      xdot(i)  =  fac*(xdot(i)*b1 - delt*gvec(i))
      xvec(i)  =  xvec(i) + xdot(i)*delt
   end do

end subroutine evolve
