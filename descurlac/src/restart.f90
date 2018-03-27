subroutine restart(ok, first)
   use mod_konst
   use mod_curve3d
   implicit none

   !***********************************************************************
   !  This routine either stores an accepted value of the local solution
   !  (first=1) or reset the present solution to a previous value (first=2)
   !***********************************************************************

   integer :: ok, first, l

   ok = 0
   if (first == 1) then
      xstore = xvec
   else
      do l = 1, n2
         xdot(l) = 0.
         xvec(l) = xstore(l)
      end do
      delt = .95 * delt
      first = 1
      nresets = nresets + 1
      if (nresets .ge. 100) then
         write(fdout,*) ' Time step reduced 100 times without convergence'
         ok = -1
      end if
   end if

end subroutine restart
