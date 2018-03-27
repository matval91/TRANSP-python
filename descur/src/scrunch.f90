subroutine scrunch(ok, rbc, zbs, rbs, zbc, rmnaxis, zmnaxis)

use mod_konst
use mod_curve3d
implicit none

integer :: ok
double precision, intent(out) :: rbc(0:mpol1,-nphi2:nphi2)
double precision, intent(out) :: rbs(0:mpol1,-nphi2:nphi2)
double precision, intent(out) :: zbc(0:mpol1,-nphi2:nphi2)
double precision, intent(out) :: zbs(0:mpol1,-nphi2:nphi2)
double precision, intent(out) :: rmnaxis(0:nphi)
double precision, intent(out) :: zmnaxis(0:nphi)


integer          :: first, ier, nplane, m, n, imodes, iter, modeno
integer          :: ntype, ntoff
double precision :: tottime, r10, gmin, g11, gout, fsq

tottime = 0
first = 1
deltf = 1

ok = 0

!*******************************************************************
!     INITIALIZE FIXED m,n ARRAYS AND FILL INPUT COMMON BLOCK
!********************************************************************
call fixaray 

!*******************************************************************
!    COMPUTE INITIAL GUESSES (MUST ORDER ANGLES MONOTONICALLY FOR
!    NON-STARLIKE DOMAINS)
!*******************************************************************
call inguess(ok)
if (ok /= 0) return

!*******************************************************************
!    BEGIN MAIN INTERATION LOOP
!*******************************************************************
!write(fdout,*) ' ITERATIONS    RMS ERROR    FORCE GRADIENT    <M>', &
!               '    m #''s <= ON'
!write(fdoc,*) ' ITERATIONS    RMS ERROR    FORCE GRADIENT    <M>', &
!               '    m #''s <= ON'
!write(fdlog,*) ' ITERATIONS    RMS ERROR    FORCE GRADIENT    <M>', &
!               '    m #''s <= ON'


do nplane = 1,nphi2
!   write(fdout,*) '                  Fitting toroidal plane # ', nplane
!   write(fdoc,*) '                  Fitting toroidal plane # ', nplane
!   write(fdlog,*) '                  Fitting toroidal plane # ', nplane

   !**************************************************************
   ! INITIALIZE M=0 and M=1 MODE AMPLITUDES
   !
   ! STACKING OF XVEC ARRAY (FOR FIXED TOROIDAL PLANE)
   ! XVEC(1,mpol):Rmcos            XVEC(1+mpol,2*mpol):Rmsin  
   ! XVEC(1+2*mpol,3*mpol):Zmcos   XVEC(1+3*mpol,4*mpol):Zmsin
   ! XVEC(4*mpol,n2): Theta angle
   !***************************************************************

   do n = 1, n2
      xstore(n)= 0.
      xdot(n)  = 0.
      xvec(n)  = 0.
   end do

   call amplitud(r0n(nplane),z0n(nplane),angle(1,nplane), &
                 xvec,xvec(1+mpol),xvec(1+mpol2),xvec(1+mpol3), &
                 xvec(1+mpol4),rin(1,nplane),zin(1,nplane) )


   r10 = .5*(abs(xvec(2      )) + abs(xvec(2+ mpol))  &
         +     abs(xvec(2+mpol2)) + abs(xvec(2+mpol3)))


   imodes = mpol
   delt = deltf

   do iter=1,niter
      call funct(xvec,xvec(1+mpol),xvec(1+mpol2),xvec(1+mpol3),  &
                 xvec(1+mpol4),gvec,gvec(1+mpol),gvec(1+mpol2),  &
                 gvec(1+mpol3),gvec(1+mpol4),fsq,r10,rin(1,nplane),  &
                 zin(1,nplane),imodes)

      if (iter == 1) then
         gmin = gnorm
         imodes = 3
         gout = sqrt(gnorm)
         modeno = mpol-1
!         write(fdout,"(i8,1p2e16.3,0pf10.2,i8)")   &
!                         iter,fsq,gout,specw,modeno
!         write(fdoc,"(i8,1p2e16.3,0pf10.2,i8)")  &
!                         iter,fsq,gout,specw,modeno
!         write(fdlog,"(i8,1p2e16.3,0pf10.2,i8)")  &
!                         iter,fsq,gout,specw,modeno

         if ( ((gnorm.lt.ftol**2).and.(imodes.eq.mpol)) &
              .or. (fsq.lt.ftol) )    exit
      elseif (iter == 2) then
          g11 = gnorm
          cycle
      else
         gmin = min(gmin,gnorm)
         call evolve(g11)
         !**********************************************************
         !  RUDIMENTARY TIME STEP CONTROL
         !**********************************************************
         if ( (gnorm/gmin).gt.1.e6 ) first = 2
         if ( (first.eq.2) .or. (gmin.eq.gnorm) ) then
            call restart(ok, first)
            if (ok /= 0) return
         end if
 
         if ( (gnorm.lt.1.e-4).and.(imodes.lt.mpol) )then
            imodes = imodes + 1
            call restart(ok, first)
            if (ok /= 0) return
            first = 2
            delt = delt/.95
         endif
         if ( (mod(iter,nstep).eq.0).or.(gnorm.lt.ftol**2) ) then
            gout = sqrt(gnorm)
            modeno = imodes - 1
!            write(fdout,"(i8,1p2e16.3,0pf10.2,i8)") &
!                         iter,fsq,gout,specw,modeno
!            write(fdoc,"(i8,1p2e16.3,0pf10.2,i8)") &
!                         iter,fsq,gout,specw,modeno
!            write(fdlog,"(i8,1p2e16.3,0pf10.2,i8)") &
!                         iter,fsq,gout,specw,modeno

            if ( ((gnorm.lt.ftol**2).and.(imodes.eq.mpol)) &
                 .or. (fsq.lt.ftol) )    exit
         end if
      end if
   end do

   !****************************************************************
   !  STORE FINAL ANSWER FOR FINAL PHI-TRANSFORM
   !****************************************************************
   do ntype = 1,4
      ntoff = mpol*(ntype-1)
      do m = 1,mpol
         result(nplane,m,ntype) = xvec(m+ntoff)
      end do
   end do

end do

!*******************************************************************
!    OUTPUT LOOP
!*******************************************************************
!write(fdout, 10) pexp, qexp
!write(fdoc, 10) pexp, qexp
!write(fdlog, 10) pexp, qexp
!10 format(/' ANGLE CONSTRAINTS WERE APPLIED ',/, &
!        ' BASED ON RM**2 + ZM**2 SPECTRUM WITH P = ', &
!        f8.2,' AND Q = ',f8.2,/)

!*******************************************************************
!    PERFORM PHI-FOURIER TRANSFORM
!*******************************************************************
call fftrans(result(1,1,1),result(1,1,2),result(1,1,3), &
             result(1,1,4),rbc,zbs,rbs,zbc,rmnaxis,zmnaxis)

end subroutine scrunch
