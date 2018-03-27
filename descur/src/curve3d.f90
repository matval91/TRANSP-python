subroutine curve3d(mpol_ext, ntheta_ext, nphi_ext, niter_ext, nstep_ext,&
                   nfp_ext, ftol_ext, pexp_ext, qexp_ext, rin_ext, zin_ext,&
                   ok, rbc_ext, zbs_ext, rbs_ext, zbc_ext,&
                   rmnaxis_ext, zmnaxis_ext)

use mod_konst
use mod_curve3d
implicit none

!  THIS IS PROGRAM DESCUR - WHICH USES A STEEPEST DESCENT
!  ALGORITHM TO FIND A LEAST SQUARES APPROXIMATION TO AN 
!  ARBITRARY 3-D SPACE CURVE. ANGLE CONSTRAINTS, BASED ON A
!  MINIMUM SPECTRAL WIDTH CRITERION, ARE APPLIED.
!  THE CONSTRAINT IS SATISFIED BY TANGENTIAL VARIATIONS
!  ALONG THE CURVE.
!
!  THE MAIN PROGRAM SETS UP THE INITIAL POINTS TO BE FITS AND
!  THEN CALLS THE SUBROUTINE SCRUNCH, WHICH DOES THE ACTUAL
!  CURVE-FITTING.
!********************************************************************
!  REVISION 1:  January 26, 1989
!               Compute Fit to individual toroidal planes separately
!               and then Fourier decompose in phi, using optimized 
!               angle representation
!********************************************************************
!
!
!*******************************************************************
!  CREATE ARRAYS OF CURVE POINTS, RIN(THETA,PHI) , ZIN(THETA,PHI),
!  TO BE FIT TO AN OPTIMIZED FOURIER SERIES BY CALLING "SCRUNCH".
!  NOTE: N IS IN UNITS NORMED TO NFP (NUMBER OF FIELD PERIODS) IN 20 LOOP
!*******************************************************************
!
!  ALTERNATIVE INPUT (NTHETA = NO. DATA POINTS PER PHI PLANE)

integer, intent(in) :: mpol_ext, ntheta_ext, nphi_ext, niter_ext,&
                       nstep_ext, nfp_ext
double precision, intent(in):: ftol_ext, pexp_ext, qexp_ext
double precision, allocatable :: rbc(:,:), rbs(:,:), zbc(:,:), zbs(:,:)
double precision, allocatable :: rmnaxis(:), zmnaxis(:)
double precision, intent(in):: rin_ext(1:ntheta_ext,1:nphi_ext),&
                               zin_ext(1:ntheta_ext,1:nphi_ext)
integer, intent(out) :: ok
double precision, intent(out):: &
                         rbc_ext(0:mpol_ext-1,-1-nphi_ext/2:1+nphi_ext/2),&
                         rbs_ext(0:mpol_ext-1,-1-nphi_ext/2:1+nphi_ext/2),&
                         zbc_ext(0:mpol_ext-1,-1-nphi_ext/2:1+nphi_ext/2),&
                         zbs_ext(0:mpol_ext-1,-1-nphi_ext/2:1+nphi_ext/2)
double precision, intent(out):: rmnaxis_ext(0:nphi_ext),&
                                zmnaxis_ext(0:nphi_ext)


integer           :: i, j, ic(6)
call readin(mpol_ext, ntheta_ext, nphi_ext, niter_ext, nstep_ext, nfp_ext,&
            ftol_ext, pexp_ext, qexp_ext, ok)
if (ok /= 0) stop
do i = 1, nphi
  do j=1,ntheta
    rin(j,i)=rin_ext(j,i)
    zin(j,i)=zin_ext(j,i)
  end do
end do


allocate (rbc(0:mpol1,-nphi2:nphi2), stat = ic(1))
allocate (rbs(0:mpol1,-nphi2:nphi2), stat = ic(2))
allocate (zbc(0:mpol1,-nphi2:nphi2), stat = ic(3))
allocate (zbs(0:mpol1,-nphi2:nphi2), stat = ic(4))
allocate (rmnaxis(0:nphi), stat = ic(5))
allocate (zmnaxis(0:nphi), stat = ic(6))
if (ANY(ic /= 0)) then
   write(fdout,*) "allocation error: ", ic
   write(fdlog,*) "allocation error: ", ic
   stop
end if
rbc=0.
rbs=0.
zbc=0.
zbs=0.
rmnaxis=0.
zmnaxis=0.


!*******************************************************************
!    PERFORM OPTIMIZING CURVE-FIT
!*******************************************************************
call scrunch(ok, rbc, zbs, rbs, zbc, rmnaxis, zmnaxis)
if (ok /= 0) return


!*******************************************************************
!    PRINT OPTIMIZED TRANSFORM COEFFICIENTS TO "OUTCURVE" FILE
!    AND STORE RIN, ZIN, RBC, ZBS DATA IN PLOTOUT FILE
!
!    PLOT FIT TO DATA POINTS (XMGR FILE)
!*******************************************************************
do i=0,nphi
  rmnaxis_ext(i)=rmnaxis(i)
  zmnaxis_ext(i)=zmnaxis(i)
end do
do i=0,nphi2
  do j=0,mpol-1
    rbc_ext(j,i)=rbc(j,i)
    rbs_ext(j,i)=rbs(j,i)
    zbc_ext(j,i)=zbc(j,i)
    zbs_ext(j,i)=zbs(j,i)
    rbc_ext(j,-i)=rbc(j,-i)
    rbs_ext(j,-i)=rbs(j,-i)
    zbc_ext(j,-i)=zbc(j,-i)
    zbs_ext(j,-i)=zbs(j,-i)
  end do
end do

deallocate(rbc)
deallocate(rbs)
deallocate(zbc)
deallocate(zbs)
deallocate(rmnaxis)
deallocate(zmnaxis)

call done

end subroutine curve3d
