subroutine readin(mpol_ext, ntheta_ext, nphi_ext, niter_ext, nstep_ext,&
                  nfp_ext, ftol_ext, pexp_ext, qexp_ext, ok)
use mod_konst
use mod_curve3d
implicit none

integer, intent(in) :: mpol_ext, ntheta_ext, nphi_ext, niter_ext, nstep_ext, nfp_ext
double precision, intent(in):: ftol_ext, pexp_ext, qexp_ext


integer        :: i, j, ok

ok = 0

!   read(5, *) mpol, ntheta, nphi
mpol = mpol_ext
ntheta = ntheta_ext
nphi=nphi_ext
write(*,"(i6,i7,i6)") mpol, ntheta, nphi
if( (mod(nphi,2).ne.0) .and. nphi.ne.1 .or. nphi.eq.0 )  then     
   write(*,*) ' NPHI > 0 must be EVEN for non-symmetric systems'
   ok = -1
   return
end if

niter=niter_ext
nstep=nstep_ext
nfp=nfp_ext
write(*,"(3i7)") niter, nstep, nfp
if ( nfp.le.0 ) then
   write(*,*) ' NFP > 0 must be positive and exceed zero'
   ok = -1
   return
end if

ftol = ftol_ext
pexp = pexp_ext
qexp = qexp_ext
write(*,"(e10.2,2f8.2)") ftol, pexp, qexp

call init_konst
call init(ok)
if (ok /= 0) return

end subroutine readin
