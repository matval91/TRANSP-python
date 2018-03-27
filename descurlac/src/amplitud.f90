subroutine amplitud(rcenter,zcenter,angin,rmc,rms,zmc,zms,xpts,xin,yin)
   use mod_konst
   use mod_curve3d
   implicit none

   double precision, intent(in)  :: rcenter, zcenter
   double precision, intent(in)  :: angin(ntheta)
   double precision, intent(out) :: rmc(mpol), rms(mpol)
   double precision, intent(out) :: zmc(mpol), zms(mpol), xpts(ntheta)
   double precision, intent(in)  :: xin(ntheta), yin(ntheta)

   !***************************************************************
   !     This subroutine assigns initial guesses for angles and
   !     Fourier mode amplitudes to the appropriate components of
   !     the xvec array
   !**************************************************************

   integer :: j, m
   real    :: xmult, arg, xi, yi

   rmc =0; rms = 0; zmc = 0; zms = 0

   rmc(1) = rcenter
   zmc(1) = zcenter
   xpts = angin               ! call scopy (ntheta,angin,1,xpts,1)
   xmult = 2./ real (ntheta)
   do j = 1,ntheta
      arg = angin(j)
      xi = xmult*(xin(j) - rcenter)
      yi = xmult*(yin(j) - zcenter)
      do m = 2,mpol
         rmc(m) = rmc(m) + cos((m-1)*arg)*xi
         rms(m) = rms(m) + sin((m-1)*arg)*xi
         zmc(m) = zmc(m) + cos((m-1)*arg)*yi
         zms(m) = zms(m) + sin((m-1)*arg)*yi
      end do
   end do

   !
   !    INITIALIZE RMS = ZMC FOR M=1
   !
   ! rms(1) = zmc(1)    ????????????????????????
   rms(2) = 0.5*(rms(2) + zmc(2))
   zmc(2) = rms(2)

end subroutine amplitud
