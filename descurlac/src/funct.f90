subroutine funct(rmc,rms,zmc,zms,xpts,grc,grs,gzc,gzs,gpts, &
                 fsq,r10,xin,yin,mmax)
   use mod_konst
   use mod_curve3d
   implicit none

   double precision, intent(in)  :: rmc(mpol), rms(mpol)
   double precision, intent(in)  :: zmc(mpol), zms(mpol), xpts(ntheta)
   double precision, intent(out) :: grc(mpol), grs(mpol)
   double precision, intent(out) :: gzc(mpol), gzs(mpol), gpts(ntheta)
   double precision, intent(out) :: fsq
   double precision, intent(in)  :: r10
   double precision, intent(in)  :: xin(ntheta), yin(ntheta) 
   integer, intent(in)           :: mmax

   double precision :: yc(ntheta), ys(ntheta), gcon(ntheta), gtt(ntheta)
   double precision :: r1(ntheta), z1(ntheta), rt1(ntheta)
   double precision :: zt1(ntheta), arg(ntheta), rcon(ntheta)
   double precision :: zcon(ntheta), cosa(ntheta,mpol), sina(ntheta,mpol)
   double precision :: denom, t1, t2, tcon, fac
   integer          :: i, j, l, m


   ! FORCES: dW/dRmn = SUM(i)[ R(i) - RIN(i)]*cos(m*u[i]) ....
   !  dW/dZmn = SUM(i)[ Z(i) - ZIN(i)]*sin(m*u[i]) ....
   !  dW/du[i]=    rt(i)*[R(i)-RIN(i)] + zt(i)*[Z(i) - ZIN(i)]
   ! THE NORM ON THE ANGLE FORCE (dW/du[i]) FOLLOWS FROM NEWTONS
   ! LAW AND IS APPROXIMATELY GTT = RT**2 + ZT**2

 
   denom = 0.
   gnorm = 0.
   specw = 0.
   r1(:) = -xin(:)
   z1(:) = -yin(:)
   rcon = 0.
   zcon = 0.
   rt1  = 0.
   zt1  = 0.
   grc = 0.
   grs = 0.
   gzc = 0.
   gzs = 0.
   gpts = 0.

   !*******************************************************************
   !    COMPUTE SPECTRAL WIDTH OF CURVE
   !*******************************************************************

   do  m = 2,mmax
      t2 = (rmc(m)*rmc(m)+zmc(m)*zmc(m) + rms(m)*rms(m)+zms(m)*zms(m))*xmpq(m,1)
      denom = denom+t2
      specw = specw+xmpq(m,2)*t2
   end do
   specw = specw/denom

   !*******************************************************************
   !    COMPUTE CURVE AND CONSTRAINT FORCES
   !*******************************************************************

   do  m = 1,mmax
      do  l = 1,ntheta
         arg(l)    = dm1(m) * xpts(l)
         cosa(l,m) = cos(arg(l))
         sina(l,m) = sin(arg(l))
         gtt(l)  = rmc(m)*cosa(l,m) + rms(m)*sina(l,m)
         gcon(l) = zmc(m)*cosa(l,m) + zms(m)*sina(l,m)
         r1(l)   = r1(l)  + gtt(l)
         z1(l)   = z1(l)  + gcon(l)
         rt1(l)  = rt1(l) + dm1(m)*(rms(m)*cosa(l,m) - rmc(m)*sina(l,m))
         zt1(l)  = zt1(l) + dm1(m)*(zms(m)*cosa(l,m) - zmc(m)*sina(l,m))
         rcon(l) = rcon(l) + xmpq(m,4)*gtt(l)
         zcon(l) = zcon(l) + xmpq(m,4)*gcon(l)
      end do
   end do

   do l = 1,ntheta
      gtt(l)  = rt1(l)*rt1(l) + zt1(l)*zt1(l)
      gpts(l) = r1(l)*rt1(l)+z1(l)*zt1(l)
      gcon(l) = rcon(l)*rt1(l) + zcon(l)*zt1(l)
   end do

   t1 = -1.0
   do l = 1,ntheta
      t1 = max(t1,gtt(l))
   enddo

   !*******************************************************************
   !    COMPUTE MEAN-SQUARE DEVIATION BETWEEN POINTS AND FIT
   !*******************************************************************

   fsq = 0.
   do l = 1,ntheta
      gpts(l) = gpts(l)/t1
      fsq = fsq + .5*(r1(l)**2 + z1(l)**2)
   enddo
   fsq = sqrt(dnorm*fsq)/r10

   !*******************************************************************
   !    FILTER CONSTRAINT FORCE TO REMOVE ALIASING
   !*******************************************************************

   tcon = 1./(t1*sqrt(1.+xmpq(mmax,3)))
   do m = 2,mmax-1
      yc(m) = dot_product (cosa(:,m),gcon)*faccon(m)*tcon
      ys(m) = dot_product (sina(:,m),gcon)*faccon(m)*tcon
   enddo

   gcon = 0.
   do m = 2,mmax-1
      do l = 1,ntheta
         gcon(l) = gcon(l) + yc(m)*cosa(l,m) + ys(m)*sina(l,m)
      enddo
   enddo

   !*******************************************************************
   !    ADD CURVE AND CONSTRAINT FORCES
   !*******************************************************************

   do m = 1,mmax
      rcon(:) = r1(:)  + gcon(:)*rt1(:)*xmpq(m,3)
      zcon(:) = z1(:)  + gcon(:)*zt1(:)*xmpq(m,3)
   enddo

   do m = 1,mmax
      do l = 1,ntheta
         grc(m)  = grc(m) + cosa(l,m)*rcon(l)
         gzc(m)  = gzc(m) + cosa(l,m)*zcon(l)
         grs(m)  = grs(m) + sina(l,m)*rcon(l)
         gzs(m)  = gzs(m) + sina(l,m)*zcon(l)
      enddo
   enddo

   !*******************************************************************
   !    COMPUTE m=1 CONSTRAINT (ZC(m=1)=RS(m=1)) to the Forces
   !*******************************************************************

   gzc(2) = grs(2) + gzc(2)
   grs(2) = gzc(2)
   grc(1) = 0.5*grc(1)

   grc = dnorm*grc
   grs = dnorm*grs
   gzc = dnorm*gzc
   gzs = dnorm*gzs

   do j=1,mmax
      gnorm = gnorm + grc(j)*grc(j) + gzc(j)*gzc(j)    &
                    + grs(j)*grs(j) + gzs(j)*gzs(j)
   enddo
   gnorm = gnorm/r10**2

   do j = 1,ntheta
      gnorm = gnorm + dnorm*gpts(j)*gpts(j)
   enddo
   do m = 2,mmax
      fac = 1./(1. + tcon*xmpq(m,3))
      grc(m) = grc(m)*fac
      gzc(m) = gzc(m)*fac
      grs(m) = grs(m)*fac
      gzs(m) = gzs(m)*fac
   enddo

end subroutine funct
