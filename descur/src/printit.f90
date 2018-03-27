subroutine printit (ok,rbc,zbs,rbs,zbc,rmnaxis,zmnaxis)
   use mod_konst
   use mod_curve3d
   implicit none

   double precision, intent(in) :: rbc(0:mpol1,-nphi2:nphi2)
   double precision, intent(in) :: zbs(0:mpol1,-nphi2:nphi2)
   double precision, intent(in) :: rbs(0:mpol1,-nphi2:nphi2)
   double precision, intent(in) :: zbc(0:mpol1,-nphi2:nphi2)
   double precision, intent(in) :: rmnaxis(0:nphi2)
   double precision, intent(in) :: zmnaxis(0:nphi2)

   double precision :: r1(0:ncurve), z1(0:ncurve)

   integer          :: ok, ier, i, j, m, n, kphi, kt, k
   double precision :: tol, zeta1, pit, arg, theta
   character(100)   :: filename
   character(3)     :: zahl
   character(11)    :: klammer


   ok = 0
   open(fdpo, file=trim(plotout), action="write", form="formatted", iostat=ier)
   if (ier /= 0) then
      write(fdout,*) "error in opening ", trim(plotout), ": ", ier
      write(fdlog,*) "error in opening ", trim(plotout), ": ", ier
      ok = -1
      return
   end if

   open(fdoc1, file=trim(outcurve2), action="write", form="formatted", iostat=ier)
   if (ier /= 0) then
      write(fdout,*) "error in opening ", trim(outcurve2), ": ", ier
      write(fdlog,*) "error in opening ", trim(outcurve2), ": ", ier
      ok = -1
      return
   end if

!   write (fdpo,"(7i10)") mpol, ntheta, nphi, mpol-1, nphi2, nfp, mpnt
   do j = 1, nphi
      do i = 1, ntheta
         write (fdpo, "(1p4e14.6)") rin(i,j), zin(i,j)
      end do
   end do
   write(fdout,*) ' OUTPUTTING FOURIER COEFFICIENTS TO OUTCURVE FILE'
   write(fdlog,*) ' OUTPUTTING FOURIER COEFFICIENTS TO OUTCURVE FILE'

   !*******************************************************************
   !    This subroutine prints out data into the file "outcurve"
   !    The format of the modes is compatible with input to VMEC
   !*******************************************************************

   write (fdoc,10)
10 format(/'   MB  NB      RBC         RBS         ', &
               'ZBC         ZBS        RAXIS       ZAXIS')

  tol = 1.e-8 * rbc(0,0)

   do m = 0, mpol-1
      do n = -nphi2, nphi2
!         write (fdpo, "(1p4e14.6)") rbc(m,n),zbs(m,n),rbs(m,n),zbc(m,n)
         if (abs(rbc(m,n)) .lt. tol .and. abs(zbs(m,n)) .lt. tol)  cycle
         if (m .eq. 0 .and. n .ge. 0) then
            write(fdoc1,"(a13,1pe12.4,a)") "RAXIS       =", rmnaxis(n), ","
            write(fdoc1,"(a13,1pe12.4,a)") "ZAXIS       =", zmnaxis(n), ","
            write(fdoc,"(i5,i4,1p6e12.4)")  &
                 m,n,rbc(m,n),rbs(m,n),zbc(m,n),zbs(m,n),rmnaxis(n),zmnaxis(n)
         else
            write(fdoc,"(i5,i4,1p6e12.4)")   &
                 m,n,rbc(m,n),rbs(m,n),zbc(m,n),zbs(m,n)
         end if
         klammer = " "
         klammer(1:1) = "("
         k = 2
         write (zahl, "(i3)") n
         zahl = adjustl(zahl)
         klammer(k:k+len_trim(zahl)) = zahl(1:len_trim(zahl))
         k = k+len_trim(zahl)
         write (klammer(k:k),"(a)") ","
         k = k + 1
         write (zahl, "(i3)") m
         zahl = adjustl(zahl)
         klammer(k:k+len_trim(zahl)) = zahl(1:len_trim(zahl))
         k = k+len_trim(zahl)
         write (klammer(k:k),"(a)") ")"
         write (klammer(10:10),"(a)") "="
         write(fdoc1, &
            "(a3,a10,1pe12.4,a5,a10,1pe12.4,a5,a10,1pe12.4,a5,a10,1pe12.4,a1)") &
            "rbc", klammer, rbc(m,n), ", rbs", klammer, rbs(m,n), &
            ", zbc", klammer, zbc(m,n), ", zbs", klammer, zbs(m,n), ","
      end do
   end do

   !--------------------------------------------------------------
   ! ... 05.05.2000:cas from  subroutine totz

   write (fdpo, *) "&&"

   do kphi=1,nphi
      ! build fit output name
      write (zahl, "(i3)") kphi
      if (kphi < 10) then
         zahl(1:2) = "00"
      else if (kphi < 100) then
         zahl(1:1) = "0"
      end if
      do i = 1, 100
         filename(i:i) = fitdat(i:i)
         if (fitdat(i:i) == ".") then
            j = i
            filename(j:j+2) = zahl(1:3)
            j = j+3
            filename(j:j) = "."
            exit
         else if (fitdat(i:i) == " ") then
            filename(i:i+2) = zahl(1:3)
            j = i+2
            exit
         end if
      end do
      j = j + 1
      i = i + 1
      do while (fitdat(i:i) /= " ")
         filename(j:j) = fitdat(i:i)
         j = j + 1
         if (j > 100) then
            write(fdout,*) "filename buffer too short - name truncated"
            write(fdlog,*) "filename buffer too short - name truncated"
            exit
         end if
         i = i + 1
      end do

      open(fdfit, file=trim(filename), action="write", form="formatted", &
                  iostat=ier)
      if (ier /= 0) then
         write(fdout,*) "error in opening ", trim(filename), ": ", ier
         write(fdlog,*) "error in opening ", trim(filename), ": ", ier
         ok = -1
         return
      end if

      r1=0.
      z1=0.
      zeta1=twopi*(kphi-1)/ real (nfp*nphi)
      pit  =1./(ncurve-1.)
      do kt=1,ncurve
         theta=twopi*pit*(kt-1)
         do m=0,mpol-1
            do n=-nphi2,nphi2
               arg = m*theta - n*zeta1
               r1(kt-1)=r1(kt-1)+rbc(m,n)*cos(arg)+rbs(m,n)*sin(arg)
               z1(kt-1)=z1(kt-1)+zbs(m,n)*sin(arg)+zbc(m,n)*cos(arg)
            end do
         end do
      end do
      r1(ncurve)=r1(0)
      z1(ncurve)=z1(0)
      do kt=0,ncurve
         write(fdfit,*) r1(kt), z1(kt)
         write (fdpo, "(1p4e14.6)") r1(kt), z1(kt)
      end do
      write (fdpo, *) "&&"

      close(fdfit)
   end do

   close(fdpo)
   close(fdoc1)

end subroutine printit
