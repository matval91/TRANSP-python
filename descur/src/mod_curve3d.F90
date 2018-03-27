MODULE mod_curve3d

implicit none

double precision, allocatable :: r0n(:), z0n(:), raxis(:), zaxis(:)
integer, allocatable :: mm(:), nn(:), m1(:), n1(:)

double precision, allocatable :: dm1(:), faccon(:), xmpq(:,:)
double precision :: gnorm, specw, delt, deltf
double precision, allocatable :: angle(:,:)
double precision              :: dnorm, elongate

double precision, allocatable :: xvec(:), gvec(:), xdot(:)
double precision, allocatable :: xstore(:), result(:,:,:)
double precision, allocatable :: rin(:,:), zin(:,:)

contains

   subroutine init(ok)
   use mod_konst

   integer :: ok, ic(19)

   ok = 0

   ic = 0
   allocate (r0n(nphi), stat = ic(1))
   allocate (z0n(nphi), stat = ic(2))
   allocate (raxis(nphi), stat = ic(3))
   allocate (zaxis(nphi), stat = ic(4))

   allocate (mm(mpol), stat = ic(5))
   allocate (nn(nphi), stat = ic(6))
   allocate (m1(mnd), stat = ic(7))
   allocate (n1(mnd), stat = ic(8))

   allocate (dm1(mpol), stat = ic(9))
   allocate (faccon(mpol), stat = ic(10))
   allocate (xmpq(mpol,4), stat = ic(11))

   allocate (angle(ntheta,nphi), stat = ic(12))

   allocate (xvec(n2), stat = ic(13))
   allocate (gvec(n2), stat = ic(14))
   allocate (xdot(n2), stat = ic(15))
   allocate (xstore(n2), stat = ic(16))
   allocate (result(nphi2,mpol,4), stat = ic(17))

   allocate (rin(ntheta,nphi), stat = ic(18))
   allocate (zin(ntheta,nphi), stat = ic(19))

   if (ANY(ic /= 0)) then
      write(fdout,*) "allocation error in init: ", ic
      write(fdlog,*) "allocation error in init: ", ic
      ok = -1
      return
   end if

   r0n = 0; z0n = 0
   raxis = 0; zaxis = 0
   mm = 0; nn = 0; m1 = 0; n1 = 0
   dm1 = 0; faccon = 0; xmpq = 0
   angle = 0
   xvec = 0; gvec = 0; xdot = 0; xstore = 0; result = 0
   rin = 0; zin = 0

end subroutine init


subroutine done
   if (allocated(r0n)) deallocate(r0n)
   if (allocated(z0n)) deallocate(z0n)
   if (allocated(raxis)) deallocate(raxis)
   if (allocated(zaxis)) deallocate(zaxis)

   if (allocated(mm)) deallocate(mm)
   if (allocated(nn)) deallocate(nn)
   if (allocated(m1)) deallocate(m1)
   if (allocated(n1)) deallocate(n1)

   if (allocated(dm1)) deallocate(dm1)
   if (allocated(faccon)) deallocate(faccon)
   if (allocated(xmpq)) deallocate(xmpq)

   if (allocated(angle)) deallocate(angle)

   if (allocated(xvec)) deallocate(xvec)
   if (allocated(gvec)) deallocate(gvec)
   if (allocated(xdot)) deallocate(xdot)
   if (allocated(xstore)) deallocate(xstore)
   if (allocated(result)) deallocate(result)

   if (allocated(rin)) deallocate(rin)
   if (allocated(zin)) deallocate(zin)

   end subroutine done

END MODULE mod_curve3d
