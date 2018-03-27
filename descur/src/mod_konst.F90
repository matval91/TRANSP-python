MODULE mod_konst

implicit none
!jks fileunits changed to 0 (except standard out)

! file units
integer, parameter :: fdin = 0
integer, parameter :: fdout = 6
integer, parameter :: fdlog = 0  ! logfile
integer, parameter :: fdoc = 0   ! outcurve
integer, parameter :: fdoc1 = 0  ! outcurve
integer, parameter :: fdpo = 0  ! plotout
integer, parameter :: fdfit = 0 ! fitdat

! filenames
character(100) :: indat
character(100) :: logfile
character(100) :: outcurve1
character(100) :: outcurve2
character(100) :: plotout
character(100) :: fitdat

! input parameter
integer          :: mpol, ntheta, nphi, niter, nstep, nfp
double precision :: ftol, pexp, qexp

! parameter derived from input
integer :: mpol1, nphi2, mnd
integer :: mpol2, mpol3, mpol4, n2 

integer, parameter :: ncurve = 199      ! from printit

integer, parameter :: max_iterate = 5   ! getangle, evtl. einlesen ???

integer :: mpnt
integer :: nresets

double precision :: twopi, bmax

contains

subroutine init_konst

   ! set input parameters
   nresets = 0
   bmax = 0.15
   twopi = 8.*atan( real (1))

   mpol1 = mpol - 1
   nphi2 = 1 + nphi/2
   mnd = mpol * (nphi+1)
   mpol2 = 2*mpol
   mpol3 = mpol + mpol2
   mpol4 = mpol + mpol3
   n2 = 4*mpol + ntheta

end subroutine init_konst

END MODULE mod_konst
