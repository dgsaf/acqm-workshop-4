!  This file is an outline of a code to perform a calculation of a
!    charged particle scattered by a spherically-symmetric short-range local
!    potential.
!  There are a number of comments throughout the file beginning with !>>> which
!    indicate sections of code you will need to complete.

! Last modified: May 25 2021
!                Liam Scarlett

program main

  use constants
  implicit none

  real*8, allocatable :: &
    Vmat(:,:),       & !V-matrix elements Vmat(kf,ki)
    V(:),            & !radial potential function V(r)
    rgrid(:),        & !radial grid
    rweights(:),     & !radial integration weights
    kgrid(:),        & !momentum-space grid
    kweights(:),     & !momentum-space integration weights (and Green's func)
    contwaves(:,:),  & !projectile radial continuum waves contwaves(k,r)
    DCS(:),          & !array to hold differential cross section - DCS(theta)
    theta(:),        & !array to hold values of theta - in degrees
    ICS(:),          & !integrated cross section per l

  real*8 :: &
    rmax,   & !max value of radial grid
    dr,     & !radial grid step size
    energy, & !projectile energy
    k,      & !projectile momentum
    kg_A, kg_B, kg_P !some parameters for setting up kgrid

  complex*16, allocatable :: Ton(:) !on-shell T-matrix element per l

  integer :: &
    nrmax,      & !number of rgrid points
    nkmax,      & !number of kgrid points
    zproj,      & !projectile charge 
    l,          & !partial-wave angular momentum
    lmin, lmax, & !min and max values of l
    iounit,     & !a unit number for input/ouput
    ntheta,     & !an index to iterate over theta
    nthetamax,  & !max number of theta
    kg_Na, kg_Nb, kg_Np !some more parameters for setting up kgrid
  
  !set kgrid parameters - leave this as is
    kg_Na = 30; kg_Nb = 30; kg_Np = 10
    kg_a = 0.85; kg_b = 2.5; kg_p = 4.0
    nkmax=kg_Na+kg_Nb+kg_Np+1

  !>>> open data.in file and read input parameters
  !    note: energy should be read in electron volts
  !      and grid parameters in atomic units

  !>>> do any input validation you think is necessary here

  !>>> convert the energy to atomic units and calculate the
  !      projectile momentum

  !>>> determine number of rgrid points nrmax
  !    note: nrmax should be even for simpson's integration
  !          to take into account that the r=0 point has been omitted 

  !allocate memory
    allocate(rgrid(nrmax),rweights(nrmax))
    allocate(kgrid(nkmax),kweights(nkmax))
    allocate(contwaves(nkmax,nrmax))
    allocate(V(nrmax))
    allocate(Ton(lmin:lmax),ICS(lmin:lmax))
    allocate(Vmat(nkmax,nkmax))

  !setup grids
    call setup_rgrid(nrmax, dr, rgrid, rweights)
    call setup_kgrid(k, nkmax, kg_Na, kg_a, kg_Nb, kg_b, kg_Np, kg_p, kgrid, kweights)
  
  !>>> define short-range potential V(r)

  !begin loop over angular momenta
  do l=lmin, lmax
    !populate contwaves matrix with a continuum wave for each off-shell k
      call setup_contwaves(nkmax,kgrid,l,nrmax,rgrid,contwaves)
    
    !evaluate the V-matrix elements  
      call calculate_Vmatrix(nkmax,kgrid,contwaves,nrmax,rgrid,rweights,V,Vmat)

    !solve the Lippman-Schwinger equation for the on-shell T-matrix
      call tmatrix_solver(nkmax,kgrid,kweights,Vmat,Ton(l))
  enddo

  !populate theta - theta = 0, 1,..., 180 degrees
  nthetamax = 181
  allocate(theta(nthetamax),DCS(nthetamax))
  do ntheta = 1, nthetamax
    theta(ntheta) = dble(ntheta-1)
  enddo

  !call subroutines to calculate DCS and ICS
    call compute_dcs(nthetamax, theta, lmin, lmax, Ton, k, DCS)
    call compute_ics(lmin, lmax, Ton, k, ICS)

  !>>> output the DCS and ICS to files
  !    to easily study the convergence you can write the ics to file as a function of l
  !    along with a running total so that the running total in the final line is your
  !    total ICS summed over l


end program main

subroutine compute_ics(lmin, lmax, Ton, k, ICS)
  use constants
  implicit none
  integer, intent(in) :: lmin, lmax
  complex*16, intent(in) :: Ton(lmin:lmax)
  real*8, intent(in) :: k
  real*8, intent(out) :: ICS(lmin:lmax)
  integer :: l

  !>>> populate the ICS array with the partial-wave ICS per l
end subroutine compute_ics

subroutine compute_dcs(nthetamax, theta, lmin, lmax, Ton, k, DCS)
  use constants
  implicit none
  integer, intent(in) :: nthetamax, lmin, lmax
  real*8, intent(in) :: theta(nthetamax), k
  complex*16, intent(in) :: Ton(0:lmax)
  real*8, intent(out) :: DCS(nthetamax)
  integer :: l, ntheta !loop indices
  real*8:: PL !Legendre polynomials - from file plql.f
  real*8 :: costheta !use this to store cos(theta in radians)
  complex*16 :: f(nthetamax) !scattering amplitude

  !>>> calculate the scattering amplitude f(theta) for each theta
  !    by iterating over l and using the partial-wave
  !    expansion of f

  !>>> obtain the DCS from the scattering amplitude

end subroutine compute_dcs

subroutine setup_rgrid(nrmax, dr, rgrid, rweights)
  implicit none
  integer, intent(in) :: nrmax
  real*8, intent(in) :: dr
  real*8, intent(out) :: rgrid(nrmax), rweights(nrmax)
  integer :: ir !index to iterate over r

  !>>> iterate over r and populate the rgrid and rweights arrays
  !      - rweights should contain Simpson's integration weights:
  !        (4, 2, 4, 2, ..., 2, 4) * dr / 3.0
  !      - you can make use of the intrinsic MOD function for the 
  !        alternating 4, 2 terms
  !      - note we have neglected the terms with a coefficient of 1 (rather than 4 or 2) since
  !        the first term (r=0) is skipped and the last term corresponds to the end of the
  !        radial grid where we assume all functions should be zero (and if not the grid is not large enough)

end subroutine setup_rgrid

subroutine setup_contwaves(nkmax, kgrid, l, nrmax, rgrid, contwaves)
  implicit none
  integer, intent(in) :: nkmax, l, nrmax
  real*8, intent(in) :: kgrid(nkmax), rgrid(nrmax)
  real*8, intent(out) :: contwaves(nkmax,nrmax)
  real*8 :: ncontwaves(nkmax,nrmax)
  integer :: nk, nr !indices to loop over k and r
  real*8 :: E
  
  !>>> iterate over k, populating the contwaves matrix                                 

end subroutine setup_contwaves

   
!>>> your forwards Numerov subroutine can go here




subroutine calculate_Vmatrix(nkmax,kgrid,contwaves,nrmax,rgrid,rweights,V,Vmat)
  use constants
  implicit none
  integer, intent(in) :: nkmax, nrmax
  real*8, intent(in) :: kgrid(nkmax), contwaves(nkmax,nrmax), rgrid(nrmax), rweights(nrmax), V(nrmax)
  real*8, intent(out) :: Vmat(nkmax,nkmax)
  integer :: nkf,nki !indices for looping over on- and off-shell k

  !>>> evaluate the V-matrix elements and store in the Vmat matrix
  !    note: the V-matrix is symmetric, make use of this fact to reduce the 
  !          amount of time spent in this subroutine

end subroutine calculate_Vmatrix
    
subroutine tmatrix_solver(nkmax,kgrid,kweights,Vmat,Ton)
  use constants
  implicit none
  integer, intent(in) :: nkmax
  real*8, intent(in) :: kgrid(nkmax), kweights(nkmax), Vmat(nkmax,nkmax)
  complex*16, intent(out) :: Ton !on-shell T-matrix element
  real*8 :: &
    Koff(nkmax-1), & !half-off-shell K-matrix elements
    Kon,           & !on-shell K-matrix element
    Von,           & !on-shell V-matrix element
    A(nkmax-1,nkmax-1) !Coefficient matrix for the linear system Ax=b
  integer :: j, ipiv(nkmax-1), info

  !>>> store the on-shell V-matrix element in Von

  !>>> populate the matrix A according to Eq (113) in the slides
  
  !>>> populate the vector Koff with the half-on-shell V-matrix elements (RHS of Eq (112)) 
 
  !Here is the call to DGESV
  call dgesv( nkmax-1, 1, A, nkmax-1, ipiv, Koff, nkmax-1, info )
  if(info /= 0) then
    print*, 'ERROR in dgesv: info = ', info
  endif

  !>>> Now use the half-on-shell K matrix which has been stored in Koff to get the on-shell K-matrix element Kon
 
  !>>> And then use Kon to get the on-shell T-matrix element Ton

end subroutine tmatrix_solver

!A subroutine provided for you to set up the kgrid and kweights
!note: the kgrid is setup with the on-shell point in the first element
!      and the corresponding kweights include the integration weights
!      AND the Green's function
subroutine setup_kgrid(k,nkmax,Na,a,Nb,b,Np,p,kgrid,kweights)
  implicit none
  real*8, intent(in) :: k
  integer, intent(in) :: Na, Nb, Np
  integer, intent(in) :: nkmax
  real*8, intent(out) :: kgrid(nkmax), kweights(nkmax)
  integer :: nk
  real*8, intent(in) ::  a, b, p
  real*8 :: grid1(nkmax-1), weight1(nkmax-1)

  call kgrid_igor(0.0,k,a,Na,b,Nb,p,Np,nkmax-1,grid1,weight1)
  
  kgrid(1) = k
  kgrid(2:nkmax) = grid1
  kweights(1) = 0.0d0
  kweights(2:nkmax) = weight1
  do nk=2, nkmax
    kweights(nk) = 2.0d0* kweights(nk) / (k**2 - kgrid(nk)**2)
  enddo
end subroutine setup_kgrid
