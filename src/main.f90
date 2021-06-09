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
      ICS(:)             !integrated cross section per l

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

  !! additional variables
  character(len=200) :: filename !! for file io
  integer :: status !! for error checking
  integer :: i_r, i_k !! for iterating through r, k grids

  !set kgrid parameters - leave this as is
  kg_Na = 30; kg_Nb = 30; kg_Np = 10
  kg_a = 0.85; kg_b = 2.5; kg_p = 4.0
  nkmax=kg_Na+kg_Nb+kg_Np+1

!>>> open data.in file and read input parameters
!>>>    note: energy should be read in electron volts
!>>>      and grid parameters in atomic units
  !! open input file for reading
  filename="data/input/data.in"
  open(unit=iounit, file=trim(adjustl(filename)), action="read", &
      iostat=status)

  !! handle invalid open
  if (status /= 0) then
    write (*, *) "[error] ", trim(adjustl(filename)), " could not be opened"
    call exit(status)
  end if

  !! read energy
  read(iounit, *) energy
  read(iounit, *) nrmax, dr
  read(iounit, *) zproj, lmin, lmax

  !! close input file
  close(iounit)

!>>> do any input validation you think is necessary here
  !! check input validity
  status = 0

  if (energy < 0d0) status = 1
  if (nrmax < 1) status = 1
  if (dr < 0d0) status = 1
  if (lmin < 0) status = 1
  if (lmax < 0) status = 1
  if (lmin > lmax) status = 1

  !! handle invalid input
  if (status /= 0) then
    write (*, *) "[error] input data is invalid"
    call exit(status)
  end if

!>>> convert the energy to atomic units and calculate the
!>>>      projectile momentum
  !! scale energy
  energy = energy / 27.2113862459d0

  !! set momentum
  k = sqrt(2d0*energy)

!>>> determine number of rgrid points nrmax
!>>>    note: nrmax should be even for simpson's integration
!>>>          to take into account that the r=0 point has been omitted
  !! set nrmax such that nrmax*dr >= rmax
  nrmax = ceiling(rmax / dr)

  !! increment nrmax if it is odd
  if (mod(nrmax, 2) == 1) nrmax = nrmax + 1

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
  !! v(r) = z*(1 + (1/r))*exp(-2r)
  do i_r = 1, nrmax
    V(i_r) = zproj*(1d0 + (1d0/rgrid(i_r)))*exp(-2d0*rgrid(i_r))
  end do

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
!>>>    to easily study the convergence you can write the ics to file as a function of l
!>>>    along with a running total so that the running total in the final line is your
!>>>    total ICS summed over l
  !! open dcs output file for writing
  filename="data/output/dcs.txt"
  open(unit=iounit, file=trim(adjustl(filename)), action="write", &
      iostat=status)

  !! handle invalid open
  if (status /= 0) then
    write (*, *) "[error] ", trim(adjustl(filename)), " could not be opened"
    call exit(status)
  end if

  !! write dcs to file
  write (iounit, *) "# theta, dcs(theta)"
  do ntheta = 1, nthetamax
    write (iounit, *) theta(ntheta), DCS(ntheta)
  end do

  !! close dcs output file
  close (iounit)

  !! open ics output file for writing
  filename="data/output/ics.txt"
  open(unit=iounit, file=trim(adjustl(filename)), action="write", &
      iostat=status)

  !! handle invalid open
  if (status /= 0) then
    write (*, *) "[error] ", trim(adjustl(filename)), " could not be opened"
    call exit(status)
  end if

  !! write ics to file
  write (iounit, *) "# l, ics(l), sum(ics(1:l))"
  do l = lmin, lmax
    write (iounit, *) l, ics(l), sum(ics(1:l))
  end do

  !! close ics output file
  close (iounit)

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
  !! calculate ics unscaled
  do l = lmin, lmax
    ics(l) = (2d0*l + 1d0)*((abs(Ton(l))) ** 2)
  end do

  !! scale by constant term
  ics(:) = ics(:) * ((4d0*(pi ** 3))/(k ** 4))

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
!>>>    by iterating over l and using the partial-wave
!>>>    expansion of f
  do ntheta = 1, nthetamax
    costheta = cos(theta * (pi/180d0))

    f(ntheta) = 0d0
    do l = lmin, lmax
      f(ntheta) = f(ntheta) + ((2d0*l + 1d0)*Ton(l)*PL(l, costheta))
    end do
    f(ntheta) = f(ntheta) * (-pi/(k ** 2))

  end do

!>>> obtain the DCS from the scattering amplitude
  do ntheta = 1, nthetamax
    dcs(ntheta) = (abs(f(ntheta))) ** 2
  end do

end subroutine compute_dcs

subroutine setup_rgrid(nrmax, dr, rgrid, rweights)
  implicit none
  integer, intent(in) :: nrmax
  real*8, intent(in) :: dr
  real*8, intent(out) :: rgrid(nrmax), rweights(nrmax)
  integer :: ir !index to iterate over r

!>>> iterate over r and populate the rgrid and rweights arrays
!>>>      - rweights should contain Simpson's integration weights:
!>>>        (4, 2, 4, 2, ..., 2, 4) * dr / 3.0
!>>>      - you can make use of the intrinsic MOD function for the
!>>>        alternating 4, 2 terms
!>>>      - note we have neglected the terms with a coefficient of 1 (rather
!>>>        than 4 or 2) since the first term (r=0) is skipped and the last term
!>>>        corresponds to the end of the radial grid where we assume all
!>>>        functions should be zero (and if not the grid is not large enough)
  !! calculate radial grid
  do ir = 1, nrmax
    rgrid(ir) = dr * ir
  end do

  !! calculate simpson's integration weights
  rweights(1:nrmax:2) = 4d0
  rweights(2:nrmax:2) = 2d0
  !! debug
  write (*, *) "[debug] <rweights>"
  do ir = 1, nrmax
    write (*, *) ii, rweights(ii)
  end do
  !! scale weights appropriately
  rweights(:) = rweights(:) * (dr/3d0)

end subroutine setup_rgrid

subroutine setup_contwaves(nkmax, kgrid, l, nrmax, rgrid, contwaves)
  implicit none
  integer, intent(in) :: nkmax, l, nrmax
  real*8, intent(in) :: kgrid(nkmax), rgrid(nrmax)
  real*8, intent(out) :: contwaves(nkmax,nrmax)
  real*8 :: ncontwaves(nkmax,nrmax)
  integer :: nk, nr !indices to loop over k and r
  real*8 :: E
  !! numerov variables
  real*8 :: s_grid(nrmax), g_grid(nrmax), v_grid(nrmax)
  integer :: fact_term
  integer :: status
  integer :: ii

!>>> iterate over k, populating the contwaves matrix
  !! initialise numerov grids
  s_grid(:) = 0d0
  v_grid(:) = ((l)*(l + 1d0))/(rgrid(:) ** 2)

  !! calculate double factorial term
  fact_term = 1
  do ii = 1, (2*l + 1), 2
    fact_term = fact_term * ii

    !! debug
    write (*, *) "[debug] <ii> = ", ii, "<fact_term> = ", fact_term
  end do

  !! iterate over k
  do nk = 1, nkmax
    !! todo: should inline this
    g_grid(:) = (kgrid(nk) ** 2) - v_grid(:)

    !! calculate boundary conditions
    contwaves(nk, 1:2) = ((rgrid(1:2)*kgrid(nk)) ** (l + 1))/dble(fact_term)

    !! perform numerov method
    call numerov_f(nrmax, rgrid(2)-rgrid(1), s_grid, g_grid, &
        contwaves(nk, :), status)

    !! handle numerov_f failing
    if (status /= 0) then
      write (*, *) "[error] numerov_f failed"
      call exit(status)
    end if
  end do

end subroutine setup_contwaves


!>>> your forwards Numerov subroutine can go here
!! numerov_f
!!
!! Brief:
!! Calculates the values of $y(x)$, defined by the differential equation
!! $$ y''(x) = -g(x)y(x) + s(x) $$ ,
!! on a grid to fourth-order accuracy using the forward Numerov method.
!!
!! Summary:
!! Suppose, $X = \{x_{1}, \dotsc, x_{n_{x}}\}$, is a grid with $n_{x}$ points,
!! such that $x_{i+1} - x_{i} = h$ for all $i = 1, \dotsc, n_{x}-1$,
!! and that $g(x), s(x) : \mathbb{R} \to \mathbb{R}$ have been evaluated on
!! this grid.
!! Suppose further that $y : \mathbb{R} \to \mathbb{R}$ is defined by the
!! differential equation,
!! $$ y''(x) = -g(x)y(x) + s(x) $$ ,
!! and the values of $y(x_{1}), y(x_{2})$ are known.
!! The values of $y(x_{i})$ for $i = 3, \dotsc, n_{x}$ are calculated to
!! fourth-order accuracy using the forward Numerov method.
!!
!! Input:
!! - `n_x` is the number of grid points.
!! - `step_size` is the distance between consecutive points on the grid $X$.
!! - `s_grid` is the evaluation of $s(x)$ on the grid $X$.
!! - `g_grid` is the evaluation of $g(x)$ on the grid $X$.
!!
!! Output:
!! - `y_grid` is the evaluation of $y(x)$ on the grid $X$, calculated via the
!!   forward Numerov method.
!! - `status` integer status code which takes the following values:
!!   - `status == 0` indicates successful execution;
!!   - `status > 0` indicates that the arguments were invalid;
!!   - `status == -1` indicates that a numerical error (NaN or infinity)
!!     occured during execution.
subroutine numerov_f (n_x, step_size, s_grid, g_grid, y_grid, status)
  integer , intent(in) :: n_x
  double precision , intent(in) :: step_size
  double precision , intent(in) :: s_grid(n_x)
  double precision , intent(in) :: g_grid(n_x)
  double precision , intent(out) :: y_grid(n_x)
  integer , intent(out) :: status
  double precision :: step_s_grid(n_x)
  double precision :: step_g_grid(n_x)
  integer :: ii

  ! check if arguments are valid
  status = 0

  if (n_x < 1) then
    ! no grid points
    status = 1
  end if

  if (step_size < 0d0) then
    ! non-positive step_size
    status = 2
  end if

  ! terminate subroutine if arguments are invalid, otherwise proceed
  if (status /= 0) then
    y_grid(:) = 0d0
    return
  end if

  ! perform forward-numerov method
  step_s_grid(:) = s_grid(:) * (step_size ** 2) / 12d0
  step_g_grid(:) = g_grid(:) * (step_size ** 2) / 12d0

  if (n_x >= 3) then
    do ii = 2, n_x - 1
      y_grid(ii+1) = &
          ((2*y_grid(ii)*(1d0 - 5d0*step_g_grid(ii))) &
          - (y_grid(ii-1)*(1d0 + step_g_grid(ii-1))) &
          + (step_s_grid(ii+1) + 10d0*step_s_grid(ii) + step_s_grid(ii-1)) &
          ) / (1d0 + step_g_grid(ii+1))

      ! check for and handle numerical error (NaN or Infinity) and terminates
      if ((y_grid(ii+1) /= y_grid(ii+1)) .or. (y_grid(ii+1) > huge(0d0))) then
        status = -1

        if (ii + 2 <= n_x) then
          y_grid(ii+2:n_x) = 0d0
        end if

        return
      end if
    end do
  end if

end subroutine numerov_f

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
