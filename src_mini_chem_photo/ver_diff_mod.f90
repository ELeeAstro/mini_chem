module vert_diff_mod
  use, intrinsic :: iso_fortran_env
  implicit none


  integer, parameter :: dp = REAL64 ! Precision variable

  real(dp), parameter :: CFL = 0.40_dp
  real(dp), parameter :: kb = 1.380649e-16_dp
  real(dp), parameter :: amu = 1.66053906660e-24_dp

  public :: vert_diff
  private :: linear_interp

contains

  subroutine vert_diff(nlay, nlev, t_end, mu, grav_in, Tl, pl_in, pe_in, Kzz, nq, q, q0)
    implicit none

    integer, intent(in) :: nlay, nlev, nq
    real(dp), intent(in) :: t_end, grav_in
    real(dp), dimension(nlay), intent(in) :: Tl, pl_in, Kzz, mu
    real(dp), dimension(nlev), intent(in) :: pe_in
    real(dp), dimension(nq), intent(in) :: q0

    real(dp), dimension(nlay,nq), intent(inout) :: q

    integer :: k
    real(dp) :: h1, h2, d1Kzz, grav
    real(dp), dimension(nlev) :: alte, lpe, Kzz_e, Te, nde, pe
    real(dp), dimension(nlay) :: delz, delz_mid, pl

    real(dp), dimension(nlay) ::  lpl, lTl, nd
    real(dp), dimension(nlay,nq) :: flx
    real(dp), dimension(nq) :: phi_jp, phi_jm

    integer :: n_it, n
    real(dp) :: dt, t_now

    pl(:) = pl_in(:) * 10.0_dp   ! Convert pascal to dyne cm-2
    pe(:) = pe_in(:) * 10.0_dp   ! Convert pascal to dyne cm-2


    grav = grav_in * 100.0_dp

    !! First interpolate Kzz at the layer to Kzz at the edges
    !! Log pressure grid
    lpe(:) = log10(pe(:))
    lpl(:) = log10(pl(:))
    lTl(:) = log10(Tl(:))

    ! Find Kzz at levels
    Kzz_e(1) = Kzz(1)
    do k = 2, nlay
      Kzz_e(k) = (Kzz(k-1) + Kzz(k))/2.0_dp
    end do
    Kzz_e(nlev) = Kzz(nlay)

    do k = 2, nlay
      call linear_interp(lpe(k), lpl(k-1), lpl(k), lTl(k-1), lTl(k), Te(k))
    end do
    ! Edges are linearly interpolated
    Te(1) = (log10(Tl(1)) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * log10(Tl(1)/10.0_dp**Te(2)))
    Te(nlev) = (log10(Tl(nlay)) &
    &  + (log10(pe(nlev)/pe(nlay))/log10(pl(nlay)/pe(nlay))) * log10(Tl(nlay)/10.0_dp**Te(nlay)))
    ! De-log the temperature levels (edges)
    Te(:) = 10.0_dp**Te(:)

    ! Number density at layers and levels
    nd(:) = pl(:)/(kb * Tl(:))
    nde(:) = pe(:)/(kb * Te(:)) 

    !! First calculate the vertical height (cm) assuming hydrostatic equilibrium and differences
    alte(nlev) = 0.0_dp
    do k = nlev-1, 1, -1
      alte(k) = alte(k+1) + (kb*Tl(k))/(mu(k) * amu * grav) * log(pe(k+1)/pe(k))
      delz(k) = alte(k) - alte(k+1)
    end do

    !! Find differences between layers directly
    do k = 1, nlay-1
      delz_mid(k) = (alte(k) + alte(k+1))/2.0_dp - (alte(k+1) + alte(k+2))/2.0_dp
    end do
    delz_mid(nlay) = delz(nlay)

    !! We now follow the 1st order VULCAN scheme

    !! Prepare timestepping routine
    !! Find minimum timestep that satifies the CFL condition
    dt = t_end
    do k = 1, nlay
      dt = min(dt,CFL*(delz_mid(k))**2/Kzz_e(k+1))
    end do

    !! Begin timestepping routine
    t_now = 0.0_dp
    n_it = 1

    do while ((t_now < t_end) .and. (n_it < 100000))

      !! If next time step overshoots - last time step is equal tend
      if ((t_now + dt) > t_end) then
        dt = t_end - t_now
      end if

      !! Apply tracer lower boundary conditions
      q(nlay,:) = q0(:)
      q(1,:) = q(2,:)

      do n = 1, nq
        q(:,n) = max(q(:,n),1e-30_dp)
      end do

      !! Find flux between layers
      do k = 2, nlay-1
        phi_jp(:) = -Kzz_e(k) * nde(k) * (q(k-1,:) - q(k,:))/delz_mid(k-1)
        phi_jm(:) = -Kzz_e(k+1) * nde(k+1) * (q(k,:) - q(k+1,:))/delz_mid(k)
        flx(k,:) = -(phi_jp(:) - phi_jm(:))/delz(k)
        !print*, k,   phi_jp(:),  phi_jm(:), flx(k,:), dh(k)
      end do

      !! Perform tracer timestepping
      do k = 2, nlay-1
        q(k,:) = q(k,:) + flx(k,:)/nd(k)*dt
      end do

      !! Increment time and iteration number
      t_now = t_now + dt
      n_it = n_it + 1

    end do

    ! Apply boundary conditions
    q(1,:) = q(2,:)
    q(nlay,:) = q0(:)

    do n = 1, nq
      q(:,n) = max(q(:,n),1e-30_dp)
    end do


  end subroutine vert_diff

  ! Perform linear interpolation
  subroutine linear_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    norm = 1.0_dp / (x2 - x1)

    yval = (y1 * (x2 - xval) + y2 * (xval - x1)) * norm

  end subroutine linear_interp

end module vert_diff_mod

