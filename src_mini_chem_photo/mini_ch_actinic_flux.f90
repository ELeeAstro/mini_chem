module mini_ch_actinic_flux
  use mini_ch_precision
  use mini_ch_class
  implicit none


  public :: calc_actinic_flux
  !private ::

contains


  subroutine calc_actinic_flux(nlay, n_sp, Tl, pl_in, pe_in, mu_in, VMR, grav_in, mu_z, Rs, sm_ax)
    implicit none

    integer, intent(in) :: nlay, n_sp
    real(dp), intent(in) :: grav_in, mu_z, Rs, sm_ax
    real(dp), dimension(nlay), intent(in) :: Tl, pl_in, mu_in
    real(dp), dimension(nlay+1), intent(in) :: pe_in
    real(dp), dimension(nlay,n_sp), intent(in) :: VMR

    integer :: i, l, n
    real(dp) :: grav, delta_tau, slant_delta_tau, layer_factor, photon_energy
    real(dp), dimension(nlay) :: nd_atm, rho, pl, k_tot, dpe
    real(dp), dimension(nlay+1) :: pe, tau_tot
    real(dp), dimension(nwl) :: Finc

    if (mu_z < 1e-2) then
      a_flux(:,:) = 0.0_dp
      return
    end if

    !! Find the number density of the atmosphere
    pl(:) = pl_in(:) * 10.0_dp   ! Convert pascal to dyne cm-2
    pe(:) = pe_in(:) * 10.0_dp    ! Convert pascal to dyne cm-2
    nd_atm(:) = pl(:)/(kb*Tl(:))  ! Find initial number density [cm-3] of atmosphere
    rho(:) = (pl(:) * mu_in(:) * amu)/(kb * Tl(:)) ! Mass density of atmosphere [g cm-3]
    grav = grav_in * 100.0_dp

    ! Calculate pressure difference in levels
    do i = 1, nlay
      dpe(i) = pe(i+1) - pe(i)
    end do

    !! Scale stellar flux to flux at planetary atmosphere
    Finc(:) = ((Rs * Rsun)/(sm_ax * au))**2 * s_flux(:)

    ! Do direct beam calculation for each wavelength
    do l = 1, nwl
      ! Find total opacity [cm2 g-1] in each layer
      do i = 1, nlay
        k_tot(i) = 0.0_dp
        do n = 1, n_sp
          k_tot(i) = k_tot(i) + &
            & (g_sp(n)%Ray_xsec(l) * nd_atm(i) * VMR(i,n))/rho(i) + (g_sp(n)%ph_axsec(l) * nd_atm(i) * VMR(i,n))/rho(i)
        end do
      end do

      ! Find total optical depth at each level - assume ghost later to zero pressure equal to uppermost opacity value
      tau_tot(1) = (k_tot(1) * pe(1)) / grav
      do i = 1, nlay
        tau_tot(i+1) = tau_tot(i) + (k_tot(i) * dpe(i)) / grav
      end do

      ! Exact layer-averaged direct-beam flux for piecewise-constant opacity.
      ! Use a series expansion for optically thin layers to avoid cancellation
      ! in (1 - exp(-slant_delta_tau))/slant_delta_tau.
      photon_energy = (h_p*c_s)/wl_grid(l)
      do i = 1, nlay
        delta_tau = tau_tot(i+1) - tau_tot(i)
        slant_delta_tau = delta_tau/mu_z
        if (abs(slant_delta_tau) < 1.0e-6_dp) then
          layer_factor = 1.0_dp - slant_delta_tau/2.0_dp + &
            & slant_delta_tau**2/6.0_dp - slant_delta_tau**3/24.0_dp
        else
          layer_factor = (1.0_dp - exp(-slant_delta_tau))/slant_delta_tau
        end if
        a_flux(i,l) = Finc(l) * exp(-tau_tot(i)/mu_z) * layer_factor / photon_energy
      end do

      !a_flux(:,l) = 0.0_dp
      !print*, l, tau_tot(:)
       
    end do

  end subroutine calc_actinic_flux
  
end module mini_ch_actinic_flux
