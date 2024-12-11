module mini_ch_init
  use mini_ch_precision
  use mini_ch_class
  implicit none

  public :: read_react_list, init_photochem
  !private ::


  ! Parameters for H
  real(dp), parameter :: wl_ly = 121.567_dp * 1.0e-7_dp ! Lyman alpha wavelength [cm]
  real(dp), parameter :: f_ly = c_s/wl_ly
  real(dp), parameter :: w_l = (2.0_dp * pi * f_ly) / 0.75_dp
  real(dp), dimension(10), parameter :: cp = (/1.26537_dp,3.73766_dp,8.8127_dp,19.1515_dp, &
  &  39.919_dp,81.1018_dp,161.896_dp,319.001_dp,622.229_dp,1203.82_dp/)

  real(dp), parameter :: sigT = 6.6524587051e-25_dp

contains

  subroutine read_react_list(data_file, sp_file, net_dir, met)
    implicit none

    character(len=200), intent(in) :: data_file, sp_file, net_dir, met
    integer :: i, u, j, k, u2

    !! Read in the reaction list
    open(newunit=u,file=trim(data_file),status='old',action='read',form='formatted')

    do i = 1, 7
      read(u,*)
    end do

    read(u,*) n_reac

    allocate(re(n_reac))

    do i = 1, n_reac
      re(i)%id = i
      read(u,*)
      read(u,*) re(i)%re_t, re(i)%Tmin, re(i)%Tmax
      if (re(i)%re_t == 2) then
        read(u,*) re(i)%n_re, re(i)%n_pr, re(i)%A, re(i)%B, re(i)%C
      else if (re(i)%re_t == 3) then
        read(u,*) re(i)%n_re, re(i)%n_pr, re(i)%A0, re(i)%B0, re(i)%C0, re(i)%Ainf, re(i)%Binf, re(i)%Cinf
      else if (re(i)%re_t == 4) then
        read(u,*) re(i)%n_re, re(i)%n_pr, re(i)%fname
      else if (re(i)%re_t == 5) then
        read(u,*) re(i)%n_re, re(i)%n_pr, re(i)%br_idx
      end if

      allocate(re(i)%stoi_re(re(i)%n_re), re(i)%stoi_pr(re(i)%n_pr))
      read(u,*) (re(i)%stoi_re(j),j=1,re(i)%n_re), (re(i)%stoi_pr(k),k=1,re(i)%n_pr)

      allocate(re(i)%c_re(re(i)%n_re),re(i)%c_pr(re(i)%n_pr))
      read(u,*) (re(i)%c_re(j),j=1,re(i)%n_re), (re(i)%c_pr(k),k=1,re(i)%n_pr)

      ! Read reaction rate table if from table
      if (re(i)%re_t == 4) then
        print*, 'Reading: ', trim(met)//'_'//trim(re(i)%fname)
        open(newunit=u2,file=trim(net_dir)//trim(met)//'_'//trim(re(i)%fname),status='old',action='read',form='formatted')
        read(u2,*)
        read(u2,*) re(i)%nT, re(i)%nP, re(i)%nkf
        allocate(re(i)%T(re(i)%nT), re(i)%P(re(i)%nP), re(i)%kf(re(i)%nT,re(i)%nP))
        allocate(re(i)%lT(re(i)%nT), re(i)%lP(re(i)%nP), re(i)%lkf(re(i)%nT,re(i)%nP))
        read(u2,*)
        read(u2,*) (re(i)%T(j), j = 1, re(i)%nT)
        re(i)%lT(:) = log10(re(i)%T(:))
        read(u2,*)
        read(u2,*) (re(i)%P(k), k = 1, re(i)%nP)
        ! Convert P from bar to dyne
        re(i)%P(:) = re(i)%P(:) * 1.0e6_dp
        re(i)%lP(:) = log10(re(i)%P(:))
        read(u2,*)
        do j = 1,  re(i)%nT
          do k = 1, re(i)%nP
            read(u2,*) re(i)%kf(j,k)
            re(i)%kf(j,k) = max(re(i)%kf(j,k),1.0e-199_dp)
          end do
        end do
        re(i)%lkf(:,:) = log10(re(i)%kf(:,:))
        close(u2)
      end if

      ! print*, (re(i)%T(j), j = 1, re(i)%nT)
      ! print*, (re(i)%P(k), k = 1, re(i)%nP)

      ! print*, re(i)%id, re(i)%re_t, re(i)%Tmin, re(i)%Tmax
      ! if (re(i)%re_t == 2) then
      !   print*, re(i)%n_re, re(i)%n_pr, re(i)%A, re(i)%B, re(i)%C
      ! else if (re(i)%re_t == 3) then
      !   print*, re(i)%n_re, re(i)%n_pr, re(i)%A0, re(i)%B0, re(i)%C0, re(i)%Ainf, re(i)%Binf, re(i)%Cinf
      ! else if (re(i)%re_t == 4) then
      !   print*, re(i)%n_re, re(i)%n_pr, re(i)%fname
      ! end if
      ! print*, (re(i)%stoi_re(j),j=1,re(i)%n_re), (re(i)%stoi_pr(k),k=1,re(i)%n_pr)
      ! print*, (re(i)%c_re(j),j=1,re(i)%n_re), ' -> ', (re(i)%c_pr(k),k=1,re(i)%n_pr)

      ! We now find the delta mu.the difference in the number of products - reactants
      re(i)%dmu = real(sum(re(i)%stoi_pr(:)) - sum(re(i)%stoi_re(:)),dp)

      !print*, i, re(i)%dmu

    end do

    close(u)

    !! Get the species list from the species input file
    open(newunit=u,file=trim(sp_file),status='old',action='read',form='formatted')
    do i = 1, 6
      read(u,*)
    end do
    read(u,*) n_sp
    read(u,*)

    allocate(g_sp(n_sp))

    do i = 1, n_sp
      g_sp(i)%id = i
      read(u,*) g_sp(i)%c, g_sp(i)%mw, g_sp(i)%thresh, g_sp(i)%n_a
      allocate(g_sp(i)%a_l(g_sp(i)%n_a))
      read(u,*) g_sp(i)%a_l(:)
      allocate(g_sp(i)%a_h(g_sp(i)%n_a))
      read(u,*) g_sp(i)%a_h(:)
      ! print*, g_sp(i)%id, g_sp(i)%c, g_sp(i)%mw,  g_sp(i)%n_a
      ! print*, g_sp(i)%a_l(:)
      ! print*, g_sp(i)%a_h(:)
    end do

    close(u)

    !! For each reaction, we find the id number for each of the reacting species
    do i = 1, n_reac
      ! Allocate the gas index arrays
      allocate(re(i)%gi_re(re(i)%n_re), re(i)%gi_pr(re(i)%n_pr))
      ! Find the reactant index
      do k = 1, re(i)%n_re
        do j = 1, n_sp
          if (g_sp(j)%c == re(i)%c_re(k)) then
            re(i)%gi_re(k) = j
            exit
          end if
        end do
      end do

      ! Find the product index
      do k = 1, re(i)%n_pr
        do j = 1, n_sp
          if (g_sp(j)%c == re(i)%c_pr(k)) then
            re(i)%gi_pr(k) = j
            exit
          end if
        end do
      end do

      !print*, i, re(i)%c_re(:), re(i)%gi_re(:)
      !print*, i, re(i)%c_pr(:), re(i)%gi_pr(:)
    end do

  end subroutine read_react_list

  subroutine init_photochem(nlay, dbin1, dbin2, dbin_12trans, wl_s, wl_e, stellar_file)
    implicit none

    integer, intent(in) :: nlay
    real(dp), intent(in) :: dbin1, dbin2, dbin_12trans, wl_s, wl_e
    character(len=200), intent(in) :: stellar_file

    integer :: i, l, n1, n2, u, nlines, io
    integer :: idx, idx1
    real(dp) :: wl_now
    real(dp), allocatable, dimension(:) :: wl_f, flx_f

    real(dp) :: A, B, C, nd_stp, King, n_ref, Dpol, a_vol
    real(dp) :: freq, w, wwl, xsec, wb

    integer :: j, k, p
    real(dp), allocatable, dimension(:) :: wl_xsec, axsec, dxsec, ixsec


    !! First calculate wavelength grid for photochemistry calculations
    !! We follow the vulcan method, splitting the range into 2 parts, one at width dbin1 and dbin2
    !! that transitions at dbin_12trans - then output number of wavelengths and generate wavelength grid

    !! Count number of intervals in first length
    n1 = 1
    wl_now = wl_s
    do while (wl_now < dbin_12trans)
      wl_now = wl_now + dbin1
      n1 = n1 + 1
      !print*, n1, wl_now,  dbin_12trans
    end do

    !! Remove one in case of overshoot
    n1 = n1 - 1

    !! Now count number of intervals in second length
    n2 = 1
    wl_now = dbin_12trans
    do while (wl_now < wl_e)
      wl_now = wl_now + dbin2
      n2 = n2 + 1
      !print*, n2, wl_now, wl_e
    end do

    !! Remove one in case of overshoot
    n2 = n2 - 1

    !! Generate grid of wavelength according to total number of intervals
    nwl = n1 + n2
    allocate(wl_grid(nwl), wn_grid(nwl)) ! Global variable

    wl_grid(1) = wl_s
    do i = 2, nwl-1
      if (i <= n1) then 
        wl_grid(i) = wl_grid(i-1) + dbin1
      else 
        wl_grid(i) = wl_grid(i-1) + dbin2
      end if
    end do
    wl_grid(nwl) = wl_e

    !! Convert nm to cm for integration purposes
    wl_grid(:) = wl_grid(:)*1.0e-7_dp
    wn_grid(:) = 1.0_dp/wl_grid(:)

    !! Allocate actinic flux array and stellar flux
    allocate(a_flux(nlay,nwl), s_flux(nwl))

    !! Read in stellar flux file and interpolate to wavelength grid
    !! Scale flux to distance and radius of planet

    open(newunit=u,file=trim(stellar_file),status='old',action='read',form='formatted')

    ! Read header
    read(u,*)

    ! Find number of lines in file
    nlines = 0
    do
      read(u,*,iostat=io)
      if (io /= 0) then 
        exit
      end if
      nlines = nlines + 1
    end do

    ! Allocate values for stellar flux file
    allocate(wl_f(nlines),flx_f(nlines)) 
    
    ! Rewind file
    rewind(u)
    ! Read header again
    read(u,*)

    ! Read file data
    do i = 1, nlines
      read(u,*) wl_f(i), flx_f(i)
      !print*, i, wl_f(i), flx_f(i)
    end do
    wl_f(:) = wl_f(:)*1.0e-7_dp

    close(u)


    ! Now we must interpolate the file to the decided wavelength grid
    do i = 1, nwl

      call locate(wl_f(:), nlines, wl_grid(i), idx)

      if (idx < 1) then
        s_flux(i) = 0.0_dp
      else if (idx >= nlines) then
        s_flux(i) = 0.0_dp
      else
        idx1 = idx + 1
        call linear_interp(wl_grid(i), wl_f(idx), wl_f(idx1), flx_f(idx), flx_f(idx1), s_flux(i))
      end if

      !print*, i, wl_grid(i)*1e7_dp, s_flux(i)
    end do

    ! Convert nm-1 to cm-1
    s_flux(:) = s_flux(:) * 1e7_dp

    deallocate(wl_f, flx_f)

    !! Now we must read in photochemical cross sections for each species and interpolate as
    !! well as calculate Rayleigh scattering for each species 

    ! Loop over each species and calculate Rayleigh xsec
    ! - we need smooth cross-sections
    ! So typically used functions in the field sometimes cannot be used
    do i = 1, n_sp

      allocate(g_sp(i)%Ray_xsec(nwl))

      !! Calculate Rayleigh scattering data if available
      select case(trim(g_sp(i)%c))

      case('OH')

        !! Use polarisability measurement
        a_vol = 6.965_dp / (1e8_dp)**3
        King = 1.0_dp
        do l = 1, nwl
          g_sp(i)%Ray_xsec(l) = 128.0_dp/3.0_dp * pi**5 * a_vol**2 * wn_grid(l)**4 * King
        end do

      case('H2')

        !! Use Irwin (2009)+ parameters
        A = 13.58e-5_dp ; B = 7.52e-3_dp
        nd_stp = 2.65163e19_dp
        King = 1.0_dp

        do l = 1, nwl
          n_ref = A * (1.0_dp + B/(wl_grid(l)*1e4_dp)**2) + 1.0_dp
          g_sp(i)%Ray_xsec(l) = ((24.0_dp * pi**3 * wn_grid(l)**4)/(nd_stp**2)) &
            & * ((n_ref**2 - 1.0_dp)/(n_ref**2 + 2.0_dp))**2  * King
        end do

      case('H2O')  

        g_sp(i)%Ray_xsec(:) = 0.0_dp

      case('H')

        do l = 1, nwl

          freq = c_s/wl_grid(l)
          w = 2.0_dp * pi * freq
          wwl = w/w_l

          ! Lee and Kim (2004)
          if (wwl <= 0.6_dp) then
            ! Low energy limit
            xsec = 0.0_dp
            do p = 0, 9
              xsec = xsec + (cp(p+1) * wwl**(2 * p))
            end do
            xsec = xsec * wwl**4
          else
            ! High energy limit (approaching Lyman alpha wavelengths)
            wb = (w - 0.75_dp*w_l)/(0.75_dp*w_l)
            xsec = (0.0433056_dp/wb**2)*(1.0_dp - 1.792_dp*wb - 23.637_dp*wb**2 - 83.1393_dp*wb**3 &
            & - 244.1453_dp*wb**4 - 699.473_dp*wb**5)
          end if
          ! Multiply by Thomson x-section
          g_sp(i)%Ray_xsec(l) = xsec * sigT
        end do

      case('CO')  

        !! Use Irwin (2009)+ parameters
        A = 32.7e-5_dp ; B = 8.1e-3_dp
        nd_stp = 2.65163e19_dp
        King = 1.0_dp

        do l = 1, nwl
          n_ref = A * (1.0_dp + B/(wl_grid(l)*1e4_dp)**2) + 1.0_dp
          g_sp(i)%Ray_xsec(l) = ((24.0_dp * pi**3 * wn_grid(l)**4)/(nd_stp**2)) &
            & * ((n_ref**2 - 1.0_dp)/(n_ref**2 + 2.0_dp))**2  * King
        end do
        
      case('CO2')

        !! Use Irwin (2009)+ parameters
        A = 43.9e-5_dp ; B = 6.4e-3_dp
        nd_stp = 2.65163e19_dp
        King = 1.0_dp

        do l = 1, nwl
          n_ref = A * (1.0_dp + B/(wl_grid(l)*1e4_dp)**2) + 1.0_dp
          g_sp(i)%Ray_xsec(l) = ((24.0_dp * pi**3 * wn_grid(l)**4)/(nd_stp**2)) &
            & * ((n_ref**2 - 1.0_dp)/(n_ref**2 + 2.0_dp))**2  * King
          !print*, l, wl_grid(l)*1e7_dp, n_ref, g_sp(i)%Ray_xsec(l)
        end do

      case('O') 

        !! Use polarisability measurement
        a_vol = 0.802_dp / (1e8_dp)**3
        King = 1.0_dp
        do l = 1, nwl
          g_sp(i)%Ray_xsec(l) = 128.0_dp/3.0_dp * pi**5 * a_vol**2 * wn_grid(l)**4 * King
        end do

      case('CH4')

        ! Use He et al. (2021) expression
        A = 3603.09_dp ; B = 4.40362e14_dp ; C = 1.1741e10_dp
        nd_stp = 2.546899e19_dp
        King = 1.0_dp

        do l = 1, nwl
          if (wl_grid(l)*1e7_dp < 92.5_dp) then
            g_sp(i)%Ray_xsec(l) = 0.0_dp
            cycle
          end if
          n_ref = max((A + (B / (C - wn_grid(l)**2)))/1e8_dp + 1.0_dp,1.0_dp+1e-7_dp)
          g_sp(i)%Ray_xsec(l) = ((24.0_dp * pi**3 * wn_grid(l)**4)/(nd_stp**2)) &
            & * ((n_ref**2 - 1.0_dp)/(n_ref**2 + 2.0_dp))**2  * King
         ! print*, l, wl_grid(l)*1e7_dp, n_ref, g_sp(i)%Ray_xsec(l)
        end do

      case('C2H2')

        !! Use polarisability measurement
        a_vol = 3.487_dp / (1e8_dp)**3
        King = 1.0_dp
        do l = 1, nwl
          g_sp(i)%Ray_xsec(l) = 128.0_dp/3.0_dp * pi**5 * a_vol**2 * wn_grid(l)**4 * King
        end do

      case('NH3')
        
        !! Use Irwin (2009)+ parameters
        A = 37.0e-5_dp ; B = 12.0e-3_dp; Dpol = 0.0922_dp
        nd_stp = 2.65163e19_dp
        King = (6.0_dp+3.0_dp*Dpol)/(6.0_dp-7.0_dp*Dpol)

        do l = 1, nwl
          n_ref = A * (1.0_dp + B/(wl_grid(l)*1e4_dp)**2) + 1.0_dp
          g_sp(i)%Ray_xsec(l) = ((24.0_dp * pi**3 * wn_grid(l)**4)/(nd_stp**2)) &
            & * ((n_ref**2 - 1.0_dp)/(n_ref**2 + 2.0_dp))**2  * King
        end do

      case('N2')

        !! Use Irwin (2009)+ parameters
        A = 29.06e-5_dp ; B = 7.7e-3_dp; Dpol = 0.030_dp
        nd_stp = 2.65163e19_dp

        King = (6.0_dp+3.0_dp*Dpol)/(6.0_dp-7.0_dp*Dpol)

        do l = 1, nwl
          n_ref = A * (1.0_dp + B/(wl_grid(l)*1e4_dp)**2) + 1.0_dp
          g_sp(i)%Ray_xsec(l) = ((24.0_dp * pi**3 * wn_grid(l)**4)/(nd_stp**2)) &
            & * ((n_ref**2 - 1.0_dp)/(n_ref**2 + 2.0_dp))**2  * King
        end do

      case('HCN')

        !! Use polarisability measurement
        a_vol = 2.593_dp / (1e8_dp)**3
        King = 1.0_dp
        do l = 1, nwl
          g_sp(i)%Ray_xsec(l) = 128.0_dp/3.0_dp * pi**5 * a_vol**2 * wn_grid(l)**4 * King
        end do 

      case('He')

        !! Use Irwin (2009)+ parameters
        A = 3.48e-5_dp ; B = 2.3e-3_dp
        nd_stp = 2.65163e19_dp
        King = 1.0_dp

        do l = 1, nwl
          n_ref = A * (1.0_dp + B/(wl_grid(l)*1e4_dp)**2) + 1.0_dp
          g_sp(i)%Ray_xsec(l) = ((24.0_dp * pi**3 * wn_grid(l)**4)/(nd_stp**2)) &
            & * ((n_ref**2 - 1.0_dp)/(n_ref**2 + 2.0_dp))**2  * King
        end do

      case default
        print*, 'Species not found: ', trim(g_sp(i)%c), ' - Setting zero Rayleigh xsec'
        g_sp(i)%Ray_xsec(:) = 0.0_dp
      end select

    end do


    !! Now read in and interpolate species photo cross-sections
    do i = 1, n_sp

      !! Read x-section data for species if available
      select case(trim(g_sp(i)%c))

      case('OH')

        g_sp(i)%nbr = 1; g_sp(i)%nbr_wl = 2

      case('H2')

        g_sp(i)%nbr = 1; g_sp(i)%nbr_wl = 2

      case('H2O')

        g_sp(i)%nbr = 3; g_sp(i)%nbr_wl = 4

      case('H')

        g_sp(i)%nbr = 1; g_sp(i)%nbr_wl = 2

      case('CO')

        g_sp(i)%nbr = 1; g_sp(i)%nbr_wl = 2

      case('CO2')

        g_sp(i)%nbr = 2; g_sp(i)%nbr_wl = 2

      case('O')

        g_sp(i)%nbr = 1; g_sp(i)%nbr_wl = 2

      case('CH4')

        g_sp(i)%nbr = 4; g_sp(i)%nbr_wl = 26

      case('C2H2')

        g_sp(i)%nbr = 1; g_sp(i)%nbr_wl = 2

      case('NH3')

        g_sp(i)%nbr = 2; g_sp(i)%nbr_wl = 201

      case('N2')

        g_sp(i)%nbr = 1; g_sp(i)%nbr_wl = 2

      case('HCN')

        g_sp(i)%nbr = 1; g_sp(i)%nbr_wl = 2

      case('He')

        g_sp(i)%nbr = 1; g_sp(i)%nbr_wl = 2

      case default        
        print*, 'Species not found: ', trim(g_sp(i)%c), ' - Setting zero photo xsec'
        g_sp(i)%nbr = 1; g_sp(i)%nbr_wl = 1
        allocate(g_sp(i)%ph_axsec(nwl),g_sp(i)%ph_dxsec(nwl,g_sp(i)%nbr),g_sp(i)%ph_ixsec(nwl))
        g_sp(i)%ph_axsec(:) = 0.0_dp
        g_sp(i)%ph_dxsec(:,:) = 0.0_dp
        g_sp(i)%ph_ixsec(:) = 0.0_dp
        cycle
      end select

      allocate(g_sp(i)%br_wl(g_sp(i)%nbr_wl), g_sp(i)%br(g_sp(i)%nbr_wl, g_sp(i)%nbr))
      allocate(g_sp(i)%ph_axsec(nwl), g_sp(i)%ph_dxsec(nwl,g_sp(i)%nbr), g_sp(i)%ph_ixsec(nwl))

      !! Read cross sections from file
      open(newunit=u,file='chem_data/ph_xsecs/'//trim(g_sp(i)%c)//'_cross.csv' & 
        & ,status='old',action='read',form='formatted')
      ! Read header
      read(u,*)

      ! Find number of lines in file
      nlines = 0
      do
        read(u,*,iostat=io)
        if (io /= 0) then 
          exit
        end if
        nlines = nlines + 1
      end do

      ! Allocate values for stellar flux file
      allocate(wl_xsec(nlines),axsec(nlines),dxsec(nlines),ixsec(nlines)) 
      
      ! Rewind file
      rewind(u)
      ! Read header again
      read(u,*)

      ! Read in cross sections
      do l = 1, nlines
        read(u,*) wl_xsec(l), axsec(l), dxsec(l), ixsec(l)
      end do
      wl_xsec(:) = wl_xsec(:)*1e-7_dp

      ! Interpolate to wl grid for simulations
      do l = 1, nwl

        call locate(wl_xsec(:), nlines, wl_grid(l), idx)

        if (idx < 1) then
          g_sp(i)%ph_axsec(l) = 0.0_dp
          g_sp(i)%ph_dxsec(l,:) = 0.0_dp
          g_sp(i)%ph_ixsec(l) = 0.0_dp
        else if (idx >= nlines) then
          g_sp(i)%ph_axsec(l) = 0.0_dp
          g_sp(i)%ph_dxsec(l,:) = 0.0_dp
          g_sp(i)%ph_ixsec(l) = 0.0_dp
        else
          idx1 = idx + 1
          call linear_interp(wl_grid(l), wl_xsec(idx), wl_xsec(idx1), axsec(idx), axsec(idx1), g_sp(i)%ph_axsec(l))
          call linear_interp(wl_grid(l), wl_xsec(idx), wl_xsec(idx1), dxsec(idx), dxsec(idx1), g_sp(i)%ph_dxsec(l,1))
          g_sp(i)%ph_dxsec(l,:) = g_sp(i)%ph_dxsec(l,1)
          call linear_interp(wl_grid(l), wl_xsec(idx), wl_xsec(idx1), ixsec(idx), ixsec(idx1), g_sp(i)%ph_ixsec(l))
        end if

        !print*, trim(g_sp(i)%c), i, wl_grid(l)*1e7_dp, g_sp(i)%ph_axsec(l), g_sp(i)%ph_dxsec(l), g_sp(i)%ph_ixsec(l)

      end do


      close(u)

      if (trim(g_sp(i)%c) == 'H' .or. trim(g_sp(i)%c) == 'O' .or. trim(g_sp(i)%c) == 'He') then
        g_sp(i)%br_wl(1) = wl_grid(1); g_sp(i)%br_wl(2) = wl_grid(nwl);
        g_sp(i)%br(:,:) = 1.0_dp
        deallocate(wl_xsec, axsec, dxsec, ixsec)
        cycle
      end if

      !! Read branch ratios
      open(newunit=u,file='chem_data/ph_xsecs/'//trim(g_sp(i)%c)//'_branch.csv' & 
        & ,status='old',action='read',form='formatted')

      read(u,*) ; read(u,*)

      do l = 1, g_sp(i)%nbr_wl
        read(u,*) g_sp(i)%br_wl(l),  g_sp(i)%br(l,:)
        !print*, trim(g_sp(i)%c), i, l,  g_sp(i)%br_wl(l),  g_sp(i)%br(l,:)
      end do
      g_sp(i)%br_wl(:) = g_sp(i)%br_wl(:)*1e-7_dp


      close(u)

      deallocate(wl_xsec, axsec, dxsec, ixsec)

      !! Now we weight the dissosiation cross-section with the branching ratio for each species
      do k = 1, g_sp(i)%nbr_wl-1
        do l = 1, nwl
          if ((wl_grid(l) >= g_sp(i)%br_wl(k)) .and. (wl_grid(l) <= g_sp(i)%br_wl(k+1))) then
            g_sp(i)%ph_dxsec(l,:) = g_sp(i)%ph_dxsec(l,:) * g_sp(i)%br(k,:)
          end if
        end do
      end do

    end do

    !! Finally, find index for wavelength grid for photodissociation threshold
    do i = 1, n_sp
      call locate(wl_grid(:), nwl, g_sp(i)%thresh*1e-7_dp, g_sp(i)%th_idx)
      if (g_sp(i)%th_idx < 1) then
        g_sp(i)%th_idx = 0
      else if (g_sp(i)%th_idx >= nwl) then
        g_sp(i)%th_idx = nwl
      end if
        
      !print*, i, g_sp(i)%th_idx

    end do

  end subroutine init_photochem

end module mini_ch_init
