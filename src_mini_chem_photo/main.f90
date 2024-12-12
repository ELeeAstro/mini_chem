program mini_chem_main
  use vert_diff_mod, only : vert_diff
  use mini_ch_precision
  use mini_ch_class, only : g_sp, locate, linear_interp
  use mini_ch_ce_interp, only : interp_ce_table
  use mini_ch_init, only : read_react_list, init_photochem
  use mini_ch_actinic_flux, only : calc_actinic_flux
  use mini_ch_i_dlsode_photo, only : mini_ch_dlsode_photo
  implicit none

  integer :: n, n_step, u_nml
  real(dp) :: t_step, t_now
  integer :: n_sp

  character(len=200) :: data_file, sp_file, network,  net_dir, met

  character(len=200) :: IC_file

  integer :: u
  character(len=200) :: integrator

  integer :: i, nlay, nlev
  real(dp) :: mu_z, Tirr, Tint, k_IR, k_V, grav, p_top, p_bot, gam, tau0
  real(dp), allocatable, dimension(:) :: Tl, pl, pe, mu, Kzz, tau_IRl, VMR0
  real(dp), allocatable, dimension(:,:) :: VMR, VMR_IC

  character(len=200) :: stellar_file
  real(dp) :: dbin1, dbin2, dbin_12trans, wl_s, wl_e
  real(dp) :: Rs, sm_ax

  integer :: nlines, io, idx, idx1
  real(dp), allocatable, dimension(:) :: T_f, p_f, Kzz_f

  namelist /mini_chem_photo/ t_step, n_step, n_sp, data_file, sp_file, network, net_dir, met, &
    & nlay, IC_file, stellar_file, Rs, sm_ax

  !! Read input variables from namelist
  open(newunit=u_nml, file='mini_chem.nml', status='old', action='read')
  read(u_nml, nml=mini_chem_photo)
  close(u_nml)

  print*, 't_step, n_step, n_sp, nlay :'
  print*, t_step, n_step, n_sp, nlay

  nlev = nlay + 1

  mu_z = 0.67_dp
  Tirr = 1400.0_dp
  Tint = 400.0_dp
  grav = 21.40_dp
  k_V = 6e-4_dp
  k_IR = 1e-3_dp
  p_top = 1e-7_dp * 1e5_dp
  p_bot = 1000.0_dp * 1e5_dp

  ! Generate vertical T-p-Kzz profile
  allocate(Tl(nlay),pl(nlay),pe(nlev),Kzz(nlay),tau_IRl(nlay))


  p_top = log10(p_top)
  p_bot = log10(p_bot)
  do i = 1, nlev
     pe(i) = 10.0_dp**((p_bot-p_top) * real(i-1,dp) / real(nlev-1,dp) + p_top) 
  end do
  do i = 1, nlay
  pl(i) = (pe(i+1) - pe(i)) / log(pe(i+1)/pe(i))
  end do
  p_top = 10.0_dp**p_top
  p_bot = 10.0_dp**p_bot
  ! gam = k_V/k_IR
  ! tau0 = k_IR/grav * p_bot
  ! tau_IRl(:) = (pl(:)/p_bot * tau0)
  ! mu_z = max(1e-6_dp,mu_z)
  ! Tl(:) = ((3.0_dp/4.0_dp) * Tint**4 * (tau_IRl(:) + 2.0_dp/3.0_dp))
  ! Tl(:) = Tl(:) + (mu_z * 3.0_dp * Tirr**4)/4.0_dp *  &
  !   & (2.0_dp/3.0_dp + mu_z/gam + ((gam/(3.0_dp*mu_z)) - mu_z/gam) * exp(-gam*tau_IRl(:)/mu_z))
  ! Tl(:) = Tl(:)**(1.0_dp/4.0_dp)
  ! Kzz(:) = 1e10_dp

  !! Read T-p_Kzz file and interpolate T and Kzz
  open(newunit=u,file='T_p_Kzz_profiles/atm_HD189_Kzz.txt',action='read')
  ! Read header
  read(u,*) ; read(u,*)
 ! Find number of lines in file
  nlines = 0
  do
    read(u,*,iostat=io)
     if (io /= 0) then 
        exit
    end if
    nlines = nlines + 1
  end do
  ! Allocate values for T-p profile
  allocate(T_f(nlines),p_f(nlines),Kzz_f(nlines))
  ! Rewind file
  rewind(u)
  ! Read header again
  read(u,*) ; read(u,*)
  do i = nlines, 1, -1
    read(u,*) p_f(i), T_f(i), Kzz_f(i)
  end do
  p_f(:) = p_f(:)/10.0_dp ! Convert to pa
  close(u)
  ! Interpolate to pressure grid
  do i = 1, nlay
    call locate(p_f(:), nlines, pl(i), idx)
    if (idx < 1) then
      Tl(i) = T_f(1)
      Kzz(i) = Kzz_f(1)
    else if (idx >= nlines) then
      Tl(i) = T_f(nlines)
      Kzz(i) = Kzz_f(nlines)
    else
      idx1 = idx + 1
      call linear_interp(pl(i), p_f(idx), p_f(idx1), T_f(idx), T_f(idx1), Tl(i))
      call linear_interp(pl(i), p_f(idx), p_f(idx1), Kzz_f(idx), Kzz_f(idx1), Kzz(i))
    end if
  end do


  !! Print T-p-Kzz profile
  print*, 'i, pl [bar], T[k], Kzz [cm2 s-1]'
  do i = 1, nlay
    print*, i, pl(i)/1e5_dp, Tl(i), Kzz(i)
  end do

  ! Output T-p-Kzz profile 
  open(newunit=u,file='outputs_photo/T_p_Kzz.txt',action='readwrite')
  do i = 1, nlay
    write(u,*) i, pl(i)/1e5_dp, Tl(i), Kzz(i)
  end do
  close(u)

  ! Give initial conditions to VMR array
  allocate(VMR(nlay,n_sp), VMR_IC(nlay,n_sp), mu(nlay), VMR0(n_sp))
  do i = 1, nlay
    call interp_ce_table(n_sp, Tl(i), pl(i), VMR_IC(i,:), mu(i), IC_file)
    !print*, i, pl(i)/1e5_dp, mu(i), VMR_IC(i,:)
  end do

  ! Read the reaction and species list
  call read_react_list(data_file, sp_file, net_dir, met)

  ! Save the initial conditions to file
  ! Rescale IC to 1
  do i = 1, nlay
    VMR_IC(i,:) = VMR_IC(i,:)/sum(VMR_IC(i,:))
  end do

  integrator = 'dlsode'
  open(newunit=u,file='outputs_photo/'//trim(integrator)//'.txt',action='readwrite')
  write(u,*) 'layer ', 'time ', g_sp(:)%c
  do i = 1, nlay
    write(u,*) i, 0.0, VMR_IC(i,:)
  end do

  do i = 1, nlay
    VMR(i,:) = VMR_IC(i,:)
  end do

  dbin1 = 0.1_dp  ! the uniform bin width < dbin_12trans (nm)
  dbin2 = 2.0_dp   ! the uniform bin width > dbin_12trans (nm)
  dbin_12trans = 240.0_dp 
  wl_s = 0.1_dp  ! Start wavelength
  wl_e = 300.0_dp ! End wavelength

  ! Initialise photochemistry - reads stellar flux, calculates wavelength grid and interpolates cross-sections
  call init_photochem(nlay, dbin1, dbin2, dbin_12trans, wl_s, wl_e, stellar_file)

  ! Initial time
  t_now = 0.0_dp

  !! Do time marching loop
  ! - this loop emulates what a call to the model is like in the GCM
  do n = 1, n_step

    !! Calculate the actinic flux in the column
    call calc_actinic_flux(nlay, n_sp, Tl(:), pl(:), pe(:), mu(:), VMR(:,:), grav, mu_z, Rs, sm_ax)

    ! Call dlsode - bdf method - don't send He to integrator so dimensions are n_sp-1
    do i = 1, nlay
      call mini_ch_dlsode_photo(i, Tl(i), pl(i), t_step, VMR(i,:), network)
    end do

    !! Perform vertical diffusion
    VMR0(:) = VMR(nlay,:)

    call vert_diff(nlay, nlay+1, t_step, mu(:), grav, & 
      & Tl(:), pl(:), pe(:), Kzz(:), n_sp, VMR(:,:), VMR0(:))

    do i = 1, nlay
      !! Scale VMR to 1
      VMR(i,:) = max(VMR(i,:)/sum(VMR(i,:)),1e-30_dp)
    end do

    ! Update time
    t_now = t_now + t_step

    ! Time now
    print*, n, n_step, t_now

    do i = 1, nlay
      write(u,*) i, t_now, VMR(i,:)
    end do

  end do

  close(u)

end program mini_chem_main
