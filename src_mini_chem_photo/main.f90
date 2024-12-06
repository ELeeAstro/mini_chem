program mini_chem_main
  use mini_ch_precision
  use mini_ch_class, only: g_sp
  use mini_ch_ce_interp, only : interp_ce_table
  use mini_ch_read_reac_list, only : read_react_list
  use mini_ch_i_dlsode, only : mini_ch_dlsode
  implicit none

  integer :: n, n_step, u_nml
  real(dp) :: t_step, t_now
  integer :: n_sp

  character(len=200) :: data_file, sp_file, network,  net_dir, met

  character(len=200) :: IC_file

  integer :: u
  character(len=200) :: integrator

  integer :: i, nlay, nlev
  real(dp) :: mu_z
  real(dp), allocatable, dimension(:) :: T_lay, P_lay, P_lev, mu, Kzz
  real(dp), allocatable, dimension(:,:) :: VMR, VMR_IC

  integer :: nwl
  character(len=200) :: stellar_file
  real(dp), allocatable, dimension(:,:) :: a_flux

  namelist /mini_chem_photo/ t_step, n_step, n_sp, data_file, sp_file, network, net_dir, met, &
    & nlay, IC_file, stellar_file

  !! Read input variables from namelist
  open(newunit=u_nml, file='mini_chem.nml', status='old', action='read')
  read(u_nml, nml=mini_chem_photo)
  close(u_nml)

  print*, 't_step, n_step, n_sp, nlay :'
  print*, t_step, n_step, n_sp, nlay

  nlev = nlay + 1

  ! Generate vertical T-p-Kzz profile

  ! Give intial conditions to VMR array
  allocate(VMR(nlay,n_sp),VMR_IC(nlay,n_sp))
  do i = 1, nlay
    call interp_ce_table(n_sp, T_lay(i), P_lay(i), VMR_IC(i,:), mu(i), IC_file)
  end do

  ! Read the reaction and species list
  call read_react_list(data_file, sp_file, net_dir, met)

  stop

  ! Save the inital conditions to file
  ! Rescale IC to 1
  VMR_IC(:) = VMR_IC(:)/sum(VMR_IC(:))
  print*, 'integrator: ', g_sp(:)%c, ' He', ' VMR sum'
  print*, 'IC: ', VMR_IC(:), sum(VMR_IC(:))

  integrator = 'dlsode'
  open(newunit=u,file='outputs_photo/'//trim(integrator)//'.txt',action='readwrite')
  write(u,*) 'n', 'time', g_sp(:)%c, ' He'
  write(u,*) 0, 0.0, VMR_IC(:)

  VMR(:) = max(VMR_IC(:),1e-30_dp)

  nwl = 1000
  allocate(a_flux(nwl))
  a_flux(:) = 0.0_dp
 

  ! Initial time
  t_now = 0.0_dp

  !! Do time marching loop
  ! - this loop emulates what a call to the model is like in the GCM
  do n = 1, n_step

    ! Update time
    t_now = t_now + t_step

    ! Time now
    print*, n, n_step, t_now

    !! Scale VMR to 1
    VMR(:) = max(VMR(:)/sum(VMR(:)),1e-30_dp)

    ! Call dlsode - bdf method - don't send He to integrator so dimensions are n_sp-1
    call mini_ch_dlsode(T_in, P_in, t_step, VMR(1:n_sp-1), nwl, a_flux(:), network)
    print*, 'dlsode: ', VMR(:), sum(VMR(:))
    write(u,*) n, t_now, VMR(:)

    !! Scale VMR to 1
    VMR(:) = max(VMR(:)/sum(VMR(:)),1e-30_dp)

  end do

end program mini_chem_main
