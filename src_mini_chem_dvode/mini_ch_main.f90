program mini_chem_main
  use mini_ch_precision
  use mini_ch_class, only: g_sp
  use mini_ch_read_reac_list, only : read_react_list
  use mini_ch_i_dvode, only : mini_ch_dvode
  implicit none

  integer :: n, n_step, u_nml
  real(dp) :: T_in, P_in
  real(dp) :: t_step, t_now
  integer :: n_sp
  integer, parameter :: n_solver = 1
  real(dp), allocatable, dimension(:) :: VMR, VMR_IC, nd_out
  real(dp), allocatable, dimension(:,:) :: VMR_cp
  character(len=200) :: data_file, sp_file, network, net_dir, met

  integer, dimension(n_solver) :: u
  character(len=200), dimension(n_solver) :: integrator

  namelist /mini_chem/ T_in, P_in, t_step, n_step, n_sp, data_file, sp_file, network, net_dir, met
  namelist /mini_chem_VMR/ VMR

  ! Input Temperature [K], pressure [Pa] and stepping variables
  !! Read input variables from namelist
  open(newunit=u_nml, file='mini_chem.nml', status='old', action='read')
  read(u_nml, nml=mini_chem)
  !! Read VMRs from namelist
  allocate(VMR(n_sp),VMR_IC(n_sp), VMR_cp(n_solver, n_sp), nd_out(n_sp))
  read(u_nml, nml=mini_chem_VMR)
  close(u_nml)

  !! Commented out hard coded for testing
  ! T_in = 1500.0_dp
  ! P_in = 1.0e5_dp
  ! t_step = 30.0_dp
  ! n_step = 1
  !
  ! !! OH IC VMR of each species (NOTE: must be same order as mini_ch_sp.txt species)
  ! VMR(1) = 0.0_dp
  ! VMR(2) = 0.8_dp
  ! VMR(3) = 0.0_dp
  ! VMR(4) = 0.0_dp
  ! VMR(5) = 0.2_dp

  !! COH IC VMR of each species (NOTE: must be same order as mini_ch_sp.txt species)
  ! VMR(1) = 0.0_dp
  ! VMR(2) = 0.8_dp
  ! VMR(3) = 0.0_dp
  ! VMR(4) = 0.0_dp
  ! VMR(5) = 0.0_dp
  ! VMR(6) = 0.0_dp
  ! VMR(7) = 0.1_dp
  ! VMR(8) = 0.1_dp
  ! VMR(9) = 0.0_dp

  print*, 'T [K], P [bar], t_step, n_step, n_sp :'
  print*, T_in, P_in/1e5_dp, t_step, n_step, n_sp


  ! Save initial conditions to VMR_IC
  VMR_IC(:) = VMR(:)

  ! Initial time
  t_now = 0.0_dp

  ! Read the reaction and species list
  call read_react_list(data_file, sp_file, net_dir, met)

  print*, 'integrator: ', g_sp(:)%c, 'VMR sum'
  print*, 'IC: ', VMR_IC(:), sum(VMR_IC(:))

  VMR(:) = VMR_IC(:)
  do n = 1, n_solver
    VMR_cp(n,:) = VMR(:)
  end do


  integrator = (/'dvode     '/)
  do n = 1, n_solver
    open(newunit=u(n),file='outputs_dvode/'//trim(integrator(n))//'.txt',action='readwrite')
    write(u(n),*) 'n', 'time', g_sp(:)%c
    write(u(n),*) 0, 0.0, VMR_IC(:)
  end do

  !! Do time marching loop - in this test we call each solver to test them and their paramaters
  do n = 1, n_step

    ! Update time
    t_now = t_now + t_step

    ! Time now
    print*, n, n_step, t_now

    ! Call dvode - bdf method
    call mini_ch_dvode(T_in, P_in, t_step, VMR_cp(1,:), network)
    print*, 'dvode: ', VMR_cp(1,:), sum(VMR_cp(1,:))
    write(u(1),*) n, t_now, VMR_cp(1,:)

  end do

end program mini_chem_main
