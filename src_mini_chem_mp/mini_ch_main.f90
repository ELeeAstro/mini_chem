program mini_chem_main
  use mini_ch_precision
  use mini_ch_class, only: g_sp
  use mini_ch_read_reac_list, only : read_react_list
  use mini_ch_i_seulex, only : mini_ch_seulex
  use mini_ch_i_rodas, only : mini_ch_rodas
  use mini_ch_i_radau5, only : mini_ch_radau5
  use omp_lib
  implicit none

  integer :: n, n_step, u_nml, k
  real(dp) :: T_in
  real(dp) :: P_in
  real(dp) :: t_step, t_now
  integer :: n_sp
  integer, parameter :: n_solver = 3
  real(dp), allocatable, dimension(:) :: VMR, VMR_IC, nd_out
  real(dp), allocatable, dimension(:,:) :: VMR_cp
  character(len=200) :: data_file, sp_file, network,  net_dir, met

  integer, dimension(n_solver) :: u
  character(len=200), dimension(n_solver) :: integrator

  namelist /mini_chem/ T_in, P_in, t_step, n_step, n_sp, data_file, sp_file, network,  net_dir, met
  namelist /mini_chem_VMR/ VMR

  ! Input Temperature [K], pressure [Pa] and stepping variables
  !! Read input variables from namelist
  open(newunit=u_nml, file='mini_chem.nml', status='old', action='read')
  read(u_nml, nml=mini_chem)
  !! Read VMRs from namelist
  allocate(VMR(n_sp),VMR_IC(n_sp), VMR_cp(n_solver, n_sp), nd_out(n_sp))
  read(u_nml, nml=mini_chem_VMR)
  close(u_nml)

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


  integrator = (/'seulex    ','rodas     ','radau5    '/)
  do n = 1, n_solver
    open(newunit=u(n),file='outputs_mp/'//trim(integrator(n))//'.txt',action='readwrite')
    write(u(n),*) 'n', 'time', g_sp(:)%c
    write(u(n),*) 0, 0.0, VMR_IC(:)
  end do

  !! Do time marching loop - in this test we call each solver to test them and their paramaters
  !$omp parallel default(shared)
  do n = 1, n_step

    !$omp master

    ! Update time
    t_now = t_now + t_step
    ! Time now
    print*, n, n_step, t_now

    !$omp task
    ! Call radau5 O(5) - implicit Runge-Kutta method
    call mini_ch_radau5(T_in, P_in, t_step, VMR_cp(3,:), network)
    write(u(3),*) n, t_now, VMR_cp(3,:)
    print*, omp_get_thread_num(), 'radau5: ', VMR_cp(3,:), sum(VMR_cp(3,:))
    !$omp end task

    !$omp task
    ! Call rodas O(4) Rosenbrock method
    call mini_ch_rodas(T_in, P_in, t_step, VMR_cp(2,:), network)
    write(u(2),*) n, t_now, VMR_cp(2,:)
    print*, omp_get_thread_num(), 'rodas: ', VMR_cp(2,:), sum(VMR_cp(2,:))
    !$omp end task

    !$omp task
    ! Call seulex method
    call mini_ch_seulex(T_in, P_in, t_step, VMR_cp(1,:), network)
    write(u(1),*) n, t_now, VMR_cp(1,:)
    print*, omp_get_thread_num(), 'seulex: ', VMR_cp(1,:), sum(VMR_cp(1,:))
    !$omp end task

    !$omp taskwait

    !$omp end master

  end do
  !$omp end parallel


end program mini_chem_main
