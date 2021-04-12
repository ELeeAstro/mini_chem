program mini_chem_main
  use mini_ch_precision
  use mini_ch_class, only: g_sp
  use mini_ch_read_reac_list, only : read_react_list
  use mini_ch_i_seulex, only : mini_ch_seulex
  use mini_ch_i_rodas, only : mini_ch_rodas
  use mini_ch_i_radau5, only : mini_ch_radau5
  use mini_ch_i_dvode, only : mini_ch_dvode
  use mini_ch_i_dlsode, only : mini_ch_dlsode
  implicit none

  integer :: n, n_step
  real(dp) :: T_in, P_in
  real(dp) :: t_step, t_now
  integer, parameter :: n_sp = 5
  ! integer, parameter :: n_sp = 11
  real(dp), dimension(n_sp) :: VMR, VMR_IC

  ! Input Temperature [K] and pressure [Pa]
  T_in = 1000.0_dp
  P_in = 1.0e5_dp

  ! Time step and number of steps
  t_step = 30.0_dp
  n_step = 1

  !! OH IC VMR of each species (NOTE: must be same order as mini_ch_sp.txt species)
  VMR(1) = 0.0_dp
  VMR(2) = 0.8_dp
  VMR(3) = 0.0_dp
  VMR(4) = 0.0_dp
  VMR(5) = 0.2_dp

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
  ! VMR(10) = 0.0_dp
  ! VMR(11) = 0.0_dp


  ! Save initial conditions to VMR_IC
  VMR_IC(:) = VMR(:)

  ! Initial time
  t_now = 0.0_dp

  ! Read the reaction and species list
  call read_react_list()

  ! Subroutine that can produce IC for the VMR (ggCHEM etc) - NOT USED
  !call CE_IC()

  print*, 'integrator: ', g_sp(:)%c, 'VMR sum'
  print*, 'IC: ', VMR_IC(:), sum(VMR_IC(:))


  !! Do time marching loop - in this test we call each solver to test them and their paramaters
  do n = 1, n_step

    ! Call seulex - implicit Euler solver
    call mini_ch_seulex(T_in, P_in, t_step, VMR)
    print*, 'seulex: ', VMR(:), sum(VMR(:))
    VMR(:) = VMR_IC(:)

    ! Call rodas O(4) Rosenbrock method
    call mini_ch_rodas(T_in, P_in, t_step, VMR)
    print*, 'rodas: ', VMR(:), sum(VMR(:))
    VMR(:) = VMR_IC(:)

    ! Call radau5 O(5) - implict Runge-Kutta method
    call mini_ch_radau5(T_in, P_in, t_step, VMR)
    print*, 'radau5: ', VMR(:), sum(VMR(:))
    VMR(:) = VMR_IC(:)

    ! Call dvode - bdf method
    call mini_ch_dvode(T_in, P_in, t_step, VMR)
    print*, 'dvode: ', VMR(:), sum(VMR(:))
    VMR(:) = VMR_IC(:)

    ! Call dlsode - bdf method
    call mini_ch_dlsode(T_in, P_in, t_step, VMR)
    print*, 'dlsode: ', VMR(:), sum(VMR(:))
    VMR(:) = VMR_IC(:)

    ! Update time
    t_now = t_now + t_step

    ! Time now
    print*, t_now

  end do

end program mini_chem_main
