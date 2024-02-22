module mini_ch_i_dlsodes
  use mini_ch_precision
  use mini_ch_class
  use mini_ch_chem
  implicit none

  logical, parameter :: use_stiff = .True.
  real(dp) :: nd_atm

  public ::  mini_ch_dlsodes, RHS_update, jac_dummy, &
    & jac_HO, jac_CHO, jac_NCHO

contains

  subroutine mini_ch_dlsodes(T_in, P_in, t_end, VMR, network)
    implicit none

    real(dp), intent(in) :: T_in, P_in, t_end
    real(dp), dimension(n_sp), intent(inout) :: VMR
    character(len=200), intent(in) :: network

    integer :: ncall
    real(dp) :: P_cgs

    ! Time controls
    real(dp) :: t_begin, t_now, t_old, t_goal
    logical :: con = .False.

    ! DLSODES variables
    real(dp) :: rtol, atol
    real(dp), dimension(n_sp) :: y, y_old
    real(dp), allocatable, dimension(:) :: rwork
    integer, allocatable, dimension(:) :: iwork
    integer :: itol, itask, istate, iopt, mf
    integer :: rworkdim, iworkdim

    !! Find the number density of the atmosphere
    P_cgs = P_in * 10.0_dp   ! Convert pascal to dyne cm-2
    nd_atm = P_cgs/(kb*T_in)  ! Find initial number density [cm-3] of atmosphere

    allocate(Keq(n_reac), re_f(n_reac), re_r(n_reac))

    ! First find the reverse reaction coefficents (Keq)
    call reverse_reactions(T_in)
    ! Find the forward, backward and net reaction rates
    call reaction_rates(T_in, P_cgs, nd_atm)

    !! Find initial number density of all species from VMR
    y(:) = nd_atm * VMR(:)

    ! -----------------------------------------
    ! ***  parameters for the DLSODES solver  ***
    ! -----------------------------------------

    itask = 1
    istate = 1
    iopt = 1

    ! Method flag
    if (use_stiff .eqv. .True.) then
      ! Problem is stiff (usual)
      ! mf = 121 - full jacobian matrix with jacobian save
      ! mf = 222 - internal calculated jacobian
      mf = 222
      rworkdim = 20 + int((2.0_dp + 1.0_dp/2.0_dp)*n_sp**2 + (11.0_dp + 9.0_dp/2.0_dp)*n_sp)
      iworkdim = 30
      allocate(rwork(rworkdim), iwork(iworkdim))

      itol = 1
      rtol = 1.0e-3_dp           ! Relative tolerances for each scalar
      atol = 1.0e-99_dp               ! Absolute tolerance for each scalar (floor value)

      rwork(:) = 0.0_dp
      iwork(:) = 0

      iwork(5) = 2
      iwork(7) = 1

    else
      ! Problem is not too stiff (not typical)
      ! mf = 10 - full jacobian matrix with jacobian save
      mf = 10
      rworkdim = 20 + 16*n_sp
      iworkdim = 30
      allocate(rwork(rworkdim), iwork(iworkdim))
      itol = 4
      rtol = 1.0e-3_dp
      atol = 1.0e-99_dp

    end if

    t_begin = 0.0_dp
    t_now = t_begin

    ! Set the printing flag
    ! 0 = no printing, 1 = printing
    call xsetf(1)

    ncall = 0

    do while (t_now < t_end)

      y_old(:) = y(:)
      t_old = t_now

      select case(network)
      case('HO')
        call DLSODES (RHS_update, n_sp, y, t_now, t_end, itol, rtol, atol, itask, &
        & istate, iopt, rwork, rworkdim, iwork, iworkdim, jac_HO, mf)
      case('CHO')
        call DLSODES (RHS_update, n_sp, y, t_now, t_end, itol, rtol, atol, itask, &
        & istate, iopt, rwork, rworkdim, iwork, iworkdim, jac_CHO, mf)
      case('NCHO')
        call DLSODES (RHS_update, n_sp, y, t_now, t_end, itol, rtol, atol, itask, &
        & istate, iopt, rwork, rworkdim, iwork, iworkdim, jac_dummy, mf)
      case default
        print*, 'Invalid network provided: ', trim(network)
        stop
      end select

      ! call check_con(n_sp,y(:),y_old(:),t_now,t_old,con)
      ! if (con .eqv. .True.) then
      !   exit
      ! end if

      ncall = ncall + 1

      if (mod(ncall,50) == 0) then
        istate = 1
      else  if (istate == -1) then
        istate = 2
      else if (istate < -1) then
        print*, istate
        exit
      end if

    end do

    VMR(:) = y(:)/nd_atm

    deallocate(Keq, re_r, re_f, rwork, iwork)

  end subroutine mini_ch_dlsodes

  subroutine RHS_update(NEQ, time, y, f, rpar, ipar)
    implicit none

    integer, intent(in) ::  NEQ
    real(dp), intent(inout) :: time
    real(dp), dimension(NEQ), intent(inout) :: y
    real(dp), dimension(NEQ), intent(inout) :: f
    real(dp), intent(inout) :: rpar
    integer, intent(inout) :: ipar

    integer :: i, k, j2, j3, j4, j5
    real(dp) :: msum, msum2, frate, rrate
    real(dp), dimension(n_reac) :: net_pr, net_re
    real(dp), dimension(NEQ) :: f1_pr, f2_pr, f3_pr, f4_pr, f5_pr
    real(dp), dimension(NEQ) :: f1_re, f2_re, f3_re, f4_re, f5_re

    ! Calculate the rate of change of number density for all species [cm-3/s]
    ! this is the f vector

    ! Loop through reactions add rates to the f array
    do i = 1, n_reac
      ! Do the forward and backward flux calculation for each speices in the reaction

      ! Find number density multiple for reactants in reaction
      msum = y(re(i)%gi_re(1))
      do k = 2, re(i)%n_re
         msum = msum * y(re(i)%gi_re(k))
      end do

      ! Find number density multiple for products in reaction
      msum2 = y(re(i)%gi_pr(1))
      do k = 2, re(i)%n_pr
         msum2 = msum2 * y(re(i)%gi_pr(k))
      end do

      if (re(i)%re_t == 3) then
        ! Mutliply both msum and msum2 by atmosphere nd if neutral body involved
        msum = msum * nd_atm
        msum2 = msum2 * nd_atm
      end if

      frate = msum * re_f(i)
      rrate = msum2 * re_r(i)

      ! Find flux for products - regular addition
      !f(re(i)%gi_pr(:)) = f(re(i)%gi_pr(:)) + frate - rrate
      !f(re(i)%gi_re(:)) = f(re(i)%gi_re(:)) + rrate - frate

      net_pr(i) = frate - rrate
      net_re(i) = rrate - frate

    end do

    !! Perform peicewise summation over the arrays
    !! here we just assume split into 5 blocks since n_reac is quite small
    f1_pr(:) = 0.0_dp
    f2_pr(:) = 0.0_dp
    f3_pr(:) = 0.0_dp
    f4_pr(:) = 0.0_dp
    f5_pr(:) = 0.0_dp

    f1_re(:) = 0.0_dp
    f2_re(:) = 0.0_dp
    f3_re(:) = 0.0_dp
    f4_re(:) = 0.0_dp
    f5_re(:) = 0.0_dp

    do i = 1, n_reac/5

      j2 = i + n_reac/5
      j3 = i + 2*n_reac/5
      j4 = i + 3*n_reac/5
      j5 = i + 4*n_reac/5

      f1_pr(re(i)%gi_pr(:)) = f1_pr(re(i)%gi_pr(:)) + net_pr(i)
      f2_pr(re(j2)%gi_pr(:)) = f2_pr(re(j2)%gi_pr(:)) + net_pr(j2)
      f3_pr(re(j3)%gi_pr(:)) = f3_pr(re(j3)%gi_pr(:)) + net_pr(j3)
      f4_pr(re(j4)%gi_pr(:)) = f4_pr(re(j4)%gi_pr(:)) + net_pr(j4)
      f5_pr(re(j5)%gi_pr(:)) = f5_pr(re(j5)%gi_pr(:)) + net_pr(j5)     
     
      f1_re(re(i)%gi_re(:)) = f1_re(re(i)%gi_re(:)) + net_re(i)
      f2_re(re(j2)%gi_re(:)) = f2_re(re(j2)%gi_re(:)) + net_re(j2)
      f3_re(re(j3)%gi_re(:)) = f3_re(re(j3)%gi_re(:)) + net_re(j3)
      f4_re(re(j4)%gi_re(:)) = f4_re(re(j4)%gi_re(:)) + net_re(j4)
      f5_re(re(j5)%gi_re(:)) = f5_re(re(j5)%gi_re(:)) + net_re(j5)
 
    end do

    f(:) = (f1_pr(:) + f2_pr(:) + f3_pr(:) + f4_pr(:) + f5_pr(:)) + &
      & (f1_re(:) + f2_re(:) + f3_re(:) + f4_re(:) + f5_re(:))

  end subroutine RHS_update

  subroutine jac_dummy(NEQ, T, Y, J, IAN, JAN, PDJ)
    integer, intent(in) :: NEQ, J
    real(dp), intent(in) :: T
    real(dp), dimension(NEQ), intent(in) :: Y, IAN, JAN
    real(dp), dimension(NEQ), intent(inout) :: PDJ
  end subroutine jac_dummy

  subroutine jac_NCHO(NEQ, T, Y, J, IAN, JAN, PDJ)
    implicit none
    integer, intent(in) :: NEQ, J
    real(dp), intent(in) :: T
    real(dp), dimension(NEQ), intent(in) :: Y, IAN, JAN
    real(dp), dimension(NEQ), intent(inout) :: PDJ

    pdj(1,1) = -re_f(1)*y(2) - re_f(2)*y(5) - re_r(3)*y(4)
    pdj(1,2) = -re_f(1)*y(1) + re_f(3)*y(7)
    pdj(1,3) = re_r(1)*y(4)
    pdj(1,4) = re_r(1)*y(3) + re_r(2)*y(6) - re_r(3)*y(1)
    pdjj(1,5) = -re_f(2)*y(1)
    pdj(1,6) = re_r(2)*y(4)
    pdj(1,7) = re_f(3)*y(2)
    pdj(1,8) = 0.0_dp
    pdj(1,9) = 0.0_dp
    pdj(1,10) = 0.0_dp
    pdj(1,11) = 0.0_dp
    pdj(1,12) = 0.0_dp
    pdj(2,1) = -re_f(1)*y(2) + re_r(3)*y(4)
    pdj(2,2) = -nd_atm*re_r(4) - 9.0_dp*re_r(5)*y(2)**2*y(5) - &
      & 9.0_dp*re_r(6)*y(2)**2*y(9) - 9.0_dp*re_r(8)*y(11)*y(2)**2 - &
      & 9.0_dp*re_r(9)*y(12)*y(2)**2 - re_f(1)*y(1) - re_f(3)*y(7)
    pdj(2,3) = re_r(1)*y(4) + 3.0_dp*re_f(5)*y(8)
    pdj(2,4) = 2.0_dp*nd_atm*re_f(4)*y(4) + re_r(1)*y(3) + re_r(3)*y(1)
    pdj(2,5) = -3.0_dp*re_r(5)*y(2)**3
    pdj(2,6) = 0.0_dp
    pdj(2,7) = -re_f(3)*y(2)
    pdj(2,8) = 6.0_dp*re_f(6)*y(8) + 3.0_dp*re_f(9)*y(10) + 3.0_dp*re_f(5)*y(3)
    pdj(2,9) = -3.0_dp*re_r(6)*y(2)**3
    pdj(2,10) = 6.0_dp*re_f(8)*y(10) + 3.0_dp*re_f(9)*y(8)
    pdj(2,11) = -3.0_dp*re_r(8)*y(2)**3
    pdj(2,12) = -3.0_dp*re_r(9)*y(2)**3
    pdj(3,1) = re_f(1)*y(2)
    pdj(3,2) = 3.0_dp*re_r(5)*y(2)**2*y(5) + re_f(1)*y(1)
    pdj(3,3) = -re_r(7)*y(9) - re_r(10)*y(12) - re_r(1)*y(4) - re_f(5)*y(8)
    pdj(3,4) = -re_r(1)*y(3)
    pdj(3,5) = re_r(5)*y(2)**3 + re_f(7)*y(8) + re_f(10)*y(10)
    pdj(3,6) = 0.0_dp
    pdj(3,7) = 0.0_dp
    pdj(3,8) = re_f(7)*y(5) - re_f(5)*y(3)
    pdj(3,9) = -re_r(7)*y(3)
    pdj(3,10) = re_f(10)*y(5)
    pdj(3,11) = 0.0_dp
    pdj(3,12) = -re_r(10)*y(3)
    pdj(4,1) = re_f(1)*y(2) + re_f(2)*y(5) - re_r(3)*y(4)
    pdj(4,2) = 2.0_dp*nd_atm*re_r(4) + re_f(1)*y(1) + re_f(3)*y(7)
    pdj(4,3) = -re_r(1)*y(4)
    pdj(4,4) = -4.0_dp*nd_atm*re_f(4)*y(4) - re_r(1)*y(3) - &
      & re_r(2)*y(6) - re_r(3)*y(1)
    pdj(4,5) = re_f(2)*y(1)
    pdj(4,6) = -re_r(2)*y(4)
    pdj(4,7) = re_f(3)*y(2)
    pdj(4,8) = 0.0_dp
    pdj(4,9) = 0.0_dp
    pdj(4,10) = 0.0_dp
    pdj(4,11) = 0.0_dp
    pdj(4,12) = 0.0_dp
    pdj(5,1) = -re_f(2)*y(5)
    pdj(5,2) = -3.0_dp*re_r(5)*y(2)**2*y(5)
    pdj(5,3) = re_r(7)*y(9) + re_r(10)*y(12) + re_f(5)*y(8)
    pdj(5,4) = re_r(2)*y(6)
    pdj(5,5) = -re_r(5)*y(2)**3 - re_f(7)*y(8) - re_f(10)*y(10) - re_f(2)*y(1)
    pdj(5,6) = re_r(2)*y(4)
    pdj(5,7) = 0.0_dp
    pdj(5,8) = -re_f(7)*y(5) + re_f(5)*y(3)
    pdj(5,9) = re_r(7)*y(3)
    pdj(5,10) = -re_f(10)*y(5)
    pdj(5,11) = 0.0_dp
    pdj(5,12) = re_r(10)*y(3)
    pdj(6,1) = re_f(2)*y(5)
    pdj(6,2) = 0.0_dp
    pdj(6,3) = 0.0_dp
    pdj(6,4) = -re_r(2)*y(6)
    pdj(6,5) = re_f(2)*y(1)
    pdj(6,6) = -re_r(2)*y(4)
    pdj(6,7) = 0.0_dp
    pdj(6,8) = 0.0_dp
    pdj(6,9) = 0.0_dp
    pdj(6,10) = 0.0_dp
    pdj(6,11) = 0.0_dp
    pdj(6,12) = 0.0_dp
    pdj(7,1) = re_r(3)*y(4)
    pdj(7,2) = -re_f(3)*y(7)
    pdj(7,3) = 0.0_dp
    pdj(7,4) = re_r(3)*y(1)
    pdj(7,5) = 0.0_dp
    pdj(7,6) = 0.0_dp
    pdj(7,7) = -re_f(3)*y(2)
    pdj(7,8) = 0.0_dp
    pdj(7,9) = 0.0_dp
    pdj(7,10) = 0.0_dp
    pdj(7,11) = 0.0_dp
    pdj(7,12) = 0.0_dp
    pdj(8,1) = 0.0_dp
    pdj(8,2) = 3.0_dp*re_r(5)*y(2)**2*y(5) + 6.0_dp*re_r(6)*y(2)**2*y(9) + &
      & 3.0_dp*re_r(9)*y(12)*y(2)**2
    pdj(8,3) = re_r(7)*y(9) - re_f(5)*y(8)
    pdj(8,4) = 0.0_dp
    pdj(8,5) = re_r(5)*y(2)**3 - re_f(7)*y(8)
    pdj(8,6) = 0.0_dp
    pdj(8,7) = 0.0_dp
    pdj(8,8) = -4.0_dp*re_f(6)*y(8) - re_f(7)*y(5) - re_f(9)*y(10) - re_f(5)*y(3)
    pdj(8,9) = 2.0_dp*re_r(6)*y(2)**3 + re_r(7)*y(3)
    pdj(8,10) = -re_f(9)*y(8)
    pdj(8,11) = 0.0_dp
    pdj(8,12) = re_r(9)*y(2)**3
    pdj(9,1) = 0.0_dp
    pdj(9,2) = -3.0_dp*re_r(6)*y(2)**2*y(9)
    pdj(9,3) = -re_r(7)*y(9)
    pdj(9,4) = 0.0_dp
    pdj(9,5) = re_f(7)*y(8)
    pdj(9,6) = 0.0_dp
    pdj(9,7) = 0.0_dp
    pdj(9,8) = 2.0_dp*re_f(6)*y(8) + re_f(7)*y(5)
    pdj(9,9) = -re_r(6)*y(2)**3 - re_r(7)*y(3)
    pdj(9,10) = 0.0_dp
    pdj(9,11) = 0.0_dp
    pdj(9,12) = 0.0_dp
    pdj(10,1) = 0.0_dp
    pdj(10,2) = 6.0_dp*re_r(8)*y(11)*y(2)**2 + 3.0_dp*re_r(9)*y(12)*y(2)**2
    pdj(10,3) = re_r(10)*y(12)
    pdj(10,4) = 0.0_dp
    pdj(10,5) = -re_f(10)*y(10)
    pdj(10,6) = 0.0_dp
    pdj(10,7) = 0.0_dp
    pdj(10,8) = -re_f(9)*y(10)
    pdj(10,9) = 0.0_dp
    pdj(10,10) = -4.0_dp*re_f(8)*y(10) - re_f(9)*y(8) - re_f(10)*y(5)
    pdj(10,11) = 2.0_dp*re_r(8)*y(2)**3
    pdj(10,12) = re_r(9)*y(2)**3 + re_r(10)*y(3)
    pdj(11,1) = 0.0_dp
    pdj(11,2) = -3.0_dp*re_r(8)*y(11)*y(2)**2
    pdj(11,3) = 0.0_dp
    pdj(11,4) = 0.0_dp
    pdj(11,5) = 0.0_dp
    pdj(11,6) = 0.0_dp
    pdj(11,7) = 0.0_dp
    pdj(11,8) = 0.0_dp
    pdj(11,9) = 0.0_dp
    pdj(11,10) = 2.0_dp*re_f(8)*y(10)
    pdj(11,11) = -re_r(8)*y(2)**3
    pdj(11,12) = 0.0_dp
    pdj(12,1) = 0.0_dp
    pdj(12,2) = -3.0_dp*re_r(9)*y(12)*y(2)**2
    pdj(12,3) = -re_r(10)*y(12)
    pdj(12,4) = 0.0_dp
    pdj(12,5) = re_f(10)*y(10)
    pdj(12,6) = 0.0_dp
    pdj(12,7) = 0.0_dp
    pdj(12,8) = re_f(9)*y(10)
    pdj(12,9) = 0.0_dp
    pdj(12,10) = re_f(9)*y(8) + re_f(10)*y(5)
    pdj(12,11) = 0.0_dp
    pdj(12,12) = -re_r(9)*y(2)**3 - re_r(10)*y(3)

  end subroutine jac_NCHO

  subroutine jac_CHO(NEQ, T, Y, J, IAN, JAN, PDJ)
    implicit none
    integer, intent(in) :: NEQ, J
    real(dp), intent(in) :: T
    real(dp), dimension(NEQ), intent(in) :: Y, IAN, JAN
    real(dp), dimension(NEQ), intent(inout) :: PDJ

    pdj(1,1) = -re_f(1)*y(2) - re_f(2)*y(5) - re_r(3)*y(4)
    pdj(1,2) = -re_f(1)*y(1) + re_f(3)*y(7)
    pdj(1,3) = re_r(1)*y(4)
    pdj(1,4) = re_r(1)*y(3) + re_r(2)*y(6) - re_r(3)*y(1)
    pdj(1,5) = -re_f(2)*y(1)
    pdj(1,6) = re_r(2)*y(4)
    pdj(1,7) = re_f(3)*y(2)
    pdj(1,8) = 0.0_dp
    pdj(1,9) = 0.0_dp
    pdj(2,1) = -re_f(1)*y(2) + re_r(3)*y(4)
    pdj(2,2) = -nd_atm*re_r(4) - 9.0_dp*re_r(5)*y(2)**2*y(5) - &
    & 9.0_dp*re_r(6)*y(2)**2*y(9) - re_f(1)*y(1) - re_f(3)*y(7)
    pdj(2,3) = re_r(1)*y(4) + 3.0_dp*re_f(5)*y(8)
    pdj(2,4) = 2.0_dp*nd_atm*re_f(4)*y(4) + re_r(1)*y(3) + re_r(3)*y(1)
    pdj(2,5) = -3.0_dp*re_r(5)*y(2)**3
    pdj(2,6) = 0.0_dp
    pdj(2,7) = -re_f(3)*y(2)
    pdj(2,8) = 6.0_dp*re_f(6)*y(8) + 3.0_dp*re_f(5)*y(3)
    pdj(2,9) = -3.0_dp*re_r(6)*y(2)**3
    pdj(3,1) = re_f(1)*y(2)
    pdj(3,2) = 3.0_dp*re_r(5)*y(2)**2*y(5) + re_f(1)*y(1)
    pdj(3,3) = -re_r(7)*y(9) - re_r(1)*y(4) - re_f(5)*y(8)
    pdj(3,4) = -re_r(1)*y(3)
    pdj(3,5) = re_r(5)*y(2)**3 + re_f(7)*y(8)
    pdj(3,6) = 0.0_dp
    pdj(3,7) = 0.0_dp
    pdj(3,8) = re_f(7)*y(5) - re_f(5)*y(3)
    pdj(3,9) = -re_r(7)*y(3)
    pdj(4,1) = re_f(1)*y(2) + re_f(2)*y(5) - re_r(3)*y(4)
    pdj(4,2) = 2.0_dp*nd_atm*re_r(4) + re_f(1)*y(1) + re_f(3)*y(7)
    pdj(4,3) = -re_r(1)*y(4)
    pdj(4,4) = -4.0_dp*nd_atm*re_f(4)*y(4) - re_r(1)*y(3) - re_r(2)*y(6) - re_r(3)*y(1)
    pdj(4,5) = re_f(2)*y(1)
    pdj(4,6) = -re_r(2)*y(4)
    pdj(4,7) = re_f(3)*y(2)
    pdj(4,8) = 0.0_dp
    pdj(4,9) = 0.0_dp
    pdj(5,1) = -re_f(2)*y(5)
    pdj(5,2) = -3.0_dp*re_r(5)*y(2)**2*y(5)
    pdj(5,3) = re_r(7)*y(9) + re_f(5)*y(8)
    pdj(5,4) = re_r(2)*y(6)
    pdj(5,5) = -re_r(5)*y(2)**3 - re_f(7)*y(8) - re_f(2)*y(1)
    pdj(5,6) = re_r(2)*y(4)
    pdj(5,7) = 0.0_dp
    pdj(5,8) = -re_f(7)*y(5) + re_f(5)*y(3)
    pdj(5,9) = re_r(7)*y(3)
    pdj(6,1) = re_f(2)*y(5)
    pdj(6,2) = 0.0_dp
    pdj(6,3) = 0.0_dp
    pdj(6,4) = -re_r(2)*y(6)
    pdj(6,5) = re_f(2)*y(1)
    pdj(6,6) = -re_r(2)*y(4)
    pdj(6,7) = 0.0_dp
    pdj(6,8) = 0.0_dp
    pdj(6,9) = 0.0_dp
    pdj(7,1) = re_r(3)*y(4)
    pdj(7,2) = -re_f(3)*y(7)
    pdj(7,3) = 0.0_dp
    pdj(7,4) = re_r(3)*y(1)
    pdj(7,5) = 0.0_dp
    pdj(7,6) = 0.0_dp
    pdj(7,7) = -re_f(3)*y(2)
    pdj(7,8) = 0.0_dp
    pdj(7,9) = 0.0_dp
    pdj(8,1) = 0.0_dp
    pdj(8,2) = 3.0_dp*re_r(5)*y(2)**2*y(5) + 6.0_dp*re_r(6)*y(2)**2*y(9)
    pdj(8,3) = re_r(7)*y(9) - re_f(5)*y(8)
    pdj(8,4) = 0.0_dp
    pdj(8,5) = re_r(5)*y(2)**3 - re_f(7)*y(8)
    pdj(8,6) = 0.0_dp
    pdj(8,7) = 0.0_dp
    pdj(8,8) = -4.0_dp*re_f(6)*y(8) - re_f(7)*y(5) - re_f(5)*y(3)
    pdj(8,9) = 2.0_dp*re_r(6)*y(2)**3 + re_r(7)*y(3)
    pdj(9,1) = 0.0_dp
    pdj(9,2) = -3.0_dp*re_r(6)*y(2)**2*y(9)
    pdj(9,3) = -re_r(7)*y(9)
    pdj(9,4) = 0.0_dp
    pdj(9,5) = re_f(7)*y(8)
    pdj(9,6) = 0.0_dp
    pdj(9,7) = 0.0_dp
    pdj(9,8) = 2.0_dp*re_f(6)*y(8) + re_f(7)*y(5)
    pdj(9,9) = -re_r(6)*y(2)**3 - re_r(7)*y(3)

  end subroutine jac_CHO

  subroutine jac_HO(NEQ, T, Y, J, IAN, JAN, PDJ)
    implicit none
    integer, intent(in) :: NEQ, J
    real(dp), intent(in) :: T
    real(dp), dimension(NEQ), intent(in) :: Y, IAN, JAN
    real(dp), dimension(NEQ), intent(inout) :: PDJ

    pdj(1, 1) = -re_f(1)*y(2) - re_r(2)*y(4)
    pdj(1, 2) = -re_f(1)*y(1) + re_f(2)*y(5)
    pdj(1, 3) = re_r(1)*y(4)
    pdj(1, 4) = re_r(1)*y(3) - re_r(2)*y(1)
    pdj(1, 5) = re_f(2)*y(2)
    pdj(2, 1) = -re_f(1)*y(2) + re_r(2)*y(4)
    pdj(2, 2) = -nd_atm*re_r(3) - re_f(1)*y(1) - re_f(2)*y(5)
    pdj(2, 3) = re_r(1)*y(4)
    pdj(2, 4) = 2.0_dp*nd_atm*re_f(3)*y(4) + re_r(1)*y(3) + re_r(2)*y(1)
    pdj(2, 5) = -re_f(2)*y(2)
    pdj(3, 1) = re_f(1)*y(2)
    pdj(3, 2) = re_f(1)*y(1)
    pdj(3, 3) = -re_r(1)*y(4)
    pdj(3, 4) = -re_r(1)*y(3)
    pdj(3, 5) = 0.0_dp
    pdj(4, 1) = re_f(1)*y(2) - re_r(2)*y(4)
    pdj(4, 2) = 2.0_dp*nd_atm*re_r(3) + re_f(1)*y(1) + re_f(2)*y(5)
    pdj(4, 3) = -re_r(1)*y(4)
    pdj(4, 4) = -4.0_dp*nd_atm*re_f(3)*y(4) - re_r(1)*y(3) - re_r(2)*y(1)
    pdj(4, 5) = re_f(2)*y(2)
    pdj(5, 1) = re_r(2)*y(4)
    pdj(5, 2) = -re_f(2)*y(5)
    pdj(5, 3) = 0.0_dp
    pdj(5, 4) = re_r(2)*y(1)
    pdj(5, 5) = -re_f(2)*y(2)

  end subroutine jac_HO

end module mini_ch_i_dlsodes
