module mini_ch_i_dlsode
  use mini_ch_precision
  use mini_ch_class
  use mini_ch_chem
  implicit none

  logical, parameter :: use_stiff = .True.

  public ::  mini_ch_dlsode, RHS_update, jac_dummy, &
  &  jac_HO, jac_CHO, jac_NCHO

contains

  subroutine mini_ch_dlsode(T, P, t_end, VMR, nd_out, network)
    implicit none

    real(dp), intent(in) :: T, P, t_end
    real(dp), dimension(n_sp), intent(in) :: VMR
    character(len=200), intent(in) :: network

    real(dp), dimension(n_sp), intent(out) ::  nd_out

    integer :: ncall
    real(dp) :: P_cgs

    ! Time controls
    real(dp) :: t_begin, t_now, dt_init, t_old
    logical :: con = .False.

    ! DVODE variables
    real(dp), dimension(n_sp) :: y
    real(dp), allocatable, dimension(:) :: rwork, rtol, atol
    integer, allocatable, dimension(:) :: iwork
    integer :: itol, itask, istate, iopt, mf
    integer :: rworkdim, iworkdim
    real(dp) :: rpar
    integer :: ipar

    !! Find the number density of the atmosphere
    P_cgs = P * 10.0_dp   ! Convert pascal to dyne cm-2
    nd_atm = P_cgs/(kb*T)  ! Find initial number density [cm-3] of atmosphere

    !! Find initial number density of all species from VMR
    g_sp(:)%nd = nd_atm * VMR(:)

    ! First find the reverse reaction coefficents (Keq)
    call reverse_reactions(T, P_cgs)
    ! Find the forward, backward and net reaction rates
    call reaction_rates(T, P_cgs)

    ! -----------------------------------------
    ! ***  parameters for the DVODE solver  ***
    ! -----------------------------------------

    itask = 1
    istate = 1
    iopt = 1

    ! Method flag
    if (use_stiff .eqv. .True.) then
      ! Problem is stiff (usual)
      ! mf = 21 - full jacobian matrix with jacobian save
      mf = 21
      rworkdim = 22 +  9*n_sp + 2*n_sp**2
      iworkdim = 30 + n_sp
      allocate(rtol(n_sp), atol(n_sp), rwork(rworkdim), iwork(iworkdim))

      itol = 4
      rtol(:) = 1.0e-3_dp           ! Relative tolerances for each scalar
      atol(:) = 1.0e-99_dp          ! Absolute tolerance for each scalar (floor value)

      rwork(1) = t_end              ! Critical T value (don't integrate past time here)

      rwork(5:10) = 0.0_dp
      iwork(5:10) = 0

      rwork(5) = 1.0e-99_dp        ! Initial starting timestep (start low, will adapt in DVODE)
      rwork(6) = 0.0_dp              ! Maximum timestep (for heavy evaporation ~0.1 is required)
      rwork(7) = 0.0_dp         ! Minimum timestep
      
      iwork(5) = 0
      iwork(6) = 100000             ! Max number of internal steps

    else
      ! Problem is not too stiff (not typical)
      ! mf = 11 - full jacobian matrix with jacobian save
      mf = 11
      rworkdim = 22 + 16*n_sp + 2*n_sp**2
      iworkdim = 30 + n_sp
      allocate(rtol(n_sp), atol(n_sp), rwork(rworkdim), iwork(iworkdim))
      itol = 4
      rtol(:) = 1.0e-4_dp
      atol(:) = 1.0e-99_dp

      rwork(1) = t_end

      rwork(5:10) = 0.0_dp
      iwork(5:10) = 0

      rwork(5) = 1.0e-99_dp
      rwork(6) = t_end
      iwork(6) = 100000
    end if

    t_begin = 0.0_dp
    t_now = t_begin

    rpar = 0.0_dp
    ipar = 0

    ! Set the printing flag
    ! 0 = no printing, 1 = printing
    call xsetf(1)

    ncall = 0

    do while (t_now < t_end)

      ! Solution vector to seulex - y vector is number density of all species
      y(:) = g_sp(:)%nd

      t_old = t_now
      select case(network)
      case('HO')
        call DLSODE (RHS_update, n_sp, y, t_now, t_end, itol, rtol, atol, itask, &
        & istate, iopt, rwork, rworkdim, iwork, iworkdim, jac_HO, mf)
      case('CHO')
        call DLSODE (RHS_update, n_sp, y, t_now, t_end, itol, rtol, atol, itask, &
        & istate, iopt, rwork, rworkdim, iwork, iworkdim, jac_CHO, mf)
      case('NCHO')
        call DLSODE (RHS_update, n_sp, y, t_now, t_end, itol, rtol, atol, itask, &
        & istate, iopt, rwork, rworkdim, iwork, iworkdim, jac_NCHO, mf)
      case default
        print*, 'Invalid network provided: ', trim(network)
        stop
      end select

      if (t_now >= t_end*f_con) then
        call check_con(n_sp,g_sp(:)%nd,y(:),t_now,t_old,con)
        if (con .eqv. .True.) then
          g_sp(:)%nd = y(:)
          exit
        end if
      end if

      g_sp(:)%nd = y(:)

      ncall = ncall + 1

      if (istate == -1) then
        istate = 2
      else if (istate < -1) then
        exit
      end if

    end do

    nd_out(:) = g_sp(:)%nd

    deallocate(rtol, atol, rwork, iwork)

  end subroutine mini_ch_dlsode

  subroutine RHS_update(NEQ, time, y, f, rpar, ipar)
    implicit none

    integer, intent(in) ::  NEQ
    real(dp), intent(inout) :: time
    real(dp), dimension(NEQ), intent(inout) :: y
    real(dp), dimension(NEQ), intent(inout) :: f
    real(dp), intent(inout) :: rpar
    integer, intent(inout) :: ipar

    integer :: i, j, k
    real(dp) :: msum, msum2, frate, rrate

    nd_atm = sum(y(:))

    ! Update current number density of all species from y vector
    g_sp(:)%nd = y(:)

    ! Find the forward, backward and net reaction rates
    !call reaction_rates(T, P_cgs)

    ! Calculate the rate of change of number density for all species [cm-3/s]
    ! this is the f vector

    ! Loop through reactions add rates to the f array
    f(:) = 0.0_dp
    do i = 1, n_reac
      ! Do the forward and backward flux calculation for each speices in the reaction

      ! Find number density multiple for reactants in reaction
      msum = g_sp(re(i)%gi_re(1))%nd
      do k = 2, re(i)%n_re
         msum = msum * g_sp(re(i)%gi_re(k))%nd
      end do

      ! Find number density multiple for products in reaction
      msum2 = g_sp(re(i)%gi_pr(1))%nd
      do k = 2, re(i)%n_pr
         msum2 = msum2 * g_sp(re(i)%gi_pr(k))%nd
      end do

      if (re(i)%re_t == 3) then
        ! Mutliply both msum and msum2 by atmosphere nd if neutral body involved
        msum = msum * nd_atm
        msum2 = msum2 * nd_atm
      end if

      frate = msum * re(i)%f
      rrate = msum2 * re(i)%r

      ! Find flux for products
      f(re(i)%gi_pr(:)) = f(re(i)%gi_pr(:)) + frate - rrate

      f(re(i)%gi_re(:)) = f(re(i)%gi_re(:)) + rrate - frate

    end do

  end subroutine RHS_update

  subroutine jac_dummy (NEQ, X, Y, ML, MU, PD, NROWPD)
    integer, intent(in) :: NEQ, ML, MU, NROWPD
    real(dp), intent(in) :: X
    real(dp), dimension(NEQ), intent(in) :: Y
    real(dp), dimension(NROWPD, NEQ), intent(inout) :: PD
  end subroutine jac_dummy

  subroutine jac_NCHO(N, X, Y, ML, MU, DFY, NROWPD)
    implicit none
    integer, intent(in) :: N, ML, MU, NROWPD
    real(dp), intent(in) :: X
    real(dp), dimension(N), intent(in) :: Y
    real(dp), dimension(NROWPD, N), intent(out) :: DFY

    nd_atm = sum(y(:))

    ! Update current number density of all species from y vector
    g_sp(:)%nd = y(:)

    dfy(1,1) = -re(1)%f*y(2) - re(2)%f*y(5) - re(3)%r*y(4)
    dfy(1,2) = -re(1)%f*y(1) + re(3)%f*y(7)
    dfy(1,3) = re(1)%r*y(4)
    dfy(1,4) = re(1)%r*y(3) + re(2)%r*y(6) - re(3)%r*y(1)
    dfy(1,5) = -re(2)%f*y(1)
    dfy(1,6) = re(2)%r*y(4)
    dfy(1,7) = re(3)%f*y(2)
    dfy(1,8) = 0.0_dp
    dfy(1,9) = 0.0_dp
    dfy(1,10) = 0.0_dp
    dfy(1,11) = 0.0_dp
    dfy(1,12) = 0.0_dp
    dfy(2,1) = -re(1)%f*y(2) + re(3)%r*y(4)
    dfy(2,2) = -nd_atm*re(4)%r - 9.0_dp*re(5)%r*y(2)**2*y(5) - &
      & 9.0_dp*re(6)%r*y(2)**2*y(9) - 9.0_dp*re(8)%r*y(11)*y(2)**2 - &
      & 9.0_dp*re(9)%r*y(12)*y(2)**2 - re(1)%f*y(1) - re(3)%f*y(7)
    dfy(2,3) = re(1)%r*y(4) + 3.0_dp*re(5)%f*y(8)
    dfy(2,4) = 2.0_dp*nd_atm*re(4)%f*y(4) + re(1)%r*y(3) + re(3)%r*y(1)
    dfy(2,5) = -3.0_dp*re(5)%r*y(2)**3
    dfy(2,6) = 0.0_dp
    dfy(2,7) = -re(3)%f*y(2)
    dfy(2,8) = 6.0_dp*re(6)%f*y(8) + 3.0_dp*re(9)%f*y(10) + 3.0_dp*re(5)%f*y(3)
    dfy(2,9) = -3.0_dp*re(6)%r*y(2)**3
    dfy(2,10) = 6.0_dp*re(8)%f*y(10) + 3.0_dp*re(9)%f*y(8)
    dfy(2,11) = -3.0_dp*re(8)%r*y(2)**3
    dfy(2,12) = -3.0_dp*re(9)%r*y(2)**3
    dfy(3,1) = re(1)%f*y(2)
    dfy(3,2) = 3.0_dp*re(5)%r*y(2)**2*y(5) + re(1)%f*y(1)
    dfy(3,3) = -re(7)%r*y(9) - re(10)%r*y(12) - re(1)%r*y(4) - re(5)%f*y(8)
    dfy(3,4) = -re(1)%r*y(3)
    dfy(3,5) = re(5)%r*y(2)**3 + re(7)%f*y(8) + re(10)%f*y(10)
    dfy(3,6) = 0.0_dp
    dfy(3,7) = 0.0_dp
    dfy(3,8) = re(7)%f*y(5) - re(5)%f*y(3)
    dfy(3,9) = -re(7)%r*y(3)
    dfy(3,10) = re(10)%f*y(5)
    dfy(3,11) = 0.0_dp
    dfy(3,12) = -re(10)%r*y(3)
    dfy(4,1) = re(1)%f*y(2) + re(2)%f*y(5) - re(3)%r*y(4)
    dfy(4,2) = 2.0_dp*nd_atm*re(4)%r + re(1)%f*y(1) + re(3)%f*y(7)
    dfy(4,3) = -re(1)%r*y(4)
    dfy(4,4) = -4.0_dp*nd_atm*re(4)%f*y(4) - re(1)%r*y(3) - &
      & re(2)%r*y(6) - re(3)%r*y(1)
    dfy(4,5) = re(2)%f*y(1)
    dfy(4,6) = -re(2)%r*y(4)
    dfy(4,7) = re(3)%f*y(2)
    dfy(4,8) = 0.0_dp
    dfy(4,9) = 0.0_dp
    dfy(4,10) = 0.0_dp
    dfy(4,11) = 0.0_dp
    dfy(4,12) = 0.0_dp
    dfy(5,1) = -re(2)%f*y(5)
    dfy(5,2) = -3.0_dp*re(5)%r*y(2)**2*y(5)
    dfy(5,3) = re(7)%r*y(9) + re(10)%r*y(12) + re(5)%f*y(8)
    dfy(5,4) = re(2)%r*y(6)
    dfy(5,5) = -re(5)%r*y(2)**3 - re(7)%f*y(8) - re(10)%f*y(10) - re(2)%f*y(1)
    dfy(5,6) = re(2)%r*y(4)
    dfy(5,7) = 0.0_dp
    dfy(5,8) = -re(7)%f*y(5) + re(5)%f*y(3)
    dfy(5,9) = re(7)%r*y(3)
    dfy(5,10) = -re(10)%f*y(5)
    dfy(5,11) = 0.0_dp
    dfy(5,12) = re(10)%r*y(3)
    dfy(6,1) = re(2)%f*y(5)
    dfy(6,2) = 0.0_dp
    dfy(6,3) = 0.0_dp
    dfy(6,4) = -re(2)%r*y(6)
    dfy(6,5) = re(2)%f*y(1)
    dfy(6,6) = -re(2)%r*y(4)
    dfy(6,7) = 0.0_dp
    dfy(6,8) = 0.0_dp
    dfy(6,9) = 0.0_dp
    dfy(6,10) = 0.0_dp
    dfy(6,11) = 0.0_dp
    dfy(6,12) = 0.0_dp
    dfy(7,1) = re(3)%r*y(4)
    dfy(7,2) = -re(3)%f*y(7)
    dfy(7,3) = 0.0_dp
    dfy(7,4) = re(3)%r*y(1)
    dfy(7,5) = 0.0_dp
    dfy(7,6) = 0.0_dp
    dfy(7,7) = -re(3)%f*y(2)
    dfy(7,8) = 0.0_dp
    dfy(7,9) = 0.0_dp
    dfy(7,10) = 0.0_dp
    dfy(7,11) = 0.0_dp
    dfy(7,12) = 0.0_dp
    dfy(8,1) = 0.0_dp
    dfy(8,2) = 3.0_dp*re(5)%r*y(2)**2*y(5) + 6.0_dp*re(6)%r*y(2)**2*y(9) + &
      & 3.0_dp*re(9)%r*y(12)*y(2)**2
    dfy(8,3) = re(7)%r*y(9) - re(5)%f*y(8)
    dfy(8,4) = 0.0_dp
    dfy(8,5) = re(5)%r*y(2)**3 - re(7)%f*y(8)
    dfy(8,6) = 0.0_dp
    dfy(8,7) = 0.0_dp
    dfy(8,8) = -4.0_dp*re(6)%f*y(8) - re(7)%f*y(5) - re(9)%f*y(10) - re(5)%f*y(3)
    dfy(8,9) = 2.0_dp*re(6)%r*y(2)**3 + re(7)%r*y(3)
    dfy(8,10) = -re(9)%f*y(8)
    dfy(8,11) = 0.0_dp
    dfy(8,12) = re(9)%r*y(2)**3
    dfy(9,1) = 0.0_dp
    dfy(9,2) = -3.0_dp*re(6)%r*y(2)**2*y(9)
    dfy(9,3) = -re(7)%r*y(9)
    dfy(9,4) = 0.0_dp
    dfy(9,5) = re(7)%f*y(8)
    dfy(9,6) = 0.0_dp
    dfy(9,7) = 0.0_dp
    dfy(9,8) = 2.0_dp*re(6)%f*y(8) + re(7)%f*y(5)
    dfy(9,9) = -re(6)%r*y(2)**3 - re(7)%r*y(3)
    dfy(9,10) = 0.0_dp
    dfy(9,11) = 0.0_dp
    dfy(9,12) = 0.0_dp
    dfy(10,1) = 0.0_dp
    dfy(10,2) = 6.0_dp*re(8)%r*y(11)*y(2)**2 + 3.0_dp*re(9)%r*y(12)*y(2)**2
    dfy(10,3) = re(10)%r*y(12)
    dfy(10,4) = 0.0_dp
    dfy(10,5) = -re(10)%f*y(10)
    dfy(10,6) = 0.0_dp
    dfy(10,7) = 0.0_dp
    dfy(10,8) = -re(9)%f*y(10)
    dfy(10,9) = 0.0_dp
    dfy(10,10) = -4.0_dp*re(8)%f*y(10) - re(9)%f*y(8) - re(10)%f*y(5)
    dfy(10,11) = 2.0_dp*re(8)%r*y(2)**3
    dfy(10,12) = re(9)%r*y(2)**3 + re(10)%r*y(3)
    dfy(11,1) = 0.0_dp
    dfy(11,2) = -3.0_dp*re(8)%r*y(11)*y(2)**2
    dfy(11,3) = 0.0_dp
    dfy(11,4) = 0.0_dp
    dfy(11,5) = 0.0_dp
    dfy(11,6) = 0.0_dp
    dfy(11,7) = 0.0_dp
    dfy(11,8) = 0.0_dp
    dfy(11,9) = 0.0_dp
    dfy(11,10) = 2.0_dp*re(8)%f*y(10)
    dfy(11,11) = -re(8)%r*y(2)**3
    dfy(11,12) = 0.0_dp
    dfy(12,1) = 0.0_dp
    dfy(12,2) = -3.0_dp*re(9)%r*y(12)*y(2)**2
    dfy(12,3) = -re(10)%r*y(12)
    dfy(12,4) = 0.0_dp
    dfy(12,5) = re(10)%f*y(10)
    dfy(12,6) = 0.0_dp
    dfy(12,7) = 0.0_dp
    dfy(12,8) = re(9)%f*y(10)
    dfy(12,9) = 0.0_dp
    dfy(12,10) = re(9)%f*y(8) + re(10)%f*y(5)
    dfy(12,11) = 0.0_dp
    dfy(12,12) = -re(9)%r*y(2)**3 - re(10)%r*y(3)

  end subroutine jac_NCHO

  subroutine jac_CHO(N, X, Y, ML, MU, DFY, NROWPD)
    implicit none
    integer, intent(in) :: N, ML, MU, NROWPD
    real(dp), intent(in) :: X
    real(dp), dimension(N), intent(in) :: Y
    real(dp), dimension(NROWPD, N), intent(out) :: DFY

    nd_atm = sum(y(:))

    ! Update current number density of all species from y vector
    g_sp(:)%nd = y(:)

    ! Find the forward, backward and net reaction rates
    !call reaction_rates(T, P_cgs)

    dfy(1,1) = -re(1)%f*y(2) - re(2)%f*y(5) - re(3)%r*y(4)
    dfy(1,2) = -re(1)%f*y(1) + re(3)%f*y(7)
    dfy(1,3) = re(1)%r*y(4)
    dfy(1,4) = re(1)%r*y(3) + re(2)%r*y(6) - re(3)%r*y(1)
    dfy(1,5) = -re(2)%f*y(1)
    dfy(1,6) = re(2)%r*y(4)
    dfy(1,7) = re(3)%f*y(2)
    dfy(1,8) = 0.0_dp
    dfy(1,9) = 0.0_dp
    dfy(2,1) = -re(1)%f*y(2) + re(3)%r*y(4)
    dfy(2,2) = -nd_atm*re(4)%r - 9.0_dp*re(5)%r*y(2)**2*y(5) - &
    & 9.0_dp*re(6)%r*y(2)**2*y(9) - re(1)%f*y(1) - re(3)%f*y(7)
    dfy(2,3) = re(1)%r*y(4) + 3.0_dp*re(5)%f*y(8)
    dfy(2,4) = 2.0_dp*nd_atm*re(4)%f*y(4) + re(1)%r*y(3) + re(3)%r*y(1)
    dfy(2,5) = -3.0_dp*re(5)%r*y(2)**3
    dfy(2,6) = 0.0_dp
    dfy(2,7) = -re(3)%f*y(2)
    dfy(2,8) = 6.0_dp*re(6)%f*y(8) + 3.0_dp*re(5)%f*y(3)
    dfy(2,9) = -3.0_dp*re(6)%r*y(2)**3
    dfy(3,1) = re(1)%f*y(2)
    dfy(3,2) = 3.0_dp*re(5)%r*y(2)**2*y(5) + re(1)%f*y(1)
    dfy(3,3) = -re(7)%r*y(9) - re(1)%r*y(4) - re(5)%f*y(8)
    dfy(3,4) = -re(1)%r*y(3)
    dfy(3,5) = re(5)%r*y(2)**3 + re(7)%f*y(8)
    dfy(3,6) = 0.0_dp
    dfy(3,7) = 0.0_dp
    dfy(3,8) = re(7)%f*y(5) - re(5)%f*y(3)
    dfy(3,9) = -re(7)%r*y(3)
    dfy(4,1) = re(1)%f*y(2) + re(2)%f*y(5) - re(3)%r*y(4)
    dfy(4,2) = 2.0_dp*nd_atm*re(4)%r + re(1)%f*y(1) + re(3)%f*y(7)
    dfy(4,3) = -re(1)%r*y(4)
    dfy(4,4) = -4.0_dp*nd_atm*re(4)%f*y(4) - re(1)%r*y(3) - re(2)%r*y(6) - re(3)%r*y(1)
    dfy(4,5) = re(2)%f*y(1)
    dfy(4,6) = -re(2)%r*y(4)
    dfy(4,7) = re(3)%f*y(2)
    dfy(4,8) = 0.0_dp
    dfy(4,9) = 0.0_dp
    dfy(5,1) = -re(2)%f*y(5)
    dfy(5,2) = -3.0_dp*re(5)%r*y(2)**2*y(5)
    dfy(5,3) = re(7)%r*y(9) + re(5)%f*y(8)
    dfy(5,4) = re(2)%r*y(6)
    dfy(5,5) = -re(5)%r*y(2)**3 - re(7)%f*y(8) - re(2)%f*y(1)
    dfy(5,6) = re(2)%r*y(4)
    dfy(5,7) = 0.0_dp
    dfy(5,8) = -re(7)%f*y(5) + re(5)%f*y(3)
    dfy(5,9) = re(7)%r*y(3)
    dfy(6,1) = re(2)%f*y(5)
    dfy(6,2) = 0.0_dp
    dfy(6,3) = 0.0_dp
    dfy(6,4) = -re(2)%r*y(6)
    dfy(6,5) = re(2)%f*y(1)
    dfy(6,6) = -re(2)%r*y(4)
    dfy(6,7) = 0.0_dp
    dfy(6,8) = 0.0_dp
    dfy(6,9) = 0.0_dp
    dfy(7,1) = re(3)%r*y(4)
    dfy(7,2) = -re(3)%f*y(7)
    dfy(7,3) = 0.0_dp
    dfy(7,4) = re(3)%r*y(1)
    dfy(7,5) = 0.0_dp
    dfy(7,6) = 0.0_dp
    dfy(7,7) = -re(3)%f*y(2)
    dfy(7,8) = 0.0_dp
    dfy(7,9) = 0.0_dp
    dfy(8,1) = 0.0_dp
    dfy(8,2) = 3.0_dp*re(5)%r*y(2)**2*y(5) + 6.0_dp*re(6)%r*y(2)**2*y(9)
    dfy(8,3) = re(7)%r*y(9) - re(5)%f*y(8)
    dfy(8,4) = 0.0_dp
    dfy(8,5) = re(5)%r*y(2)**3 - re(7)%f*y(8)
    dfy(8,6) = 0.0_dp
    dfy(8,7) = 0.0_dp
    dfy(8,8) = -4.0_dp*re(6)%f*y(8) - re(7)%f*y(5) - re(5)%f*y(3)
    dfy(8,9) = 2.0_dp*re(6)%r*y(2)**3 + re(7)%r*y(3)
    dfy(9,1) = 0.0_dp
    dfy(9,2) = -3.0_dp*re(6)%r*y(2)**2*y(9)
    dfy(9,3) = -re(7)%r*y(9)
    dfy(9,4) = 0.0_dp
    dfy(9,5) = re(7)%f*y(8)
    dfy(9,6) = 0.0_dp
    dfy(9,7) = 0.0_dp
    dfy(9,8) = 2.0_dp*re(6)%f*y(8) + re(7)%f*y(5)
    dfy(9,9) = -re(6)%r*y(2)**3 - re(7)%r*y(3)

  end subroutine jac_CHO

  subroutine jac_HO(N, X, Y, ML, MU, DFY, NROWPD)
    implicit none
    integer, intent(in) :: N, ML, MU, NROWPD
    real(dp), intent(in) :: X
    real(dp), dimension(N), intent(in) :: Y
    real(dp), dimension(NROWPD, N), intent(out) :: DFY

    nd_atm = sum(y(:))

    ! Update current number density of all species from y vector
    g_sp(:)%nd = y(:)

    ! Find the forward, backward and net reaction rates
    !call reaction_rates(T, P_cgs)

    dfy(1, 1) = -re(1)%f*y(2) - re(2)%r*y(4)
    dfy(1, 2) = -re(1)%f*y(1) + re(2)%f*y(5)
    dfy(1, 3) = re(1)%r*y(4)
    dfy(1, 4) = re(1)%r*y(3) - re(2)%r*y(1)
    dfy(1, 5) = re(2)%f*y(2)
    dfy(2, 1) = -re(1)%f*y(2) + re(2)%r*y(4)
    dfy(2, 2) = -nd_atm*re(3)%r - re(1)%f*y(1) - re(2)%f*y(5)
    dfy(2, 3) = re(1)%r*y(4)
    dfy(2, 4) = 2.0_dp*nd_atm*re(3)%f*y(4) + re(1)%r*y(3) + re(2)%r*y(1)
    dfy(2, 5) = -re(2)%f*y(2)
    dfy(3, 1) = re(1)%f*y(2)
    dfy(3, 2) = re(1)%f*y(1)
    dfy(3, 3) = -re(1)%r*y(4)
    dfy(3, 4) = -re(1)%r*y(3)
    dfy(3, 5) = 0.0_dp
    dfy(4, 1) = re(1)%f*y(2) - re(2)%r*y(4)
    dfy(4, 2) = 2.0_dp*nd_atm*re(3)%r + re(1)%f*y(1) + re(2)%f*y(5)
    dfy(4, 3) = -re(1)%r*y(4)
    dfy(4, 4) = -4.0_dp*nd_atm*re(3)%f*y(4) - re(1)%r*y(3) - re(2)%r*y(1)
    dfy(4, 5) = re(2)%f*y(2)
    dfy(5, 1) = re(2)%r*y(4)
    dfy(5, 2) = -re(2)%f*y(5)
    dfy(5, 3) = 0.0_dp
    dfy(5, 4) = re(2)%r*y(1)
    dfy(5, 5) = -re(2)%f*y(2)

  end subroutine jac_HO

end module mini_ch_i_dlsode
