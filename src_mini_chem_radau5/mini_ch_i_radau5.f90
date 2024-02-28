module mini_ch_i_radau5
  use mini_ch_precision
  use mini_ch_class
  use mini_ch_chem, only: reaction_rates, reverse_reactions, check_con
  implicit none

  real(dp) :: nd_atm

  private
  public :: mini_ch_radau5, RHS_update, jac_dummy, mas_dummy, solout &
  &, jac_HO, jac_CHO, jac_NCHO

contains

  subroutine mini_ch_radau5(T_in, P_in, t_end, VMR, network)
    implicit none

    real(dp), intent(in) :: T_in, P_in, t_end
    real(dp), dimension(n_sp), intent(inout) :: VMR
    character(len=200), intent(in) :: network

    integer :: ncall
    real(dp) :: P_cgs

    ! Time controls
    real(dp) :: t_now,  t_old, dt_init

    ! seulex variables
    real(dp), dimension(n_sp) :: y, y_old
    real(dp), allocatable, dimension(:) :: rwork
    integer, allocatable, dimension(:) :: iwork
    real(dp) :: rtol, atol
    real(dp) :: rpar
    integer :: itol, ijac, mljac, mujac, imas, mlmas, mumas, iout, lrwork, liwork, ipar, idid

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
    ! ***  parameters for the RADAU5-solver  ***
    ! -----------------------------------------

    rtol = 1.0e-3_dp
    atol = 1.0e-99_dp
    itol = 0
    ijac = 1
    mljac = n_sp
    mujac = n_sp
    imas = 0
    mlmas = n_sp
    mumas = n_sp
    iout = 0
    idid = 0

    ! Real work array
    lrwork = n_sp*(n_sp*3 + 12 + 8) + 4*12 + 20
    allocate(rwork(lrwork))
    rwork(:) = 0.0_dp

    rwork(1) = 1.0e-16_dp ! Rounding unit - default 1e-16
    rwork(2) = 0.0_dp  ! Safety factor
    rwork(3) = 1.0e-4_dp  ! Jacobian recompuation rate
    rwork(4) = 0.0_dp ! Stopping critera for Netwon step
    rwork(5) = 0.0_dp  ! 
    rwork(6) = 0.0_dp  ! 
    rwork(7) = 0.0_dp  ! Max time-step
    rwork(8) = 0.0_dp ! 
    rwork(9) = 0.0_dp !

    ! Integer work array
    liwork = 2*n_sp + 12 + 20
    allocate(iwork(liwork))
    iwork(:) = 0

    iwork(1) = 0 ! Hessenberg form?
    iwork(2) = 0 ! Default step numbers - 0 = 100000
    iwork(3) = 0 ! Max columns in extrapolation table - default 12
    iwork(4) = 0 ! Switch for step size sequence - defualt 0 = 2
    iwork(5) = 0 ! Lambda parameter for dense output
    iwork(6) = 0 ! Number of components for which dense output is required
    iwork(7) = 0
    iwork(8) = 0
    iwork(9) = 0
    iwork(10) = 0

    rpar = 0.0_dp
    ipar = 0

    t_now = 0.0_dp
    dt_init = 0.0_dp

    ncall = 0

    do while((t_now < t_end))

      y_old(:) = y(:)
      t_old = t_now

      select case(network)
      case('HO')
        call RADAU5(n_sp,RHS_update,t_now,y,t_end,dt_init, &
          &                  rtol,atol,itol, &
          &                  jac_HO,ijac,mljac,mujac, &
          &                  mas_dummy,imas,mlmas,mumas, &
          &                  solout,iout, &
          &                  rwork,lrwork,iwork,liwork,rpar,ipar,idid)
      case('CHO')
        call RADAU5(n_sp,RHS_update,t_now,y,t_end,dt_init, &
          &                  rtol,atol,itol, &
          &                  jac_CHO,ijac,mljac,mujac, &
          &                  mas_dummy,imas,mlmas,mumas, &
          &                  solout,iout, &
          &                  rwork,lrwork,iwork,liwork,rpar,ipar,idid)
      case('NCHO')
        call RADAU5(n_sp,RHS_update,t_now,y,t_end,dt_init, &
          &                  rtol,atol,itol, &
          &                  jac_NCHO,ijac,mljac,mujac, &
          &                  mas_dummy,imas,mlmas,mumas, &
          &                  solout,iout, &
          &                  rwork,lrwork,iwork,liwork,rpar,ipar,idid)
      case default
        print*, 'Invalid network provided: ', trim(network)
        stop
      end select

      !call check_con(n_sp,y(:),y_old(:),t_now,t_old,con)
      !if (con .eqv. .True.) then
      !  exit
      !end if

      ncall = ncall + 1

    end do

    VMR(:) = y(:)/nd_atm

    deallocate(rwork, iwork, Keq, re_r, re_f)

  end subroutine mini_ch_radau5

  subroutine RHS_update(NEQ, time, y, f, rpar, ipar)
    implicit none

    integer, intent(in) ::  NEQ
    real(dp), intent(inout) :: time
    real(dp), dimension(NEQ), intent(inout) :: y
    real(dp), dimension(NEQ), intent(inout) :: f
    real(dp), intent(inout) :: rpar
    integer, intent(inout) :: ipar

    integer :: i, k, j
    real(dp) :: msum, msum2, frate, rrate
    real(dp), dimension(n_reac) :: net_pr, net_re
    real(dp), dimension(NEQ) :: f_pr, f_re, t_pr, t_re
    real(dp), dimension(NEQ) :: c_pr, c_re

    ! Calculate the rate of change of number density for all species [cm-3/s]
    ! this is the f vector
    f_pr(:) = 0.0_dp
    c_pr(:) = 0.0_dp
    f_re(:) = 0.0_dp
    c_re(:) = 0.0_dp

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

      net_pr(i) = frate - rrate
      net_re(i) = -net_pr(i)

      !! Perform the Kahan-Babushka-Neumaier compensation sum algorithm
      ! This is slightly slower than peicewise addition for small timesteps, but faster for larger timesteps, 
      ! and more general (should work for all networks)

      !! Add the product rates
      do j = 1, re(i)%n_pr
        t_pr(re(i)%gi_pr(j)) = f_pr(re(i)%gi_pr(j)) + net_pr(i)
        if (abs(f_pr(re(i)%gi_pr(j))) >= abs(net_pr(i))) then
          c_pr(re(i)%gi_pr(j)) = c_pr(re(i)%gi_pr(j)) + (f_pr(re(i)%gi_pr(j)) - t_pr(re(i)%gi_pr(j))) + net_pr(i)
        else
          c_pr(re(i)%gi_pr(j)) = c_pr(re(i)%gi_pr(j)) + (net_pr(i) - t_pr(re(i)%gi_pr(j))) + f_pr(re(i)%gi_pr(j))
        end if
        f_pr(re(i)%gi_pr(j)) = t_pr(re(i)%gi_pr(j))
      end do
      f_pr(re(i)%gi_pr(:)) =  f_pr(re(i)%gi_pr(:)) + c_pr(re(i)%gi_pr(:))

      !! Add the reactant rates
      do j = 1, re(i)%n_re
        t_re(re(i)%gi_re(j)) = f_re(re(i)%gi_re(j)) + net_re(i)
        if (abs(f_re(re(i)%gi_re(j))) >= abs(net_re(i))) then
          c_re(re(i)%gi_re(j)) = c_re(re(i)%gi_re(j)) + (f_re(re(i)%gi_re(j)) - t_re(re(i)%gi_re(j))) + net_re(i)
        else
          c_re(re(i)%gi_re(j)) = c_re(re(i)%gi_re(j)) + (net_re(i) - t_re(re(i)%gi_re(j))) + f_re(re(i)%gi_re(j))
        end if
        f_re(re(i)%gi_re(j)) = t_re(re(i)%gi_re(j))
      end do
      f_re(re(i)%gi_re(:)) =  f_re(re(i)%gi_re(:)) + c_re(re(i)%gi_re(:))

    end do

    !! Sum product and reactant rates to get net rate for species
    f(:) = f_pr(:) + f_re(:)

  end subroutine RHS_update

  subroutine jac_dummy(N,X,Y,DFY,LDFY,RPAR,IPAR)
    integer :: N,LDFY,IPAR
    double precision :: X,Y(N),DFY(LDFY,N),RPAR
  end subroutine jac_dummy

  subroutine jac_NCHO(N,X,Y,DFY,LDFY,RPAR,IPAR)
    implicit none
    integer, intent(in) :: N, LDFY, ipar
    real(dp), intent(in) :: X, RPAR
    real(dp), dimension(N), intent(in) :: Y
    real(dp), dimension(LDFY, N),intent(out) :: DFY

    dfy(1,1) = -re_f(1)*y(2) - re_f(2)*y(5) - re_r(3)*y(4)
    dfy(1,2) = -re_f(1)*y(1) + re_f(3)*y(7)
    dfy(1,3) = re_r(1)*y(4)
    dfy(1,4) = re_r(1)*y(3) + re_r(2)*y(6) - re_r(3)*y(1)
    dfy(1,5) = -re_f(2)*y(1)
    dfy(1,6) = re_r(2)*y(4)
    dfy(1,7) = re_f(3)*y(2)
    dfy(1,8) = 0.0_dp
    dfy(1,9) = 0.0_dp
    dfy(1,10) = 0.0_dp
    dfy(1,11) = 0.0_dp
    dfy(1,12) = 0.0_dp
    dfy(2,1) = -re_f(1)*y(2) + re_r(3)*y(4)
    dfy(2,2) = -nd_atm*re_r(4) - 9.0_dp*re_r(5)*y(2)**2*y(5) - &
      & 9.0_dp*re_r(6)*y(2)**2*y(9) - 9.0_dp*re_r(8)*y(11)*y(2)**2 - &
      & 9.0_dp*re_r(9)*y(12)*y(2)**2 - re_f(1)*y(1) - re_f(3)*y(7)
    dfy(2,3) = re_r(1)*y(4) + 3.0_dp*re_f(5)*y(8)
    dfy(2,4) = 2.0_dp*nd_atm*re_f(4)*y(4) + re_r(1)*y(3) + re_r(3)*y(1)
    dfy(2,5) = -3.0_dp*re_r(5)*y(2)**3
    dfy(2,6) = 0.0_dp
    dfy(2,7) = -re_f(3)*y(2)
    dfy(2,8) = 6.0_dp*re_f(6)*y(8) + 3.0_dp*re_f(9)*y(10) + 3.0_dp*re_f(5)*y(3)
    dfy(2,9) = -3.0_dp*re_r(6)*y(2)**3
    dfy(2,10) = 6.0_dp*re_f(8)*y(10) + 3.0_dp*re_f(9)*y(8)
    dfy(2,11) = -3.0_dp*re_r(8)*y(2)**3
    dfy(2,12) = -3.0_dp*re_r(9)*y(2)**3
    dfy(3,1) = re_f(1)*y(2)
    dfy(3,2) = 3.0_dp*re_r(5)*y(2)**2*y(5) + re_f(1)*y(1)
    dfy(3,3) = -re_r(7)*y(9) - re_r(10)*y(12) - re_r(1)*y(4) - re_f(5)*y(8)
    dfy(3,4) = -re_r(1)*y(3)
    dfy(3,5) = re_r(5)*y(2)**3 + re_f(7)*y(8) + re_f(10)*y(10)
    dfy(3,6) = 0.0_dp
    dfy(3,7) = 0.0_dp
    dfy(3,8) = re_f(7)*y(5) - re_f(5)*y(3)
    dfy(3,9) = -re_r(7)*y(3)
    dfy(3,10) = re_f(10)*y(5)
    dfy(3,11) = 0.0_dp
    dfy(3,12) = -re_r(10)*y(3)
    dfy(4,1) = re_f(1)*y(2) + re_f(2)*y(5) - re_r(3)*y(4)
    dfy(4,2) = 2.0_dp*nd_atm*re_r(4) + re_f(1)*y(1) + re_f(3)*y(7)
    dfy(4,3) = -re_r(1)*y(4)
    dfy(4,4) = -4.0_dp*nd_atm*re_f(4)*y(4) - re_r(1)*y(3) - &
      & re_r(2)*y(6) - re_r(3)*y(1)
    dfy(4,5) = re_f(2)*y(1)
    dfy(4,6) = -re_r(2)*y(4)
    dfy(4,7) = re_f(3)*y(2)
    dfy(4,8) = 0.0_dp
    dfy(4,9) = 0.0_dp
    dfy(4,10) = 0.0_dp
    dfy(4,11) = 0.0_dp
    dfy(4,12) = 0.0_dp
    dfy(5,1) = -re_f(2)*y(5)
    dfy(5,2) = -3.0_dp*re_r(5)*y(2)**2*y(5)
    dfy(5,3) = re_r(7)*y(9) + re_r(10)*y(12) + re_f(5)*y(8)
    dfy(5,4) = re_r(2)*y(6)
    dfy(5,5) = -re_r(5)*y(2)**3 - re_f(7)*y(8) - re_f(10)*y(10) - re_f(2)*y(1)
    dfy(5,6) = re_r(2)*y(4)
    dfy(5,7) = 0.0_dp
    dfy(5,8) = -re_f(7)*y(5) + re_f(5)*y(3)
    dfy(5,9) = re_r(7)*y(3)
    dfy(5,10) = -re_f(10)*y(5)
    dfy(5,11) = 0.0_dp
    dfy(5,12) = re_r(10)*y(3)
    dfy(6,1) = re_f(2)*y(5)
    dfy(6,2) = 0.0_dp
    dfy(6,3) = 0.0_dp
    dfy(6,4) = -re_r(2)*y(6)
    dfy(6,5) = re_f(2)*y(1)
    dfy(6,6) = -re_r(2)*y(4)
    dfy(6,7) = 0.0_dp
    dfy(6,8) = 0.0_dp
    dfy(6,9) = 0.0_dp
    dfy(6,10) = 0.0_dp
    dfy(6,11) = 0.0_dp
    dfy(6,12) = 0.0_dp
    dfy(7,1) = re_r(3)*y(4)
    dfy(7,2) = -re_f(3)*y(7)
    dfy(7,3) = 0.0_dp
    dfy(7,4) = re_r(3)*y(1)
    dfy(7,5) = 0.0_dp
    dfy(7,6) = 0.0_dp
    dfy(7,7) = -re_f(3)*y(2)
    dfy(7,8) = 0.0_dp
    dfy(7,9) = 0.0_dp
    dfy(7,10) = 0.0_dp
    dfy(7,11) = 0.0_dp
    dfy(7,12) = 0.0_dp
    dfy(8,1) = 0.0_dp
    dfy(8,2) = 3.0_dp*re_r(5)*y(2)**2*y(5) + 6.0_dp*re_r(6)*y(2)**2*y(9) + &
      & 3.0_dp*re_r(9)*y(12)*y(2)**2
    dfy(8,3) = re_r(7)*y(9) - re_f(5)*y(8)
    dfy(8,4) = 0.0_dp
    dfy(8,5) = re_r(5)*y(2)**3 - re_f(7)*y(8)
    dfy(8,6) = 0.0_dp
    dfy(8,7) = 0.0_dp
    dfy(8,8) = -4.0_dp*re_f(6)*y(8) - re_f(7)*y(5) - re_f(9)*y(10) - re_f(5)*y(3)
    dfy(8,9) = 2.0_dp*re_r(6)*y(2)**3 + re_r(7)*y(3)
    dfy(8,10) = -re_f(9)*y(8)
    dfy(8,11) = 0.0_dp
    dfy(8,12) = re_r(9)*y(2)**3
    dfy(9,1) = 0.0_dp
    dfy(9,2) = -3.0_dp*re_r(6)*y(2)**2*y(9)
    dfy(9,3) = -re_r(7)*y(9)
    dfy(9,4) = 0.0_dp
    dfy(9,5) = re_f(7)*y(8)
    dfy(9,6) = 0.0_dp
    dfy(9,7) = 0.0_dp
    dfy(9,8) = 2.0_dp*re_f(6)*y(8) + re_f(7)*y(5)
    dfy(9,9) = -re_r(6)*y(2)**3 - re_r(7)*y(3)
    dfy(9,10) = 0.0_dp
    dfy(9,11) = 0.0_dp
    dfy(9,12) = 0.0_dp
    dfy(10,1) = 0.0_dp
    dfy(10,2) = 6.0_dp*re_r(8)*y(11)*y(2)**2 + 3.0_dp*re_r(9)*y(12)*y(2)**2
    dfy(10,3) = re_r(10)*y(12)
    dfy(10,4) = 0.0_dp
    dfy(10,5) = -re_f(10)*y(10)
    dfy(10,6) = 0.0_dp
    dfy(10,7) = 0.0_dp
    dfy(10,8) = -re_f(9)*y(10)
    dfy(10,9) = 0.0_dp
    dfy(10,10) = -4.0_dp*re_f(8)*y(10) - re_f(9)*y(8) - re_f(10)*y(5)
    dfy(10,11) = 2.0_dp*re_r(8)*y(2)**3
    dfy(10,12) = re_r(9)*y(2)**3 + re_r(10)*y(3)
    dfy(11,1) = 0.0_dp
    dfy(11,2) = -3.0_dp*re_r(8)*y(11)*y(2)**2
    dfy(11,3) = 0.0_dp
    dfy(11,4) = 0.0_dp
    dfy(11,5) = 0.0_dp
    dfy(11,6) = 0.0_dp
    dfy(11,7) = 0.0_dp
    dfy(11,8) = 0.0_dp
    dfy(11,9) = 0.0_dp
    dfy(11,10) = 2.0_dp*re_f(8)*y(10)
    dfy(11,11) = -re_r(8)*y(2)**3
    dfy(11,12) = 0.0_dp
    dfy(12,1) = 0.0_dp
    dfy(12,2) = -3.0_dp*re_r(9)*y(12)*y(2)**2
    dfy(12,3) = -re_r(10)*y(12)
    dfy(12,4) = 0.0_dp
    dfy(12,5) = re_f(10)*y(10)
    dfy(12,6) = 0.0_dp
    dfy(12,7) = 0.0_dp
    dfy(12,8) = re_f(9)*y(10)
    dfy(12,9) = 0.0_dp
    dfy(12,10) = re_f(9)*y(8) + re_f(10)*y(5)
    dfy(12,11) = 0.0_dp
    dfy(12,12) = -re_r(9)*y(2)**3 - re_r(10)*y(3)

  end subroutine jac_NCHO

  subroutine jac_CHO(N,X,Y,DFY,LDFY,RPAR,IPAR)
    implicit none
    integer, intent(in) :: N, LDFY, ipar
    real(dp), intent(in) :: X, RPAR
    real(dp), dimension(N), intent(in) :: Y
    real(dp), dimension(LDFY, N),intent(out) :: DFY

    dfy(1,1) = -re_f(1)*y(2) - re_f(2)*y(5) - re_r(3)*y(4)
    dfy(1,2) = -re_f(1)*y(1) + re_f(3)*y(7)
    dfy(1,3) = re_r(1)*y(4)
    dfy(1,4) = re_r(1)*y(3) + re_r(2)*y(6) - re_r(3)*y(1)
    dfy(1,5) = -re_f(2)*y(1)
    dfy(1,6) = re_r(2)*y(4)
    dfy(1,7) = re_f(3)*y(2)
    dfy(1,8) = 0.0_dp
    dfy(1,9) = 0.0_dp
    dfy(2,1) = -re_f(1)*y(2) + re_r(3)*y(4)
    dfy(2,2) = -nd_atm*re_r(4) - 9.0_dp*re_r(5)*y(2)**2*y(5) - &
    & 9.0_dp*re_r(6)*y(2)**2*y(9) - re_f(1)*y(1) - re_f(3)*y(7)
    dfy(2,3) = re_r(1)*y(4) + 3.0_dp*re_f(5)*y(8)
    dfy(2,4) = 2.0_dp*nd_atm*re_f(4)*y(4) + re_r(1)*y(3) + re_r(3)*y(1)
    dfy(2,5) = -3.0_dp*re_r(5)*y(2)**3
    dfy(2,6) = 0.0_dp
    dfy(2,7) = -re_f(3)*y(2)
    dfy(2,8) = 6.0_dp*re_f(6)*y(8) + 3.0_dp*re_f(5)*y(3)
    dfy(2,9) = -3.0_dp*re_r(6)*y(2)**3
    dfy(3,1) = re_f(1)*y(2)
    dfy(3,2) = 3.0_dp*re_r(5)*y(2)**2*y(5) + re_f(1)*y(1)
    dfy(3,3) = -re_r(7)*y(9) - re_r(1)*y(4) - re_f(5)*y(8)
    dfy(3,4) = -re_r(1)*y(3)
    dfy(3,5) = re_r(5)*y(2)**3 + re_f(7)*y(8)
    dfy(3,6) = 0.0_dp
    dfy(3,7) = 0.0_dp
    dfy(3,8) = re_f(7)*y(5) - re_f(5)*y(3)
    dfy(3,9) = -re_r(7)*y(3)
    dfy(4,1) = re_f(1)*y(2) + re_f(2)*y(5) - re_r(3)*y(4)
    dfy(4,2) = 2.0_dp*nd_atm*re_r(4) + re_f(1)*y(1) + re_f(3)*y(7)
    dfy(4,3) = -re_r(1)*y(4)
    dfy(4,4) = -4.0_dp*nd_atm*re_f(4)*y(4) - re_r(1)*y(3) - re_r(2)*y(6) - re_r(3)*y(1)
    dfy(4,5) = re_f(2)*y(1)
    dfy(4,6) = -re_r(2)*y(4)
    dfy(4,7) = re_f(3)*y(2)
    dfy(4,8) = 0.0_dp
    dfy(4,9) = 0.0_dp
    dfy(5,1) = -re_f(2)*y(5)
    dfy(5,2) = -3.0_dp*re_r(5)*y(2)**2*y(5)
    dfy(5,3) = re_r(7)*y(9) + re_f(5)*y(8)
    dfy(5,4) = re_r(2)*y(6)
    dfy(5,5) = -re_r(5)*y(2)**3 - re_f(7)*y(8) - re_f(2)*y(1)
    dfy(5,6) = re_r(2)*y(4)
    dfy(5,7) = 0.0_dp
    dfy(5,8) = -re_f(7)*y(5) + re_f(5)*y(3)
    dfy(5,9) = re_r(7)*y(3)
    dfy(6,1) = re_f(2)*y(5)
    dfy(6,2) = 0.0_dp
    dfy(6,3) = 0.0_dp
    dfy(6,4) = -re_r(2)*y(6)
    dfy(6,5) = re_f(2)*y(1)
    dfy(6,6) = -re_r(2)*y(4)
    dfy(6,7) = 0.0_dp
    dfy(6,8) = 0.0_dp
    dfy(6,9) = 0.0_dp
    dfy(7,1) = re_r(3)*y(4)
    dfy(7,2) = -re_f(3)*y(7)
    dfy(7,3) = 0.0_dp
    dfy(7,4) = re_r(3)*y(1)
    dfy(7,5) = 0.0_dp
    dfy(7,6) = 0.0_dp
    dfy(7,7) = -re_f(3)*y(2)
    dfy(7,8) = 0.0_dp
    dfy(7,9) = 0.0_dp
    dfy(8,1) = 0.0_dp
    dfy(8,2) = 3.0_dp*re_r(5)*y(2)**2*y(5) + 6.0_dp*re_r(6)*y(2)**2*y(9)
    dfy(8,3) = re_r(7)*y(9) - re_f(5)*y(8)
    dfy(8,4) = 0.0_dp
    dfy(8,5) = re_r(5)*y(2)**3 - re_f(7)*y(8)
    dfy(8,6) = 0.0_dp
    dfy(8,7) = 0.0_dp
    dfy(8,8) = -4.0_dp*re_f(6)*y(8) - re_f(7)*y(5) - re_f(5)*y(3)
    dfy(8,9) = 2.0_dp*re_r(6)*y(2)**3 + re_r(7)*y(3)
    dfy(9,1) = 0.0_dp
    dfy(9,2) = -3.0_dp*re_r(6)*y(2)**2*y(9)
    dfy(9,3) = -re_r(7)*y(9)
    dfy(9,4) = 0.0_dp
    dfy(9,5) = re_f(7)*y(8)
    dfy(9,6) = 0.0_dp
    dfy(9,7) = 0.0_dp
    dfy(9,8) = 2.0_dp*re_f(6)*y(8) + re_f(7)*y(5)
    dfy(9,9) = -re_r(6)*y(2)**3 - re_r(7)*y(3)

  end subroutine jac_CHO

  subroutine jac_HO(N,X,Y,DFY,LDFY,RPAR,IPAR)
    implicit none
    integer, intent(in) :: N, LDFY, ipar
    real(dp), intent(in) :: X, RPAR
    real(dp), dimension(N), intent(in) :: Y
    real(dp), dimension(LDFY, N),intent(out) :: DFY

    dfy(1, 1) = -re_f(1)*y(2) - re_r(2)*y(4)
    dfy(1, 2) = -re_f(1)*y(1) + re_f(2)*y(5)
    dfy(1, 3) = re_r(1)*y(4)
    dfy(1, 4) = re_r(1)*y(3) - re_r(2)*y(1)
    dfy(1, 5) = re_f(2)*y(2)
    dfy(2, 1) = -re_f(1)*y(2) + re_r(2)*y(4)
    dfy(2, 2) = -nd_atm*re_r(3) - re_f(1)*y(1) - re_f(2)*y(5)
    dfy(2, 3) = re_r(1)*y(4)
    dfy(2, 4) = 2.0_dp*nd_atm*re_f(3)*y(4) + re_r(1)*y(3) + re_r(2)*y(1)
    dfy(2, 5) = -re_f(2)*y(2)
    dfy(3, 1) = re_f(1)*y(2)
    dfy(3, 2) = re_f(1)*y(1)
    dfy(3, 3) = -re_r(1)*y(4)
    dfy(3, 4) = -re_r(1)*y(3)
    dfy(3, 5) = 0.0_dp
    dfy(4, 1) = re_f(1)*y(2) - re_r(2)*y(4)
    dfy(4, 2) = 2.0_dp*nd_atm*re_r(3) + re_f(1)*y(1) + re_f(2)*y(5)
    dfy(4, 3) = -re_r(1)*y(4)
    dfy(4, 4) = -4.0_dp*nd_atm*re_f(3)*y(4) - re_r(1)*y(3) - re_r(2)*y(1)
    dfy(4, 5) = re_f(2)*y(2)
    dfy(5, 1) = re_r(2)*y(4)
    dfy(5, 2) = -re_f(2)*y(5)
    dfy(5, 3) = 0.0_dp
    dfy(5, 4) = re_r(2)*y(1)
    dfy(5, 5) = -re_f(2)*y(2)

  end subroutine jac_HO

  subroutine mas_dummy(N,AM,LMAS,RPAR,IPAR)
    integer :: N, LMAS, IPAR
    double precision :: AM(LMAS,N), RPAR
  end subroutine mas_dummy

  subroutine solout(NR,XOLD,X,Y,CONT,LRC,N,RPAR,IPAR,IRTRN)
    integer :: NR, LRC, N, IPAR, IRTRN
    double precision :: XOLD, X, Y(N), CONT(LRC), RPAR
  end subroutine solout

end module mini_ch_i_radau5
