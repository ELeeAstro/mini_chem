module mini_ch_i_sdirk
  use mini_ch_precision
  use mini_ch_class
  use mini_ch_chem
  use SDIRK_f90_Integrator, only : INTEGRATE
  implicit none

  public ::  mini_ch_sdirk, RHS_update, jac_dummy, &
  &  jac_HO, jac_CHO, jac_NCHO

contains

  subroutine mini_ch_sdirk(T, P, t_end, VMR, network)
    implicit none

    real(dp), intent(in) :: T, P, t_end
    real(dp), dimension(n_sp), intent(inout) :: VMR
    character(len=200), intent(in) :: network

    integer :: ncall
    real(dp) :: P_cgs

    ! Time controls
    real(dp) :: t_begin, t_now, dt_init, t_old
    logical :: con = .False.

    ! Rosenbrock variables
    integer :: ierr, nnzero
    real(dp), dimension(n_sp) :: rtol, atol
    real(dp), dimension(n_sp) :: y
    integer, dimension(20) :: icntrl, istatus
    real(dp), dimension(20) :: rcntrl, rstatus

    !! Find the number density of the atmosphere
    P_cgs = P * 10.0_dp   ! Convert pascal to dyne cm-2
    nd_atm = P_cgs/(kb*T)  ! Find number density (cm-3] of atmosphere

    !! Find initial number density of all species from VMR
    g_sp(:)%nd = nd_atm * VMR(:)

    ! First find the reverse reaction coefficents (Keq)
    call reverse_reactions(T, P_cgs)
    ! Find the forward, backward and net reaction rates
    call reaction_rates(T, P_cgs)


    ! -----------------------------------------
    ! ***  parameters for the Rosenbrock-solver  ***
    ! -----------------------------------------

    t_begin = 0.0_dp
    t_now = t_begin
    dt_init = 1.0e-99_dp

    rtol(:) = 1.0e-3_dp
    atol(:) = 1e-99_dp

    nnzero = n_sp

    icntrl(1:20) = 0

    icntrl(1) = 0
    icntrl(2) = 0
    icntrl(3) = 0
    icntrl(4) = 0
    icntrl(5) = 0
    icntrl(6) = 0

    rcntrl(1:20) = 0.0_dp

    istatus(1:20) = 0
    rstatus(1:20) = 0.0_dp

    ncall = 1


    do while((t_now < t_end))

      ! Solution vector to seulex - y vector is number density of all species
      y(:) = g_sp(:)%nd

      t_old = t_now

      select case(network)
      case('HO')
        call INTEGRATE(t_now, t_end, n_sp, nnzero, y, rtol, atol, RHS_update, jac_HO,&
          icntrl, rcntrl, istatus, rstatus, ierr )
      case('CHO')
        call INTEGRATE(t_now, t_end, n_sp, nnzero, y, rtol, atol, RHS_update, jac_CHO,&
          icntrl, rcntrl, istatus, rstatus, ierr )
      case('NCHO')
        call INTEGRATE(t_now, t_end, n_sp, nnzero, y, rtol, atol, RHS_update, jac_NCHO,&
          icntrl, rcntrl, istatus, rstatus, ierr )
      case default
        print*, 'Invalid network provided: ', trim(network)
        stop
      end select

      t_now = rstatus(1)
      rcntrl(3) = rstatus(2)

      if (t_now >= t_end*f_con) then
        call check_con(n_sp,g_sp(:)%nd,y(:),t_now,t_old,con)
        if (con .eqv. .True.) then
          g_sp(:)%nd = y(:)
          exit
        end if
      end if

      g_sp(:)%nd = y(:)


      ncall = ncall + 1

    end do

    VMR(:) = g_sp(:)%nd/nd_atm

  end subroutine mini_ch_sdirk

  subroutine RHS_update(n, time, y, f)
    implicit none

    integer, intent(in) :: n
    real(dp), intent(in) :: time
    real(dp), dimension(n), intent(in) :: y
    real(dp), dimension(n), intent(inout) :: f

    integer :: i, j, k
    real(dp) :: msum, msum2

    nd_atm = sum(y(:))

    ! Update current number density of all species from y vector
    g_sp(:)%nd = y(:)

    ! Calculate the rate of change of number density for all species (cm-3/s]
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

      ! Find flux contribution for products
      f(re(i)%gi_pr(:)) = f(re(i)%gi_pr(:)) + msum * re(i)%f - msum2 * re(i)%r
      ! Find flux contribution for reactants
      f(re(i)%gi_re(:)) = f(re(i)%gi_re(:)) + msum2 * re(i)%r - msum * re(i)%f

    end do

  end subroutine RHS_update

  subroutine jac_dummy(n, time, Y, DFY)
    integer, intent(in) :: n
    real(dp), intent(in) :: time
    real(dp), dimension(n), intent(in) :: y
    real(dp), dimension(n,n), intent(inout) :: DFY
  end subroutine jac_dummy

  subroutine jac_NCHO(n, time, Y, DFY)
    implicit none

    integer, intent(in) :: n
    real(dp), intent(in) :: time
    real(dp), dimension(n), intent(in) :: Y
    real(dp), dimension(n,n),intent(inout) :: DFY

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

  subroutine jac_CHO(n, time, Y, DFY)
    implicit none

    integer, intent(in) :: n
    real(dp), intent(in) :: time
    real(dp), dimension(n), intent(in) :: Y
    real(dp), dimension(n,n),intent(inout) :: DFY

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

  subroutine jac_HO(n, time, Y, DFY)
    implicit none

    integer, intent(in) :: n
    real(dp), intent(in) :: time
    real(dp), dimension(n), intent(in) :: Y
    real(dp), dimension(n,n),intent(inout) :: DFY

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

end module mini_ch_i_sdirk
