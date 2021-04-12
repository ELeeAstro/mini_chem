module mini_ch_i_rodas
  use mini_ch_precision
  use mini_ch_class
  use mini_ch_chem
  implicit none

  public ::  mini_ch_rodas, RHS_update, jac_dummy, mas_dummy, solout, dfx_dummy, jac_HO, jac_CHO

contains

  subroutine mini_ch_rodas(T, P, t_end, VMR)
    implicit none

    real(dp), intent(in) :: T, P, t_end
    real(dp), dimension(n_sp), intent(inout) :: VMR

    integer :: ncall
    real(dp) :: P_cgs

    ! Time controls
    real(dp) :: t_begin, t_now, dt_init

    ! RODAS variables
    real(dp), dimension(n_sp) :: y
    real(dp), allocatable, dimension(:) :: rwork
    integer, allocatable, dimension(:) :: iwork
    real(dp) :: rtol, atol
    real(dp) :: rpar
    integer :: itol, ijac, mljac, mujac, imas, mlmas, mumas, iout, lrwork, liwork, ipar, idid
    integer :: ifcn, idfx

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
    ! ***  parameters for the RODAS-solver  ***
    ! -----------------------------------------

    ifcn = 0
    idfx = 0
    rtol = 1.0e-4_dp
    atol = 1.0e-30_dp
    itol = 0
    ijac = 1
    mljac = n_sp
    mujac = 0
    imas = 0
    mlmas = n_sp
    mumas = 0
    iout = 0
    idid = 0

    ! Real work array
    lrwork = 2*n_sp*n_sp+14*n_sp+20
    allocate(rwork(lrwork))
    rwork(:) = 0.0_dp

    rwork(1) = 1.0e-16_dp ! Rounding unit - default 1e-16
    rwork(2) = t_end ! Max step size
    rwork(3) = 0.2_dp ! Parameter for step size selection - default 0.2
    rwork(4) = 6.0_dp ! Parameter for step size selection - deftaul 6.0
    rwork(5) = 0.9_dp ! Safety factor - default 0.9

    ! Integer work array
    liwork = n_sp + 20
    allocate(iwork(liwork))
    iwork(:) = 0

    iwork(1) = 0 ! Default step numbers - 0 = 100000
    iwork(2) = 1 ! Coefficent choice - default 1 (1,2,3)
    iwork(3) = 1 ! Step size strategy - default 1 (1,2)

    rpar = 0.0_dp
    ipar = 0


    t_begin = 0.0_dp
    t_now = t_begin
    dt_init = 1.0e-99_dp

    ncall = 1

    do while((t_now < t_end))

      ! Solution vector to seulex - y vector is number density of all species
      y(:) = g_sp(:)%nd

      call RODAS(n_sp,RHS_update,ifcn,t_now,y,t_end,dt_init, &
       &                  rtol,atol,itol, &
       &                  jac_HO,ijac,mljac,mujac, &
       &                  dfx_dummy,idfx, &
       &                  mas_dummy,imas,mlmas,mumas, &
       &                  solout,iout, &
       &                  rwork,lrwork,iwork,liwork,rpar,ipar,idid)

      g_sp(:)%nd = y(:)


      ncall = ncall + 1

    end do

    ! End result is update to VMR for all species
    VMR(:) = g_sp(:)%nd/sum(g_sp(:)%nd)

    deallocate(rwork, iwork)

  end subroutine mini_ch_rodas

  subroutine RHS_update(NEQ, time, y, f, rpar, ipar)
    implicit none

    integer, intent(in) ::  NEQ
    real(dp), intent(inout) :: time
    real(dp), dimension(NEQ), intent(inout) :: y
    real(dp), dimension(NEQ), intent(inout) :: f
    real(dp), intent(inout) :: rpar
    integer, intent(inout) :: ipar

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

  subroutine jac_dummy(N,X,Y,DFY,LDFY,RPAR,IPAR)
    integer :: N,LDFY,IPAR
    double precision :: X,Y(N),DFY(LDFY,N),RPAR
  end subroutine jac_dummy

  subroutine jac_CHO(N,X,Y,DFY,LDFY,RPAR,IPAR)
    implicit none
    integer, intent(in) :: N, LDFY, ipar
    real(dp), intent(in) :: X, RPAR
    real(dp), dimension(N), intent(in) :: Y
    real(dp), dimension(LDFY, N),intent(out) :: DFY

    nd_atm = sum(y(:))

    ! Update current number density of all species from y vector
    !g_sp(:)%nd = y(:)

    ! Find the forward, backward and net reaction rates
    !call reaction_rates(T, P_cgs)

    dfy(1, 1) = -re(1)%f*y(2) - re(2)%f*y(5) - re(3)%r*y(4) - re(4)%f*y(8)
    dfy(1, 2) = -re(1)%f*y(1) + re(3)%f*y(7)
    dfy(1, 3) = re(1)%r*y(4)
    dfy(1, 4) = re(1)%r*y(3) + re(2)%r*y(6) - re(3)%r*y(1) + re(4)%r*y(5)
    dfy(1, 5) = -re(2)%f*y(1) + re(4)%r*y(4)
    dfy(1, 6) = re(2)%r*y(4)
    dfy(1, 7) = re(3)%f*y(2)
    dfy(1, 8) = -re(4)%f*y(1)
    dfy(1, 9) = 0.0_dp
    dfy(1, 10) = 0.0_dp
    dfy(1, 11) = 0.0_dp
    dfy(2, 1) = -re(1)%f*y(2) + re(3)%r*y(4)
    dfy(2, 2) = -nd_atm*re(5)%r - 9.0_dp*re(6)%r*y(2)**2*y(5) - 9.0_dp*re(7)%r*y(2)**2*y(10) - &
      & re(8)%r*y(11) - re(9)%f*y(5)*y(8) - re(1)%f*y(1) - re(3)%f*y(7)
    dfy(2, 3) = 3.0_dp*re(6)%f*y(9) + re(1)%r*y(4)
    dfy(2, 4) = 2.0_dp*nd_atm*re(5)%f*y(4) + re(1)%r*y(3) + re(3)%r*y(1)
    dfy(2, 5) = -3.0_dp*re(6)%r*y(2)**3 - re(9)%f*y(2)*y(8)
    dfy(2, 6) = 0.0_dp
    dfy(2, 7) = re(9)%r*y(10) - re(3)%f*y(2)
    dfy(2, 8) = -re(9)%f*y(2)*y(5)
    dfy(2, 9) = 3.0_dp*re(6)%f*y(3) + 6.0_dp*re(7)%f*y(9) + 2.0_dp*re(8)%f*y(9)
    dfy(2, 10) = -3.0_dp*re(7)%r*y(2)**3 + re(9)%r*y(7)
    dfy(2, 11) = -re(8)%r*y(2)
    dfy(3, 1) = re(1)%f*y(2)
    dfy(3, 2) = 3.0_dp*re(6)%r*y(2)**2*y(5) + re(1)%f*y(1)
    dfy(3, 3) = -re(6)%f*y(9) - re(1)%r*y(4)
    dfy(3, 4) = -re(1)%r*y(3)
    dfy(3, 5) = re(6)%r*y(2)**3
    dfy(3, 6) = 0.0_dp
    dfy(3, 7) = 0.0_dp
    dfy(3, 8) = 0.0_dp
    dfy(3, 9) = -re(6)%f*y(3)
    dfy(3, 10) = 0.0_dp
    dfy(3, 11) = 0.0_dp
    dfy(4, 1) = re(1)%f*y(2) + re(2)%f*y(5) - re(3)%r*y(4) + re(4)%f*y(8)
    dfy(4, 2) = 2.0_dp*nd_atm*re(5)%r + re(1)%f*y(1) + re(3)%f*y(7)
    dfy(4, 3) = -re(1)%r*y(4)
    dfy(4, 4) = -4.0_dp*nd_atm*re(5)%f*y(4) - re(1)%r*y(3) - re(2)%r*y(6) - re(3)%r*y(1) - re(4)%r*y(5)
    dfy(4, 5) = re(2)%f*y(1) - re(4)%r*y(4)
    dfy(4, 6) = -re(2)%r*y(4)
    dfy(4, 7) = re(3)%f*y(2)
    dfy(4, 8) = re(4)%f*y(1)
    dfy(4, 9) = 0.0_dp
    dfy(4, 10) = 0.0_dp
    dfy(4, 11) = 0.0_dp
    dfy(5, 1) = -re(2)%f*y(5) + re(4)%f*y(8)
    dfy(5, 2) = -3.0_dp*re(6)%r*y(2)**2*y(5) - re(9)%f*y(5)*y(8)
    dfy(5, 3) = re(6)%f*y(9)
    dfy(5, 4) = re(2)%r*y(6) - re(4)%r*y(5)
    dfy(5, 5) = -re(6)%r*y(2)**3 - re(9)%f*y(2)*y(8) - re(2)%f*y(1) - re(4)%r*y(4)
    dfy(5, 6) = re(2)%r*y(4)
    dfy(5, 7) = re(9)%r*y(10)
    dfy(5, 8) = -re(9)%f*y(2)*y(5) + re(4)%f*y(1)
    dfy(5, 9) = re(6)%f*y(3)
    dfy(5, 10) = re(9)%r*y(7)
    dfy(5, 11) = 0.0_dp
    dfy(6, 1) = re(2)%f*y(5)
    dfy(6, 2) = 0.0_dp
    dfy(6, 3) = 0.0_dp
    dfy(6, 4) = -re(2)%r*y(6)
    dfy(6, 5) = re(2)%f*y(1)
    dfy(6, 6) = -re(2)%r*y(4)
    dfy(6, 7) = 0.0_dp
    dfy(6, 8) = 0.0_dp
    dfy(6, 9) = 0.0_dp
    dfy(6, 10) = 0.0_dp
    dfy(6, 11) = 0.0_dp
    dfy(7, 1) = re(3)%r*y(4)
    dfy(7, 2) = re(9)%f*y(5)*y(8) - re(3)%f*y(7)
    dfy(7, 3) = 0.0_dp
    dfy(7, 4) = re(3)%r*y(1)
    dfy(7, 5) = re(9)%f*y(2)*y(8)
    dfy(7, 6) = 0.0_dp
    dfy(7, 7) = -re(9)%r*y(10) - re(3)%f*y(2)
    dfy(7, 8) = re(9)%f*y(2)*y(5)
    dfy(7, 9) = 0.0_dp
    dfy(7, 10) = -re(9)%r*y(7)
    dfy(7, 11) = 0.0_dp
    dfy(8, 1) = -re(4)%f*y(8)
    dfy(8, 2) = -re(9)%f*y(5)*y(8)
    dfy(8, 3) = 0.0_dp
    dfy(8, 4) = re(4)%r*y(5)
    dfy(8, 5) = -re(9)%f*y(2)*y(8) + re(4)%r*y(4)
    dfy(8, 6) = 0.0_dp
    dfy(8, 7) = re(9)%r*y(10)
    dfy(8, 8) = -re(9)%f*y(2)*y(5) - re(4)%f*y(1)
    dfy(8, 9) = 0.0_dp
    dfy(8, 10) = re(9)%r*y(7)
    dfy(8, 11) = 0.0_dp
    dfy(9, 1) = 0.0_dp
    dfy(9, 2) = 3.0_dp*re(6)%r*y(2)**2*y(5) + 6.0_dp*re(7)%r*y(2)**2*y(10) + 2.0_dp*re(8)%r*y(11)
    dfy(9, 3) = -re(6)%f*y(9)
    dfy(9, 4) = 0.0_dp
    dfy(9, 5) = re(6)%r*y(2)**3
    dfy(9, 6) = 0.0_dp
    dfy(9, 7) = 0.0_dp
    dfy(9, 8) = 0.0_dp
    dfy(9, 9) = -re(6)%f*y(3) - 4.0_dp*re(7)%f*y(9) - 4.0_dp*re(8)%f*y(9)
    dfy(9, 10) = 2.0_dp*re(7)%r*y(2)**3
    dfy(9, 11) = 2.0_dp*re(8)%r*y(2)
    dfy(10, 1) = 0.0_dp
    dfy(10, 2) = -3.0_dp*re(7)%r*y(2)**2*y(10) + re(9)%f*y(5)*y(8)
    dfy(10, 3) = 0.0_dp
    dfy(10, 4) = 0.0_dp
    dfy(10, 5) = re(9)%f*y(2)*y(8)
    dfy(10, 6) = 0.0_dp
    dfy(10, 7) = -re(9)%r*y(10)
    dfy(10, 8) = re(9)%f*y(2)*y(5)
    dfy(10, 9) = 2.0_dp*re(7)%f*y(9)
    dfy(10, 10) = -re(7)%r*y(2)**3 - re(9)%r*y(7)
    dfy(10, 11) = 0.0_dp
    dfy(11, 1) = 0.0_dp
    dfy(11, 2) = -re(8)%r*y(11)
    dfy(11, 3) = 0.0_dp
    dfy(11, 4) = 0.0_dp
    dfy(11, 5) = 0.0_dp
    dfy(11, 6) = 0.0_dp
    dfy(11, 7) = 0.0_dp
    dfy(11, 8) = 0.0_dp
    dfy(11, 9) = 2.0_dp*re(8)%f*y(9)
    dfy(11, 10) = 0.0_dp
    dfy(11, 11) = -re(8)%r*y(2)

  end subroutine jac_CHO

  subroutine jac_HO(N,X,Y,DFY,LDFY,RPAR,IPAR)
    implicit none
    integer, intent(in) :: N, LDFY, ipar
    real(dp), intent(in) :: X, RPAR
    real(dp), dimension(N), intent(in) :: Y
    real(dp), dimension(LDFY, N),intent(out) :: DFY

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

  subroutine mas_dummy(N,AM,LMAS,RPAR,IPAR)
    integer :: N, LMAS, IPAR
    double precision :: AM(LMAS,N), RPAR
  end subroutine mas_dummy

  subroutine dfx_dummy(N,X,Y,FX,RPAR,IPAR)
    integer :: N,IPAR
    double precision :: X,Y(N),FX(N),RPAR
  end subroutine dfx_dummy

  subroutine solout(NR,XOLD,X,Y,CONT,LRC,N,RPAR,IPAR,IRTRN)
    integer :: NR, LRC, N, IPAR, IRTRN
    double precision :: XOLD, X, Y(N), CONT(LRC), RPAR
  end subroutine solout

end module mini_ch_i_rodas
