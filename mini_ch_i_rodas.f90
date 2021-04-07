module mini_ch_i_rodas
  use mini_ch_precision
  use mini_ch_class
  use mini_ch_chem
  implicit none

  public ::  mini_ch_rodas, RHS_update, jac_dummy, mas_dummy, solout, dfx_dummy

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
    nd_atm = P_cgs/(kb*T)  ! Find number density [cm-3] of atmosphere

    !! Find initial number density of all species from VMR
    g_sp(:)%nd = nd_atm * VMR(:)

    ! First find the reverse reaction coefficents (Keq)
    call reverse_reactions(T)
    ! Find the forward, backward and net reaction rates
    call reaction_rates(T)


    ! -----------------------------------------
    ! ***  parameters for the RODAS-solver  ***
    ! -----------------------------------------

    ifcn = 0
    idfx = 0
    rtol = 1.0e-4_dp
    atol = 1.0e-30_dp
    itol = 0
    ijac = 0
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
    iwork(2) = 3 ! Coefficent choice - default 1 (1,2,3)
    iwork(3) = 1 ! Step size strategy - default 1 (1,2)

    rpar = 0.0_dp
    ipar = 0


    t_begin = 0.0_dp
    t_now = t_begin
    dt_init = 1.0e-12_dp

    ncall = 1

    do while((t_now < t_end))

      ! Solution vector to seulex - y vector is number density of all species
      y(:) = g_sp(:)%nd

      call RODAS(n_sp,RHS_update,ifcn,t_now,y,t_end,dt_init, &
       &                  rtol,atol,itol, &
       &                  jac_dummy,ijac,mljac,mujac, &
       &                  dfx_dummy,idfx, &
       &                  mas_dummy,imas,mlmas,mumas, &
       &                  solout,iout, &
       &                  rwork,lrwork,iwork,liwork,rpar,ipar,idid)

      g_sp(:)%nd = y(:)

      exit

      ncall = ncall + 1

    end do

    deallocate(rwork, iwork)

    ! End result is update to VMR for all species
    VMR(:) = g_sp(:)%nd/sum(g_sp(:)%nd)

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

    ! Update current number density of all species from y vector
    g_sp(:)%nd = y(:)

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
      if (re(i)%re_t == 3) then
        ! Mutliply by atmosphere nd if neutral body involved
        msum = msum * nd_atm
      end if
      ! Find number density multiple for products in reaction
      msum2 = g_sp(re(i)%gi_pr(1))%nd
      do k = 2, re(i)%n_pr
         msum2 = msum2 * g_sp(re(i)%gi_pr(k))%nd
      end do
      if (re(i)%re_t == 3) then
        ! Mutliply by atmosphere nd if neutral body involved
        msum2 = msum2 * nd_atm
      end if

      ! Find flux for products
      f(re(i)%gi_pr(:)) = f(re(i)%gi_pr(:)) + msum * re(i)%f - msum2 * re(i)%r
      f(re(i)%gi_re(:)) = f(re(i)%gi_re(:)) + msum2 * re(i)%r - msum * re(i)%f

    end do

  end subroutine RHS_update

  subroutine jac_dummy(N,X,Y,DFY,LDFY,RPAR,IPAR)
    integer :: N,LDFY,IPAR
    double precision :: X,Y(N),DFY(LDFY,N),RPAR
  end subroutine jac_dummy

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
