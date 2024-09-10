module mini_ch_chem
  use mini_ch_precision
  use mini_ch_class
  implicit none


  private
  public :: reaction_rates, reverse_reactions, check_con

contains

  subroutine reaction_rates(T, P, nd_atm)
    implicit none

    real(dp), intent(in) :: T, P, nd_atm

    integer :: i, iT1, iT2, iT3, iP1, iP2, iP3
    real(dp) :: k0, kinf
    real(dp) :: lT, lP
    real(dp) :: kf
    real(dp), dimension(3) :: lPa, lTa, lkf, lkfa

    !! Calculate the forward and backward, and then net reaction rates for each reaction
    do i = 1, n_reac

      if (re(i)%re_t == 2) then
        ! Two body reaction
        kf = re(i)%A * T**re(i)%B * exp(-re(i)%C/T)

      else if (re(i)%re_t == 3) then
        ! Three body reaction
        k0 = re(i)%A0 * T**re(i)%B0 * exp(-re(i)%C0/T)
        kinf = re(i)%Ainf * T**re(i)%Binf * exp(-re(i)%Cinf/T)

        kf = k0 / (1.0_dp + ((k0 * nd_atm)/kinf))

        !print*, i, k0, kinf
      else if (re(i)%re_t == 4) then

        lT = log10(T)
        lP = log10(P)

        ! Interpolate using Bezier interpolation to find net reaction forward rate
        ! from the tables

        ! Find pressure index triplet
        call locate(re(i)%P(:), re(i)%nP, P, iP2)
        iP1 = iP2 - 1
        iP3 = iP2 + 1

        if (iP1 <= 0) then
          iP1 = 1
          iP2 = 2
          iP3 = 3
        else if (iP3 > re(i)%nP) then
          iP1 = re(i)%nP - 2
          iP2 = re(i)%nP - 1
          iP3 = re(i)%nP
        end if

        lPa(1) = re(i)%lP(iP1)
        lPa(2) = re(i)%lP(iP2)
        lPa(3) = re(i)%lP(iP3)

        ! Check if input temperature is within table range
        if (T <= re(i)%T(1)) then

          ! Perform Bezier interpolation at minimum table temperature
          lkf(1) = re(i)%lkf(1,iP1)
          lkf(2) = re(i)%lkf(1,iP2)
          lkf(3) = re(i)%lkf(1,iP3)
          call Bezier_interp(lPa(:), lkf(:), 3, lP, kf)
          kf = 10.0_dp**kf

        else if (T >= re(i)%T(re(i)%nT)) then

          ! Perform Bezier interpolation at maximum table temperature
          lkf(1) = re(i)%lkf(re(i)%nT,iP1)
          lkf(2) = re(i)%lkf(re(i)%nT,iP2)
          lkf(3) = re(i)%lkf(re(i)%nT,iP3)
          call Bezier_interp(lPa(:), lkf(:), 3, lP, kf)
          kf = 10.0_dp**kf

        else

          ! Pressure and temperature is within table grid
          ! Perform 2D Bezier interpolation by performing interpolation 4 times

          ! Find temperature index triplet
          call locate(re(i)%T(:), re(i)%nT, T, iT2)
          iT1 = iT2 - 1
          iT3 = iT2 + 1

          if (iT1 < 1) then
            iT1 = 1
            iT2 = 2
            iT3 = 3
          else if (iT3 > re(i)%nT) then
            iT1 = re(i)%nT - 2
            iT2 = re(i)%nT - 1
            iT3 = re(i)%nT
          end if

          lkf(1) = re(i)%lkf(iT1,iP1)
          lkf(2) = re(i)%lkf(iT1,iP2)
          lkf(3) = re(i)%lkf(iT1,iP3)
          call Bezier_interp(lPa(:), lkf(:), 3, lP, lkfa(1)) ! Result at T1, P_in
          lkf(1) = re(i)%lkf(iT2,iP1)
          lkf(2) = re(i)%lkf(iT2,iP2)
          lkf(3) = re(i)%lkf(iT2,iP3)
          call Bezier_interp(lPa(:), lkf(:), 3, lP, lkfa(2)) ! Result at T2, P_in
          lkf(1) = re(i)%lkf(iT3,iP1)
          lkf(2) = re(i)%lkf(iT3,iP2)
          lkf(3) = re(i)%lkf(iT3,iP3)
          call Bezier_interp(lPa(:), lkf(:), 3, lP, lkfa(3)) ! Result at T3, P_in
          lTa(1) = re(i)%lT(iT1)
          lTa(2) = re(i)%lT(iT2)
          lTa(3) = re(i)%lT(iT3)
          call Bezier_interp(lTa(:), lkfa(:), 3, lT, kf) ! Result at T, P
          kf = 10.0_dp**kf

          ! print*, i, kf

        end if
      
      end if

      !! Limit for very small rates
      re_f(i) = kf

      re_r(i) = re_f(i)/Keq(i) * ((kb * T)/P0)**(re(i)%dmu)

    end do

  end subroutine reaction_rates

  subroutine reverse_reactions(T)
    implicit none

    real(dp), intent(in) :: T

    integer :: i, j
    real(dp) :: Tn2, Tn1, lnT, T2, T3, T4, Tr
    real(dp), dimension(n_sp) :: H0, s0
    real(dp), dimension(n_reac) :: dH, ds

    Tr = T

    Tn2 = 1.0_dp/Tr**2
    Tn1 = 1.0_dp/Tr
    lnT = log(Tr)
    T2 = Tr**2
    T3 = Tr**3
    T4 = Tr**4

    !! First calculate H0 and s0 from the polynomials
    do i = 1, n_sp
      if (Tr <= 1000.0_dp) then
        H0(i) = -g_sp(i)%a_l(1)*Tn2 + g_sp(i)%a_l(2)*lnT*Tn1 + g_sp(i)%a_l(3) &
          & + g_sp(i)%a_l(4)*Tr/2.0_dp + g_sp(i)%a_l(5)*T2/3.0_dp + g_sp(i)%a_l(6)*T3/4.0_dp &
          & + g_sp(i)%a_l(7)*T4/5.0_dp + g_sp(i)%a_l(8)*Tn1
        s0(i) = -g_sp(i)%a_l(1)*Tn2/2.0_dp - g_sp(i)%a_l(2)*Tn1 + g_sp(i)%a_l(3)*lnT &
          & + g_sp(i)%a_l(4)*Tr + g_sp(i)%a_l(5)*T2/2.0_dp + g_sp(i)%a_l(6)*T3/3.0_dp &
          & + g_sp(i)%a_l(7)*T4/4.0_dp + g_sp(i)%a_l(9)
      else
        H0(i) = -g_sp(i)%a_h(1)*Tn2 + g_sp(i)%a_h(2)*lnT*Tn1 + g_sp(i)%a_h(3) &
          & + g_sp(i)%a_h(4)*Tr/2.0_dp + g_sp(i)%a_h(5)*T2/3.0_dp + g_sp(i)%a_h(6)*T3/4.0_dp &
          & + g_sp(i)%a_h(7)*T4/5.0_dp + g_sp(i)%a_h(8)*Tn1
        s0(i) = -g_sp(i)%a_h(1)*Tn2/2.0_dp - g_sp(i)%a_h(2)*Tn1 + g_sp(i)%a_h(3)*lnT &
          & + g_sp(i)%a_h(4)*Tr + g_sp(i)%a_h(5)*T2/2.0_dp + g_sp(i)%a_h(6)*T3/3.0_dp &
          & + g_sp(i)%a_h(7)*T4/4.0_dp + g_sp(i)%a_h(9)
      end if

      H0(i) = H0(i) * R*Tr
      s0(i) = s0(i) * R

    end do

    !! Second calculate the reverse reaction coefficent
    dH(:) = 0.0_dp
    ds(:) = 0.0_dp
    do i = 1, n_reac
      do j = 1, re(i)%n_re
        dH(i) = dH(i) - H0(re(i)%gi_re(j))
        ds(i) = ds(i) - s0(re(i)%gi_re(j))
      end do
      do j = 1, re(i)%n_pr
        dH(i) = dH(i) + H0(re(i)%gi_pr(j))
        ds(i) = ds(i) + s0(re(i)%gi_pr(j))
      end do

      Keq(i) = exp(-(dH(i) - Tr*ds(i))/(R*Tr))

    end do

  end subroutine reverse_reactions

  subroutine check_con(n_sp, n_kp, n_k, t_now, t_old, con)
    implicit none

    integer, intent(in) :: n_sp
    real(dp), dimension(n_sp), intent(in) :: n_kp, n_k
    real(dp), intent(in) :: t_now, t_old

    logical, intent(out) :: con

    integer :: n
    real(dp), dimension(n_sp) :: dn, eps
    real(dp) :: dt

    dt = t_now - t_old
    if (dt <= 0.0_dp) then
      con = .False.
      return
    end if

    do n = 1, n_sp
      dn(n) = abs(n_k(n) - n_kp(n))/n_k(n)
      eps(n) = dn(n)/dt
      if ((dn(n) > del_con) .and. (eps(n) > eps_con)) then
        con = .False.
        return
      end if
    end do

    con = .True.
    !print*, 'converged'
    !print*, t_now, t_old, dt, con, eps(:), dn(:)

  end subroutine check_con

end module mini_ch_chem
