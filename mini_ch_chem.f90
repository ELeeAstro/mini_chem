module mini_ch_chem
  use mini_ch_precision
  use mini_ch_class
  implicit none


public :: reaction_rates, reverse_reactions, check_con

contains

  subroutine reaction_rates(T, P)
    implicit none

    real(dp), intent(in) :: T, P

    integer :: i, iT, iT1, iP, iP1
    real(dp) :: k0, kinf
    real(dp) :: Tr
    real(dp) :: Tl, Tu, Pl, Pu, Tw, Pw
    real(dp) :: k11, k12, k21, k22, kf

    !! Calculate the forward and backward, and then net reaction rates for each reaction
    do i = 1, n_reac

      Tr = T

      ! if (T < re(i)%Tmin) then
      !   Tr = re(i)%Tmin
      ! else if (T > re(i)%Tmax) then
      !   Tr = re(i)%Tmax
      ! end if

      if (re(i)%re_t == 2) then
        ! Two body reaction
        re(i)%f = re(i)%A * Tr**re(i)%B * exp(-re(i)%C/Tr)

      else if (re(i)%re_t == 3) then
        ! Three body reaction
        k0 = re(i)%A0 * Tr**re(i)%B0 * exp(-re(i)%C0/Tr)
        kinf = re(i)%Ainf * Tr**re(i)%Binf * exp(-re(i)%Cinf/Tr)

        re(i)%f = (k0 * nd_atm) / (1.0_dp + (k0 * nd_atm/kinf))

        !print*, i, k0, kinf
      else if (re(i)%re_t == 4) then
        ! Interpolate from table to find net reaction forward rate

        ! Find temperature index
        Tw = Tr
        call locate(re(i)%T(:), re(i)%nT , Tw, iT)
        if (iT == 0) then
          iT = 1
          Tw = minval(re(i)%T(:))
        else if (iT == re(i)%nT) then
          iT = re(i)%nT-1
          Tw = maxval(re(i)%T(:))
        end if
        iT1 = iT + 1

        Tl = re(i)%T(iT)
        Tu = re(i)%T(iT1)

        ! Find pressure index
        Pw = P
        call locate(re(i)%P(:), re(i)%nP , Pw, iP)
        if (iP == 0) then
          iP = 1
          Pw = minval(re(i)%P(:))
        else if (iP == re(i)%nP) then
          iP = re(i)%nP-1
          Pw = maxval(re(i)%P(:))
        end if
        iP1 = iP + 1

        Pl = re(i)%P(iP)
        Pu = re(i)%P(iP1)

        ! Bi-linearly interpolate from table to find kf
        k11 = re(i)%kf(iT,iP)
        k12 = re(i)%kf(iT,iP1)
        k21 = re(i)%kf(iT1,iP)
        k22 = re(i)%kf(iT1,iP1)

       !call bilinear_interp(Tw, Pw, Tl, Tu, Pl, Pu, k11, k21, k12, k22, kf)
       call bilinear_log_interp(Tw, Pw, Tl, Tu, Pl, Pu, k11, k21, k12, k22, kf)


       ! print*, i
       ! print*, Tl, Tu
       ! print*, Pl, Pu
       ! print*, k11, k12, k21, k22, kf

       re(i)%f = kf

      end if

      re(i)%r = re(i)%f/re(i)%Keq
      re(i)%net = re(i)%f - re(i)%r

      !print*, i, T, re(i)%f, re(i)%r, re(i)%net

    end do


  end subroutine reaction_rates

  subroutine reverse_reactions(T, P)
    implicit none

    real(dp), intent(in) :: T, P

    integer :: i, j
    real(dp) :: Tn2, Tn1, lnT, T2, T3, T4, Tr

    Tr = T

    ! if (T < 200.0_dp) then
    !   Tr = 200.0_dp
    ! else if (T > 6000.0_dp) then
    !   Tr = 6000.0_dp
    ! end if

    Tn2 = 1.0_dp/Tr**2
    Tn1 = 1.0_dp/Tr
    lnT = log(Tr)
    T2 = Tr**2
    T3 = Tr**3
    T4 = Tr**4

    !! First calculate H0 and s0 from the polynomials
    do i = 1, n_sp
      if (Tr <= 1000.0_dp) then
        g_sp(i)%H0 = -g_sp(i)%a_l(1)*Tn2 + g_sp(i)%a_l(2)*lnT*Tn1 + g_sp(i)%a_l(3) &
          & + g_sp(i)%a_l(4)*Tr/2.0_dp + g_sp(i)%a_l(5)*T2/3.0_dp + g_sp(i)%a_l(6)*T3/4.0_dp &
          & + g_sp(i)%a_l(7)*T4/5.0_dp + g_sp(i)%a_l(8)*Tn1
        g_sp(i)%s0 = -g_sp(i)%a_l(1)*Tn2/2.0_dp - g_sp(i)%a_l(2)*Tn1 + g_sp(i)%a_l(3)*lnT &
          & + g_sp(i)%a_l(4)*Tr + g_sp(i)%a_l(5)*T2/2.0_dp + g_sp(i)%a_l(6)*T3/3.0_dp &
          & + g_sp(i)%a_l(7)*T4/4.0_dp + g_sp(i)%a_l(9)
      else
        g_sp(i)%H0 = -g_sp(i)%a_h(1)*Tn2 + g_sp(i)%a_h(2)*lnT*Tn1 + g_sp(i)%a_h(3) &
          & + g_sp(i)%a_h(4)*Tr/2.0_dp + g_sp(i)%a_h(5)*T2/3.0_dp + g_sp(i)%a_h(6)*T3/4.0_dp &
          & + g_sp(i)%a_h(7)*T4/5.0_dp + g_sp(i)%a_h(8)*Tn1
        g_sp(i)%s0 = -g_sp(i)%a_h(1)*Tn2/2.0_dp - g_sp(i)%a_h(2)*Tn1 + g_sp(i)%a_h(3)*lnT &
          & + g_sp(i)%a_h(4)*Tr + g_sp(i)%a_h(5)*T2/2.0_dp + g_sp(i)%a_h(6)*T3/3.0_dp &
          & + g_sp(i)%a_h(7)*T4/4.0_dp + g_sp(i)%a_h(9)
      end if

      g_sp(i)%H0 = g_sp(i)%H0 * R*Tr
      g_sp(i)%s0 = g_sp(i)%s0 * R

    end do

    !! Second calculate the reverse reaction coefficent
    re(:)%dH = 0.0_dp
    re(:)%ds = 0.0_dp
    do i = 1, n_reac
      do j = 1, re(i)%n_pr
        re(i)%dH = re(i)%dH + g_sp(re(i)%gi_pr(j))%H0
        re(i)%ds = re(i)%ds + g_sp(re(i)%gi_pr(j))%s0
      end do
      do j = 1, re(i)%n_re
        re(i)%dH = re(i)%dH - g_sp(re(i)%gi_re(j))%H0
        re(i)%ds = re(i)%ds - g_sp(re(i)%gi_re(j))%s0
      end do

      re(i)%Keq = exp(-(re(i)%dH - Tr*re(i)%ds)/(R*Tr)) * ((kb * Tr)/P0)**(-re(i)%dmu)

      !print*, i, re(i)%dH, re(i)%ds, re(i)%Keq

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
