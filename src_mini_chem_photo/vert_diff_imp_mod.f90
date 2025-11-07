module vert_diff_imp_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  integer, parameter :: dp = REAL64

  ! -------- Numerics / constants --------
  real(dp), parameter :: qmin   = 1.0e-99_dp
  real(dp), parameter :: eps    = 1.0e-300_dp
  real(dp), parameter :: R_gas  = 8.31446261815324e7_dp   ! erg/mol/K

  ! -------- Manual theta only --------
  ! THETA = 0.5_dp -> Crank–Nicolson (A-stable, may ring if very stiff)
  ! THETA = 1.0_dp -> Backward Euler (L-stable)
  real(dp), parameter :: THETA  = 0.5_dp

  public  :: vert_diff_imp
  private :: harm_mean, tridiag_factor, tridiag_solve_fact

contains

  subroutine vert_diff_imp(nlay, nlev, t_end, mu, grav_in, Tl, pl_in, pe_in, Kzz, nq, q, q0)
    implicit none
    ! ---- Inputs ----
    integer, intent(in) :: nlay, nlev, nq
    real(dp), intent(in) :: t_end, grav_in
    real(dp), intent(in) :: Tl(nlay), pl_in(nlay), Kzz(nlay), mu(nlay)
    real(dp), intent(in) :: pe_in(nlev)
    real(dp), intent(in) :: q0(nq)                ! Dirichlet value at bottom face (per tracer)

    ! ---- InOut ----
    real(dp), intent(inout) :: q(nlay,nq)         ! cell-centred tracer(s)

    ! ---- Geometry/state (invariants across tracers) ----
    real(dp) :: pl(nlay), pe(nlev)
    real(dp) :: alte(nlev), altm(nlay)
    real(dp) :: dz(nlay), dzm(nlay-1)
    real(dp) :: rho(nlay), rho_e(nlev), K_e(nlev)
    real(dp) :: Df(nlay+1), wcell(nlay)           ! face conductance and cell weight 1/(rho*dz)
    real(dp) :: inv_dt

    ! ---- Tridiagonal system (invariants across tracers) ----
    real(dp) :: a(nlay), b(nlay), c(nlay)
    real(dp) :: bp(nlay), cp(nlay)                ! factored diagonal & upper ratios

    ! ---- Per-tracer work ----
    real(dp) :: rhs(nlay)

    ! ---- Scalars ----
    integer  :: k, n
    real(dp) :: grav, Te, mue

    ! =========================
    ! Units + hydrostatic grid
    ! =========================
    pl   = 10.0_dp * pl_in                ! Pa -> dyne/cm^2
    pe   = 10.0_dp * pe_in
    grav = 100.0_dp * grav_in             ! m/s^2 -> cm/s^2

    alte(nlev) = 0.0_dp
    do k = nlev-1, 1, -1
      alte(k) = alte(k+1) + (R_gas*Tl(k))/(mu(k)*grav) * log(pe(k+1)/pe(k))
    end do
    do k = 1, nlay
      dz(k)   = alte(k) - alte(k+1)
      altm(k) = 0.5_dp*(alte(k) + alte(k+1))
    end do
    do k = 1, nlay-1
      dzm(k) = altm(k) - altm(k+1)
    end do

    ! Densities (centres/faces)
    rho  = pl / ((R_gas/mu) * Tl)
    rho_e(1) = pe(1) / ((R_gas/mu(1)) * Tl(1))
    do k = 2, nlay
      Te  = 0.5_dp*(Tl(k-1) + Tl(k))
      mue = 0.5_dp*(mu(k-1) + mu(k))
      rho_e(k) = pe(k) / ((R_gas/mue) * Te)
    end do
    rho_e(nlev) = pe(nlev) / ((R_gas/mu(nlay)) * Tl(nlay))

    ! Eddy K on faces (harmonic mean is robust across jumps)
    K_e(1) = Kzz(1)
    do k = 2, nlay
      K_e(k) = harm_mean(Kzz(k-1), Kzz(k))
    end do
    K_e(nlev) = Kzz(nlay)

    ! Face conductances Df(face) = (rho*K)_face / Δz_mid
    Df(1) = 0.0_dp                                 ! top Neumann (zero flux)
    do k = 2, nlay
      Df(k) = rho_e(k) * K_e(k) / (dzm(k-1) + eps) ! interior faces
    end do
    Df(nlay+1) = rho_e(nlev) * K_e(nlev) / (dzm(nlay-1) + eps)   ! bottom face

    ! Cell weight wcell = 1/(rho*Δz)
    do k = 1, nlay
      wcell(k) = 1.0_dp / (rho(k)*dz(k) + eps)
    end do

    inv_dt = 1.0_dp / max(t_end, 1.0e-300_dp)

    ! =========================
    ! Assemble matrix once
    ! =========================
    ! Top row
    a(1) = 0.0_dp
    c(1) = - THETA * wcell(1) * Df(2)
    b(1) =   inv_dt - (a(1) + c(1))

    ! Interior rows
    do k = 2, nlay-1
      a(k) = - THETA * wcell(k) * Df(k)
      c(k) = - THETA * wcell(k) * Df(k+1)
      b(k) =   inv_dt - (a(k) + c(k))
    end do

    ! Bottom row: include bottom face on diagonal; no 'c' term (face uses q0)
    a(nlay) = - THETA * wcell(nlay) * Df(nlay)
    c(nlay) = 0.0_dp
    b(nlay) =   inv_dt - a(nlay) + THETA * wcell(nlay) * Df(nlay+1)

    ! =========================
    ! Factor matrix once
    ! =========================
    call tridiag_factor(a, b, c, bp, cp)

    ! =========================
    ! Solve for each tracer
    ! =========================
    ! You can parallelize across n (independent RHS):
    ! !$omp parallel do default(none) shared(nq,nlay,q,q0,Df,wcell,inv_dt,THETA,bp,cp) private(n,k,rhs)
    do n = 1, nq
      ! ---- Assemble RHS (depends on tracer q(:,n) and boundary q0(n)) ----
      rhs(1) = inv_dt*q(1,n) - (1.0_dp - THETA)*wcell(1) * ( &
                 -Df(2) * ( q(2,n) - q(1,n) ) )

      do k = 2, nlay-1
        rhs(k) = inv_dt*q(k,n) - (1.0_dp - THETA)*wcell(k) * ( &
                   -Df(k+1) * ( q(k+1,n) - q(k,n) ) + Df(k) * ( q(k,n) - q(k-1,n) ) )
      end do

      rhs(nlay) = inv_dt*q(nlay,n) - (1.0_dp - THETA)*wcell(nlay) * ( &
                     -Df(nlay+1)*( q0(n)     - q(nlay,n) ) + Df(nlay)*( q(nlay,n) - q(nlay-1,n) ) ) &
                   + THETA * wcell(nlay) * Df(nlay+1) * q0(n)

      ! ---- Solve with precomputed factors ----
      call tridiag_solve_fact(a, bp, cp, rhs, q(:,n))

      ! ---- Positivity floor ----
      do k = 1, nlay
        if (q(k,n) < qmin) q(k,n) = qmin
      end do
    end do
    ! !$omp end parallel do
  end subroutine vert_diff_imp


  !==========================
  ! Tridiagonal factor & solve
  !==========================
  ! Factorization stage (once): compute modified diagonal bp and upper ratio cp
  subroutine tridiag_factor(a, b, c, bp, cp)
    implicit none
    real(dp), intent(in)  :: a(:), b(:), c(:)
    real(dp), intent(out) :: bp(:), cp(:)
    integer :: n, i
    real(dp) :: denom

    n = size(b)
    bp(1) = b(1)
    if (abs(bp(1)) <= eps) error stop "tridiag_factor: zero pivot at i=1"
    cp(1) = c(1) / bp(1)

    do i = 2, n
      denom = b(i) - a(i) * cp(i-1)
      if (abs(denom) <= eps) error stop "tridiag_factor: zero pivot in forward sweep"
      bp(i) = denom
      if (i < n) then
        cp(i) = c(i) / bp(i)
      else
        cp(i) = 0.0_dp
      end if
    end do
  end subroutine tridiag_factor

  ! Solve stage (many RHS): uses precomputed bp, cp and original 'a' for forward sweep
  subroutine tridiag_solve_fact(a, bp, cp, d, x)
    implicit none
    real(dp), intent(in)  :: a(:), bp(:), cp(:), d(:)
    real(dp), intent(out) :: x(:)
    integer :: n, i
    real(dp) :: y

    n = size(bp)
    ! Forward substitution into x(:) to save a temp
    x(1) = d(1) / bp(1)
    do i = 2, n
      x(i) = ( d(i) - a(i) * x(i-1) ) / bp(i)
    end do
    ! Back substitution
    do i = n-1, 1, -1
      x(i) = x(i) - cp(i) * x(i+1)
    end do
  end subroutine tridiag_solve_fact


  pure elemental real(dp) function harm_mean(a, b) result(hm)
    real(dp), intent(in) :: a, b
    hm = 2.0_dp*a*b / (a + b + eps)
  end function harm_mean

end module vert_diff_imp_mod
