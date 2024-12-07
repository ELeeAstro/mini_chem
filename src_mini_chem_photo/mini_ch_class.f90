module mini_ch_class
  use mini_ch_precision
  implicit none

  type reaction
    integer :: id
    integer :: re_t
    real(dp) :: Tmin, Tmax
    integer :: n_re, n_pr
    real(dp) :: A, B, C
    real(dp) :: A0, B0, C0, Ainf, Binf, Cinf
    real(dp) :: rev_coeff
    character(len=50) :: fname
    character(len=20), allocatable, dimension(:) :: c_re, c_pr
    integer, allocatable, dimension(:) :: stoi_re, stoi_pr
    integer, allocatable, dimension(:) :: gi_re, gi_pr
    real(dp) :: dH, ds, Keq, f, r, net, dmu
    ! Reaction table data
    integer :: nT, nP, nkf
    real(dp), allocatable, dimension(:) :: T, P, lT, lP
    real(dp), allocatable, dimension(:,:) :: kf, lkf

    integer :: br_idx

  end type reaction

  type species
    integer :: id, n_a
    character(len=20) :: c
    real(dp) :: mw, nd
    real(dp), allocatable, dimension(:) :: a_l, a_h

    real(dp) :: thresh
    real(dp), allocatable, dimension(:) :: ph_xsec, Ray_xsec

  end type species

  real(dp), parameter :: kb = 1.380649e-16_dp
  real(dp), parameter :: R = 8.31446261815324e7_dp
  real(dp), parameter :: c_s = 2.99792458e10_dp
  real(dp), parameter :: h_p = 6.62607015e-27_dp
  real(dp), parameter :: P0 = 1.0e6_dp

  real(dp), parameter :: f_con = 0.001_dp
  real(dp), parameter :: del_con = 0.01_dp
  real(dp), parameter :: eps_con = 1.0e-4_dp

  integer :: n_reac, n_sp
  type(reaction), allocatable, dimension(:) :: re
  type(species), allocatable, dimension(:) :: g_sp

  real(dp), allocatable, dimension(:) :: Keq
  real(dp), allocatable, dimension(:) :: re_f, re_r

  integer :: nwl
  real(dp), allocatable, dimension(:) :: wl_grid, s_flux
  real(dp), allocatable, dimension(:,:) :: a_flux
  !$omp threadprivate(Keq, re_f, re_r)

contains

  subroutine locate(arr, n, var, idx)
    implicit none

    integer, intent(in) :: n
    integer, intent(out) :: idx
    real(dp), dimension(n), intent(in) :: arr
    real(dp),intent(in) ::  var
    integer :: jl, jm, ju

    ! Search an array using bi-section (numerical methods)
    ! Then return array index that is lower than var in arr

    jl = 0
    ju = n+1
    do while (ju-jl > 1)
      jm = (ju+jl)/2
      if ((arr(n) > arr(1)).eqv.(var > arr(jm))) then
        jl=jm
      else
        ju=jm
      end if
    end do

    idx = jl

  end subroutine locate

  subroutine linear_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(kind=dp), intent(in) :: xval, y1, y2, x1, x2
    real(kind=dp), intent(out) :: yval
    real(kind=dp) :: norm

    if (x1 >= x2) then
      print*, 'Error in linear_interp: x1 >= x2 - STOPPING', x1, x2
      stop
    end if

    norm = 1.0_dp / (x2 - x1)

    yval = (y1 * (x2 - xval) + y2 * (xval - x1)) * norm

  end subroutine linear_interp

  ! Perform linear interpolation in log10 space
  subroutine linear_log_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp) :: lxval, ly1, ly2, lx1, lx2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    ly1 = log10(y1); ly2 = log10(y2)

    norm = 1.0_dp / log10(x2/x1)

    yval = 10.0_dp**((ly1 * log10(x2/xval) + ly2 * log10(xval/x1)) * norm)

  end subroutine linear_log_interp

  subroutine bilinear_interp(xval, yval, x1, x2, y1, y2, a11, a21, a12, a22, aval)
    implicit none

    real(dp), intent(in) :: xval, yval, x1, x2, y1, y2, a11, a21, a12, a22
    real(dp), intent(out) :: aval
    real(dp) :: norm

    norm = 1.0_dp / (x2 - x1) / (y2 - y1)

    aval = a11 * (x2 - xval) * (y2 - yval) * norm &
      & + a21 * (xval - x1) * (y2 - yval) * norm &
      & + a12 * (x2 - xval) * (yval - y1) * norm &
      & + a22 * (xval - x1) * (yval - y1) * norm

  end subroutine bilinear_interp

  subroutine bilinear_log_interp(xval, yval, x1, x2, y1, y2, a11, a21, a12, a22, aval)
    implicit none

    real(dp), intent(in) :: xval, yval, x1, x2, y1, y2, a11, a21, a12, a22
    real(dp), intent(out) :: aval

    real(dp) :: lxval, lyval, lx1, lx2, ly1, ly2, la11, la21, la12, la22
    real(dp) :: norm

    lxval = log10(xval); lyval = log10(yval)
    lx1 = log10(x1);  lx2 = log10(x2)
    ly1 = log10(y1);  ly2 = log10(y2)
    la11 = log10(a11); la21 = log10(a21); la12 = log10(a12); la22 = log10(a22)

    norm = 1.0_dp / (lx2 - lx1) / (ly2 - ly1)

    aval = 10.0**(la11 * (lx2 - lxval) * (ly2 - lyval) * norm &
      & + la21 * (lxval - lx1) * (ly2 - lyval) * norm &
      & + la12 * (lx2 - lxval) * (lyval - ly1) * norm &
      & + la22 * (lxval - lx1) * (lyval - ly1) * norm)

  end subroutine bilinear_log_interp

  subroutine bezier_interp(xi, yi, ni, x, y)
    implicit none

    integer, intent(in) :: ni
    real(dp), dimension(ni), intent(in) :: xi, yi
    real(dp), intent(in) :: x
    real(dp), intent(out) :: y

    real(dp) :: dx, dx1, dy, dy1, w, yc, t, wlim, wlim1

    !xc = (xi(1) + xi(2))/2.0_dp ! Control point (no needed here, implicitly included)
    dx = xi(2) - xi(1)
    dx1 = xi(3) - xi(2)
    dy = yi(2) - yi(1)
    dy1 = yi(3) - yi(2)

    if (x > xi(1) .and. x < xi(2)) then
      ! left hand side interpolation
      !print*,'left'
      w = dx1/(dx + dx1)
      wlim = 1.0_dp + 1.0_dp/(1.0_dp - (dy1/dy) * (dx/dx1))
      wlim1 = 1.0_dp/(1.0_dp - (dy/dy1) * (dx1/dx))
      if (w <= min(wlim,wlim1) .or. w >= max(wlim,wlim1)) then
        w = 1.0_dp
      end if
      yc = yi(2) - dx/2.0_dp * (w*dy/dx + (1.0_dp - w)*dy1/dx1)
      t = (x - xi(1))/dx
      y = (1.0_dp - t)**2 * yi(1) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(2)
    else ! (x > xi(2) and x < xi(3)) then
      ! right hand side interpolation
      !print*,'right'
      w = dx/(dx + dx1)
      wlim = 1.0_dp/(1.0_dp - (dy1/dy) * (dx/dx1))
      wlim1 = 1.0_dp + 1.0_dp/(1.0_dp - (dy/dy1) * (dx1/dx))
      if (w <= min(wlim,wlim1) .or. w >= max(wlim,wlim1)) then
        w = 1.0_dp
      end if
      yc = yi(2) + dx1/2.0_dp * (w*dy1/dx1 + (1.0_dp - w)*dy/dx)
      t = (x - xi(2))/(dx1)
      y = (1.0_dp - t)**2 * yi(2) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(3)
    end if

  end subroutine bezier_interp

  pure function trapz(x, y) result(r)
    !! Calculates the integral of an array y with respect to x using the trapezoid
    !! approximation. Note that the mesh spacing of x does not have to be uniform.
    real(kind=dp), intent(in)  :: x(:)         !! Variable x
    real(kind=dp), intent(in)  :: y(size(x))   !! Function y(x)
    real(kind=dp)              :: r            !! Integral ∫y(x)·dx

    ! Integrate using the trapezoidal rule
    associate(n => size(x))
      r = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2.0_dp
    end associate
  end function trapz

end module mini_ch_class
