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
    real(dp), allocatable, dimension(:) :: T, P
    real(dp), allocatable, dimension(:,:) :: kf
  end type reaction

  type species
    integer :: id, n_a
    character(len=20) :: c
    real(dp) :: mw, nd
    real(dp), allocatable, dimension(:) :: a_l, a_h
    real(dp) :: H0, s0
  end type species

  real(dp), parameter :: kb = 1.380649e-16_dp
  real(dp), parameter :: R = 8.31446261815324e7_dp
  real(dp), parameter :: P0 = 1.0e6_dp

  real(dp), parameter :: f_con = 1.0_dp
  real(dp), parameter :: del_con = 0.01_dp
  real(dp), parameter :: eps_con = 1.0e-4_dp

  integer :: n_reac, n_sp
  type(reaction), allocatable, dimension(:) :: re
  type(species), allocatable, dimension(:) :: g_sp

  real(dp) :: nd_atm, P, P_cgs, T

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



end module mini_ch_class
