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
    character(len=20), allocatable, dimension(:) :: c_re, c_pr
    integer, allocatable, dimension(:) :: stoi_re, stoi_pr
    integer, allocatable, dimension(:) :: gi_re, gi_pr
    real(dp) :: dH, ds, Keq, f, r, net, dmu
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

  integer :: n_reac, n_sp
  type(reaction), allocatable, dimension(:) :: re
  type(species), allocatable, dimension(:) :: g_sp

  real(dp) :: nd_atm, P, P_cgs, T



end module mini_ch_class
