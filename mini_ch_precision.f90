module mini_ch_precision
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  !!!
  ! Different sets of single, double and quad precision availible for DMC
  ! Try different sets should one fail to compile/give errors for any reason
  ! This module should be compiled first in the DMC chain and used in every
  ! module/subroutine when required
  !!!

  private
  public :: sp, dp, qp


  ! Fortran 2008 intrinsic precisions - reccomonded if possible
  integer, parameter :: sp = REAL32
  integer, parameter :: dp = REAL64
  integer, parameter :: qp = REAL128

  ! Selected real kind precisions: selected_real_kind(decimal places, exoponent)
  !integer, parameter :: sp = selected_real_kind(6, 37)
  !integer, parameter :: dp = selected_real_kind(15, 307)
  !integer, parameter :: qp = selected_real_kind(33, 4931)

  ! Ensure machine precision for single and double precision
  !integer, parameter :: sp = kind(1.0)
  !integer, parameter :: dp = kind(1.0d0)

  ! Ensure dp and qp are double and quadrouple single precision
  !integer, parameter ::  sp = kind(1.0)
  !integer, parameter :: dp = selected_real_kind(2*precision(1.0_sp))
  !integer, parameter :: qp = selected_real_kind(2*precision(1.0_dp))


end module mini_ch_precision
