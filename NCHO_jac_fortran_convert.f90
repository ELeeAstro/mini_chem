subroutine jac_NCHO(N,X,Y,DFY,LDFY,RPAR,IPAR)
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
