'''
## Reaction ##

# Chemical Network without Photochemistry

#R1
OH + H2 -> H2O + H
#R3
OH + CO -> H + CO2
#R5
O + H2 -> OH + H
#R7
OH + C -> CO + H
#M9
H + H + M -> H2 + M
#S11
CH4 + H2O -> CO + H2 + H2 + H2
#S13
CH4 + CH4 -> C2H2 + H2 + H2 + H2
#S15
CH4 + CH4 -> C2H6 + H2
#S17
C + CO + H2 -> C2H2 + O


## Mapping ##

OH: y[0], H2: y[1], H2O: y[2], H: y[3], CO: y[4], CO2: y[5], O: y[6], C: y[7], CH4: y[8], C2H2: y[9], C2H6: y[10],

OH      0        -1*v_1(k, M, y[1], y[0], y[2], y[3]) -1*v_3(k, M, y[3], y[0], y[5], y[4]) +1*v_5(k, M, y[1], y[0], y[6], y[3]) -1*v_7(k, M, y[3], y[0], y[7], y[4])
H2      1        -1*v_1(k, M, y[1], y[0], y[2], y[3]) -1*v_5(k, M, y[1], y[0], y[6], y[3]) +1*v_9(k, M, y[3], y[1]) +1*v_11(k, M, y[8], y[2], y[1], y[4]) +1*v_11(k, M, y[8], y[2], y[1], y[4]) +1*v_11(k, M, y[8], y[2], y[1], y[4]) +1*v_13(k, M, y[8], y[9], y[1]) +1*v_13(k, M, y[8], y[9], y[1]) +1*v_13(k, M, y[8], y[9], y[1]) +1*v_15(k, M, y[8], y[10], y[1]) -1*v_17(k, M, y[9], y[4], y[7], y[1], y[6])
H2O     2        +1*v_1(k, M, y[1], y[0], y[2], y[3]) -1*v_11(k, M, y[8], y[2], y[1], y[4])
H       3        +1*v_1(k, M, y[1], y[0], y[2], y[3]) +1*v_3(k, M, y[3], y[0], y[5], y[4]) +1*v_5(k, M, y[1], y[0], y[6], y[3]) +1*v_7(k, M, y[3], y[0], y[7], y[4]) -1*v_9(k, M, y[3], y[1]) -1*v_9(k, M, y[3], y[1])
CO      4        -1*v_3(k, M, y[3], y[0], y[5], y[4]) +1*v_7(k, M, y[3], y[0], y[7], y[4]) +1*v_11(k, M, y[8], y[2], y[1], y[4]) -1*v_17(k, M, y[9], y[4], y[7], y[1], y[6])
CO2     5        +1*v_3(k, M, y[3], y[0], y[5], y[4])
O       6        -1*v_5(k, M, y[1], y[0], y[6], y[3]) +1*v_17(k, M, y[9], y[4], y[7], y[1], y[6])
C       7        -1*v_7(k, M, y[3], y[0], y[7], y[4]) -1*v_17(k, M, y[9], y[4], y[7], y[1], y[6])
CH4     8        -1*v_11(k, M, y[8], y[2], y[1], y[4]) -1*v_13(k, M, y[8], y[9], y[1]) -1*v_13(k, M, y[8], y[9], y[1]) -1*v_15(k, M, y[8], y[10], y[1]) -1*v_15(k, M, y[8], y[10], y[1])
C2H2    9        +1*v_13(k, M, y[8], y[9], y[1]) +1*v_17(k, M, y[9], y[4], y[7], y[1], y[6])
C2H6    10       +1*v_15(k, M, y[8], y[10], y[1])


    dfdy[indx[0], indx[0]] = -k[1]*y[:,1] - k[3]*y[:,4] - k[6]*y[:,3] - k[7]*y[:,7]
    dfdy[indx[0], indx[1]] = -k[1]*y[:,0] + k[5]*y[:,6]
    dfdy[indx[0], indx[2]] = k[2]*y[:,3]
    dfdy[indx[0], indx[3]] = k[2]*y[:,2] + k[4]*y[:,5] - k[6]*y[:,0] + k[8]*y[:,4]
    dfdy[indx[0], indx[4]] = -k[3]*y[:,0] + k[8]*y[:,3]
    dfdy[indx[0], indx[5]] = k[4]*y[:,3]
    dfdy[indx[0], indx[6]] = k[5]*y[:,1]
    dfdy[indx[0], indx[7]] = -k[7]*y[:,0]
    dfdy[indx[0], indx[8]] = 0
    dfdy[indx[0], indx[9]] = 0
    dfdy[indx[0], indx[10]] = 0
    dfdy[indx[1], indx[0]] = -k[1]*y[:,1] + k[6]*y[:,3]
    dfdy[indx[1], indx[1]] = -M*k[10] - 9*k[12]*y[:,1]**2*y[:,4] - 9*k[14]*y[:,1]**2*y[:,9] - k[16]*y[:,10] - k[17]*y[:,4]*y[:,7] - k[1]*y[:,0] - k[5]*y[:,6]
    dfdy[indx[1], indx[2]] = 3*k[11]*y[:,8] + k[2]*y[:,3]
    dfdy[indx[1], indx[3]] = 2*M*k[9]*y[:,3] + k[2]*y[:,2] + k[6]*y[:,0]
    dfdy[indx[1], indx[4]] = -3*k[12]*y[:,1]**3 - k[17]*y[:,1]*y[:,7]
    dfdy[indx[1], indx[5]] = 0
    dfdy[indx[1], indx[6]] = k[18]*y[:,9] - k[5]*y[:,1]
    dfdy[indx[1], indx[7]] = -k[17]*y[:,1]*y[:,4]
    dfdy[indx[1], indx[8]] = 3*k[11]*y[:,2] + 6*k[13]*y[:,8] + 2*k[15]*y[:,8]
    dfdy[indx[1], indx[9]] = -3*k[14]*y[:,1]**3 + k[18]*y[:,6]
    dfdy[indx[1], indx[10]] = -k[16]*y[:,1]
    dfdy[indx[2], indx[0]] = k[1]*y[:,1]
    dfdy[indx[2], indx[1]] = 3*k[12]*y[:,1]**2*y[:,4] + k[1]*y[:,0]
    dfdy[indx[2], indx[2]] = -k[11]*y[:,8] - k[2]*y[:,3]
    dfdy[indx[2], indx[3]] = -k[2]*y[:,2]
    dfdy[indx[2], indx[4]] = k[12]*y[:,1]**3
    dfdy[indx[2], indx[5]] = 0
    dfdy[indx[2], indx[6]] = 0
    dfdy[indx[2], indx[7]] = 0
    dfdy[indx[2], indx[8]] = -k[11]*y[:,2]
    dfdy[indx[2], indx[9]] = 0
    dfdy[indx[2], indx[10]] = 0
    dfdy[indx[3], indx[0]] = k[1]*y[:,1] + k[3]*y[:,4] - k[6]*y[:,3] + k[7]*y[:,7]
    dfdy[indx[3], indx[1]] = 2*M*k[10] + k[1]*y[:,0] + k[5]*y[:,6]
    dfdy[indx[3], indx[2]] = -k[2]*y[:,3]
    dfdy[indx[3], indx[3]] = -4*M*k[9]*y[:,3] - k[2]*y[:,2] - k[4]*y[:,5] - k[6]*y[:,0] - k[8]*y[:,4]
    dfdy[indx[3], indx[4]] = k[3]*y[:,0] - k[8]*y[:,3]
    dfdy[indx[3], indx[5]] = -k[4]*y[:,3]
    dfdy[indx[3], indx[6]] = k[5]*y[:,1]
    dfdy[indx[3], indx[7]] = k[7]*y[:,0]
    dfdy[indx[3], indx[8]] = 0
    dfdy[indx[3], indx[9]] = 0
    dfdy[indx[3], indx[10]] = 0
    dfdy[indx[4], indx[0]] = -k[3]*y[:,4] + k[7]*y[:,7]
    dfdy[indx[4], indx[1]] = -3*k[12]*y[:,1]**2*y[:,4] - k[17]*y[:,4]*y[:,7]
    dfdy[indx[4], indx[2]] = k[11]*y[:,8]
    dfdy[indx[4], indx[3]] = k[4]*y[:,5] - k[8]*y[:,4]
    dfdy[indx[4], indx[4]] = -k[12]*y[:,1]**3 - k[17]*y[:,1]*y[:,7] - k[3]*y[:,0] - k[8]*y[:,3]
    dfdy[indx[4], indx[5]] = k[4]*y[:,3]
    dfdy[indx[4], indx[6]] = k[18]*y[:,9]
    dfdy[indx[4], indx[7]] = -k[17]*y[:,1]*y[:,4] + k[7]*y[:,0]
    dfdy[indx[4], indx[8]] = k[11]*y[:,2]
    dfdy[indx[4], indx[9]] = k[18]*y[:,6]
    dfdy[indx[4], indx[10]] = 0
    dfdy[indx[5], indx[0]] = k[3]*y[:,4]
    dfdy[indx[5], indx[1]] = 0
    dfdy[indx[5], indx[2]] = 0
    dfdy[indx[5], indx[3]] = -k[4]*y[:,5]
    dfdy[indx[5], indx[4]] = k[3]*y[:,0]
    dfdy[indx[5], indx[5]] = -k[4]*y[:,3]
    dfdy[indx[5], indx[6]] = 0
    dfdy[indx[5], indx[7]] = 0
    dfdy[indx[5], indx[8]] = 0
    dfdy[indx[5], indx[9]] = 0
    dfdy[indx[5], indx[10]] = 0
    dfdy[indx[6], indx[0]] = k[6]*y[:,3]
    dfdy[indx[6], indx[1]] = k[17]*y[:,4]*y[:,7] - k[5]*y[:,6]
    dfdy[indx[6], indx[2]] = 0
    dfdy[indx[6], indx[3]] = k[6]*y[:,0]
    dfdy[indx[6], indx[4]] = k[17]*y[:,1]*y[:,7]
    dfdy[indx[6], indx[5]] = 0
    dfdy[indx[6], indx[6]] = -k[18]*y[:,9] - k[5]*y[:,1]
    dfdy[indx[6], indx[7]] = k[17]*y[:,1]*y[:,4]
    dfdy[indx[6], indx[8]] = 0
    dfdy[indx[6], indx[9]] = -k[18]*y[:,6]
    dfdy[indx[6], indx[10]] = 0
    dfdy[indx[7], indx[0]] = -k[7]*y[:,7]
    dfdy[indx[7], indx[1]] = -k[17]*y[:,4]*y[:,7]
    dfdy[indx[7], indx[2]] = 0
    dfdy[indx[7], indx[3]] = k[8]*y[:,4]
    dfdy[indx[7], indx[4]] = -k[17]*y[:,1]*y[:,7] + k[8]*y[:,3]
    dfdy[indx[7], indx[5]] = 0
    dfdy[indx[7], indx[6]] = k[18]*y[:,9]
    dfdy[indx[7], indx[7]] = -k[17]*y[:,1]*y[:,4] - k[7]*y[:,0]
    dfdy[indx[7], indx[8]] = 0
    dfdy[indx[7], indx[9]] = k[18]*y[:,6]
    dfdy[indx[7], indx[10]] = 0
    dfdy[indx[8], indx[0]] = 0
    dfdy[indx[8], indx[1]] = 3*k[12]*y[:,1]**2*y[:,4] + 6*k[14]*y[:,1]**2*y[:,9] + 2*k[16]*y[:,10]
    dfdy[indx[8], indx[2]] = -k[11]*y[:,8]
    dfdy[indx[8], indx[3]] = 0
    dfdy[indx[8], indx[4]] = k[12]*y[:,1]**3
    dfdy[indx[8], indx[5]] = 0
    dfdy[indx[8], indx[6]] = 0
    dfdy[indx[8], indx[7]] = 0
    dfdy[indx[8], indx[8]] = -k[11]*y[:,2] - 4*k[13]*y[:,8] - 4*k[15]*y[:,8]
    dfdy[indx[8], indx[9]] = 2*k[14]*y[:,1]**3
    dfdy[indx[8], indx[10]] = 2*k[16]*y[:,1]
    dfdy[indx[9], indx[0]] = 0
    dfdy[indx[9], indx[1]] = -3*k[14]*y[:,1]**2*y[:,9] + k[17]*y[:,4]*y[:,7]
    dfdy[indx[9], indx[2]] = 0
    dfdy[indx[9], indx[3]] = 0
    dfdy[indx[9], indx[4]] = k[17]*y[:,1]*y[:,7]
    dfdy[indx[9], indx[5]] = 0
    dfdy[indx[9], indx[6]] = -k[18]*y[:,9]
    dfdy[indx[9], indx[7]] = k[17]*y[:,1]*y[:,4]
    dfdy[indx[9], indx[8]] = 2*k[13]*y[:,8]
    dfdy[indx[9], indx[9]] = -k[14]*y[:,1]**3 - k[18]*y[:,6]
    dfdy[indx[9], indx[10]] = 0
    dfdy[indx[10], indx[0]] = 0
    dfdy[indx[10], indx[1]] = -k[16]*y[:,10]
    dfdy[indx[10], indx[2]] = 0
    dfdy[indx[10], indx[3]] = 0
    dfdy[indx[10], indx[4]] = 0
    dfdy[indx[10], indx[5]] = 0
    dfdy[indx[10], indx[6]] = 0
    dfdy[indx[10], indx[7]] = 0
    dfdy[indx[10], indx[8]] = 2*k[15]*y[:,8]
    dfdy[indx[10], indx[9]] = 0
    dfdy[indx[10], indx[10]] = -k[16]*y[:,1]


    # Fortran translation

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
