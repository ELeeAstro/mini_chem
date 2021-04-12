'''
## Reaction ##

# Chemical Network without Photochemistry

#R1
OH + H2 -> H2O + H
#R3
O + H2 -> OH + H
#M5
H + H + M -> H2 + M


## Mapping ##

OH: y[0], H2: y[1], H2O: y[2], H: y[3], O: y[4],

OH      0        -1*v_1(k, M, y[0], y[2], y[3], y[1]) +1*v_3(k, M, y[4], y[0], y[3], y[1])
H2      1        -1*v_1(k, M, y[0], y[2], y[3], y[1]) -1*v_3(k, M, y[4], y[0], y[3], y[1]) +1*v_5(k, M, y[3], y[1])
H2O     2        +1*v_1(k, M, y[0], y[2], y[3], y[1])
H       3        +1*v_1(k, M, y[0], y[2], y[3], y[1]) +1*v_3(k, M, y[4], y[0], y[3], y[1]) -1*v_5(k, M, y[3], y[1]) -1*v_5(k, M, y[3], y[1])
O       4        -1*v_3(k, M, y[4], y[0], y[3], y[1])
'''

dfdy[indx[0], indx[0]] = -k[1]*y[:,1] - k[4]*y[:,3]
dfdy[indx[0], indx[1]] = -k[1]*y[:,0] + k[3]*y[:,4]
dfdy[indx[0], indx[2]] = k[2]*y[:,3]
dfdy[indx[0], indx[3]] = k[2]*y[:,2] - k[4]*y[:,0]
dfdy[indx[0], indx[4]] = k[3]*y[:,1]
dfdy[indx[1], indx[0]] = -k[1]*y[:,1] + k[4]*y[:,3]
dfdy[indx[1], indx[1]] = -M*k[6] - k[1]*y[:,0] - k[3]*y[:,4]
dfdy[indx[1], indx[2]] = k[2]*y[:,3]
dfdy[indx[1], indx[3]] = 2*M*k[5]*y[:,3] + k[2]*y[:,2] + k[4]*y[:,0]
dfdy[indx[1], indx[4]] = -k[3]*y[:,1]
dfdy[indx[2], indx[0]] = k[1]*y[:,1]
dfdy[indx[2], indx[1]] = k[1]*y[:,0]
dfdy[indx[2], indx[2]] = -k[2]*y[:,3]
dfdy[indx[2], indx[3]] = -k[2]*y[:,2]
dfdy[indx[2], indx[4]] = 0
dfdy[indx[3], indx[0]] = k[1]*y[:,1] - k[4]*y[:,3]
dfdy[indx[3], indx[1]] = 2*M*k[6] + k[1]*y[:,0] + k[3]*y[:,4]
dfdy[indx[3], indx[2]] = -k[2]*y[:,3]
dfdy[indx[3], indx[3]] = -4*M*k[5]*y[:,3] - k[2]*y[:,2] - k[4]*y[:,0]
dfdy[indx[3], indx[4]] = k[3]*y[:,1]
dfdy[indx[4], indx[0]] = k[4]*y[:,3]
dfdy[indx[4], indx[1]] = -k[3]*y[:,4]
dfdy[indx[4], indx[2]] = 0
dfdy[indx[4], indx[3]] = k[4]*y[:,0]
dfdy[indx[4], indx[4]] = -k[3]*y[:,1]

#Fortran translation

dfy(1, 1) = -re(1)%f*y(2) - re(2)%r*y(4)
dfy(1, 2) = -re(1)%f*y(1) + re(2)%f*y(5)
dfy(1, 3) = re(1)%r*y(4)
dfy(1, 4) = re(1)%r*y(3) - re(2)%r*y(1)
dfy(1, 5) = re(2)%f*y(2)
dfy(2, 1) = -re(1)%f*y(2) + re(2)%r*y(4)
dfy(2, 2) = -nd_atm *re(3)%r - re(1)%f*y(1) - re(2)%f*y(5)
dfy(2, 3) = re(1)%r*y(4)
dfy(2, 4) = 2.0_dp*nd_atm*re(3)%f*y(4) + re(1)%r*y(3) + re(2)%r*y(1)
dfy(2, 5) = -re(2)%f*y(2)
dfy(3, 1) = re(1)%f*y(2)
dfy(3, 2) = re(1)%f*y(1)
dfy(3, 3) = -re(1)%r*y(4)
dfy(3, 4) = -re(1)%r*y(3)
dfy(3, 5) = 0.0_dp
dfy(4, 1) = re(1)%f*y(2) - re(2)%r*y(4)
dfy(4, 2) = 2.0_dp*nd_atm*re(3)%r + re(1)%f*y(1) + re(2)%f*y(5)
dfy(4, 3) = -re(1)%r*y(4)
dfy(4, 4) = -4.0_dp*nd_atm*re(3)%f*y(4) - re(1)%r*y(3) - re(2)%r*y(1)
dfy(4, 5) = re(2)%f*y(2)
dfy(5, 1) = re(2)%r*y(4)
dfy(5, 2) = -re(2)%f*y(5)
dfy(5, 3) = 0.0_dp
dfy(5, 4) = re(2)%r*y(1)
dfy(5, 5) = -re(2)%f*y(2)
