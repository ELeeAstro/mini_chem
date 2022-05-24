#!/usr/bin/python

from scipy import *
import numpy as np
from phy_const import kb, Navo
import vulcan_cfg

'''
## Reaction ##

# Chemical Network without Photochemistry

#R1
OH + H2 -> H2O + H
#R3
OH + CO -> H + CO2
#R5
O + H2 -> OH + H
#M7
H + H + M -> H2 + M
#S9
CH4 + H2O -> CO + H2 + H2 + H2
#S11
CH4 + CH4 -> C2H2 + H2 + H2 + H2
#S13
CO + CH4 -> C2H2 + H2O
#S15
NH3 + NH3 -> N2 + H2 + H2 + H2
#S17
CH4 + NH3 -> HCN + H2 + H2 + H2
#S19
CO + NH3 -> HCN + H2O


## Mapping ##

OH: y[0], H2: y[1], H2O: y[2], H: y[3], CO: y[4], CO2: y[5], O: y[6], CH4: y[7], C2H2: y[8], NH3: y[9], N2: y[10], HCN: y[11],

OH	0	 -1*v_1(k, M, y[3], y[2], y[0], y[1]) -1*v_3(k, M, y[5], y[4], y[3], y[0]) +1*v_5(k, M, y[6], y[3], y[0], y[1])
H2	1	 -1*v_1(k, M, y[3], y[2], y[0], y[1]) -1*v_5(k, M, y[6], y[3], y[0], y[1]) +1*v_7(k, M, y[3], y[1]) +1*v_9(k, M, y[4], y[2], y[7], y[1]) +1*v_9(k, M, y[4], y[2], y[7], y[1]) +1*v_9(k, M, y[4], y[2], y[7], y[1]) +1*v_11(k, M, y[8], y[7], y[1]) +1*v_11(k, M, y[8], y[7], y[1]) +1*v_11(k, M, y[8], y[7], y[1]) +1*v_15(k, M, y[10], y[1], y[9]) +1*v_15(k, M, y[10], y[1], y[9]) +1*v_15(k, M, y[10], y[1], y[9]) +1*v_17(k, M, y[11], y[1], y[7], y[9]) +1*v_17(k, M, y[11], y[1], y[7], y[9]) +1*v_17(k, M, y[11], y[1], y[7], y[9])
H2O	2	 +1*v_1(k, M, y[3], y[2], y[0], y[1]) -1*v_9(k, M, y[4], y[2], y[7], y[1]) +1*v_13(k, M, y[4], y[8], y[2], y[7]) +1*v_19(k, M, y[4], y[11], y[2], y[9])
H	3	 +1*v_1(k, M, y[3], y[2], y[0], y[1]) +1*v_3(k, M, y[5], y[4], y[3], y[0]) +1*v_5(k, M, y[6], y[3], y[0], y[1]) -1*v_7(k, M, y[3], y[1]) -1*v_7(k, M, y[3], y[1])
CO	4	 -1*v_3(k, M, y[5], y[4], y[3], y[0]) +1*v_9(k, M, y[4], y[2], y[7], y[1]) -1*v_13(k, M, y[4], y[8], y[2], y[7]) -1*v_19(k, M, y[4], y[11], y[2], y[9])
CO2	5	 +1*v_3(k, M, y[5], y[4], y[3], y[0])
O	6	 -1*v_5(k, M, y[6], y[3], y[0], y[1])
CH4	7	 -1*v_9(k, M, y[4], y[2], y[7], y[1]) -1*v_11(k, M, y[8], y[7], y[1]) -1*v_11(k, M, y[8], y[7], y[1]) -1*v_13(k, M, y[4], y[8], y[2], y[7]) -1*v_17(k, M, y[11], y[1], y[7], y[9])
C2H2	8	 +1*v_11(k, M, y[8], y[7], y[1]) +1*v_13(k, M, y[4], y[8], y[2], y[7])
NH3	9	 -1*v_15(k, M, y[10], y[1], y[9]) -1*v_15(k, M, y[10], y[1], y[9]) -1*v_17(k, M, y[11], y[1], y[7], y[9]) -1*v_19(k, M, y[4], y[11], y[2], y[9])
N2	10	 +1*v_15(k, M, y[10], y[1], y[9])
HCN	11	 +1*v_17(k, M, y[11], y[1], y[7], y[9]) +1*v_19(k, M, y[4], y[11], y[2], y[9])
'''

#species list
spec_list = ['OH', 'H2', 'H2O', 'H', 'CO', 'CO2', 'O', 'CH4', 'C2H2', 'NH3', 'N2', 'HCN']
# the total number of species
ni = 12
# the total number of reactions (forward and reverse)
nr = 20

# store the products and the reactants in the 1st and the 2nd element for reaction j (without M)
re_dict = {1:[['OH', 'H2'], ['H2O', 'H']], 2:[['H2O', 'H'], ['OH', 'H2']], 3:[['OH', 'CO'], ['H', 'CO2']], 4:[['H', 'CO2'], ['OH', 'CO']], 5:[['O', 'H2'], ['OH', 'H']], 6:[['OH', 'H'], ['O', 'H2']], 7:[['H', 'H'], ['H2']], 8:[['H2'], ['H', 'H']], 9:[['CH4', 'H2O'], ['CO', 'H2', 'H2', 'H2']], 10:[['CO', 'H2', 'H2', 'H2'], ['CH4', 'H2O']], 11:[['CH4', 'CH4'], ['C2H2', 'H2', 'H2', 'H2']], 12:[['C2H2', 'H2', 'H2', 'H2'], ['CH4', 'CH4']], 13:[['CO', 'CH4'], ['C2H2', 'H2O']], 14:[['C2H2', 'H2O'], ['CO', 'CH4']], 15:[['NH3', 'NH3'], ['N2', 'H2', 'H2', 'H2']], 16:[['N2', 'H2', 'H2', 'H2'], ['NH3', 'NH3']], 17:[['CH4', 'NH3'], ['HCN', 'H2', 'H2', 'H2']], 18:[['HCN', 'H2', 'H2', 'H2'], ['CH4', 'NH3']], 19:[['CO', 'NH3'], ['HCN', 'H2O']], 20:[['HCN', 'H2O'], ['CO', 'NH3']]}


# store the products and the reactants in the 1st and the 2nd element for reaction j (with M)
re_wM_dict = {1:[['OH', 'H2'], ['H2O', 'H']], 2:[['H2O', 'H'], ['OH', 'H2']], 3:[['OH', 'CO'], ['H', 'CO2']], 4:[['H', 'CO2'], ['OH', 'CO']], 5:[['O', 'H2'], ['OH', 'H']], 6:[['OH', 'H'], ['O', 'H2']], 7:[['H', 'H', 'M'], ['H2', 'M']], 8:[['H2', 'M'], ['H', 'H', 'M']], 9:[['CH4', 'H2O'], ['CO', 'H2', 'H2', 'H2']], 10:[['CO', 'H2', 'H2', 'H2'], ['CH4', 'H2O']], 11:[['CH4', 'CH4'], ['C2H2', 'H2', 'H2', 'H2']], 12:[['C2H2', 'H2', 'H2', 'H2'], ['CH4', 'CH4']], 13:[['CO', 'CH4'], ['C2H2', 'H2O']], 14:[['C2H2', 'H2O'], ['CO', 'CH4']], 15:[['NH3', 'NH3'], ['N2', 'H2', 'H2', 'H2']], 16:[['N2', 'H2', 'H2', 'H2'], ['NH3', 'NH3']], 17:[['CH4', 'NH3'], ['HCN', 'H2', 'H2', 'H2']], 18:[['HCN', 'H2', 'H2', 'H2'], ['CH4', 'NH3']], 19:[['CO', 'NH3'], ['HCN', 'H2O']], 20:[['HCN', 'H2O'], ['CO', 'NH3']]}


def chemdf(y, M, k):
    y = np.transpose(y)
    dydt = np.zeros(shape=y.shape)
    dydt[0] =  -1*v_1(k, M, y[3], y[2], y[0], y[1]) -1*v_3(k, M, y[5], y[4], y[3], y[0]) +1*v_5(k, M, y[6], y[3], y[0], y[1])
    dydt[1] =  -1*v_1(k, M, y[3], y[2], y[0], y[1]) -1*v_5(k, M, y[6], y[3], y[0], y[1]) +1*v_7(k, M, y[3], y[1]) +1*v_9(k, M, y[4], y[2], y[7], y[1]) +1*v_9(k, M, y[4], y[2], y[7], y[1]) +1*v_9(k, M, y[4], y[2], y[7], y[1]) +1*v_11(k, M, y[8], y[7], y[1]) +1*v_11(k, M, y[8], y[7], y[1]) +1*v_11(k, M, y[8], y[7], y[1]) +1*v_15(k, M, y[10], y[1], y[9]) +1*v_15(k, M, y[10], y[1], y[9]) +1*v_15(k, M, y[10], y[1], y[9]) +1*v_17(k, M, y[11], y[1], y[7], y[9]) +1*v_17(k, M, y[11], y[1], y[7], y[9]) +1*v_17(k, M, y[11], y[1], y[7], y[9])
    dydt[2] =  +1*v_1(k, M, y[3], y[2], y[0], y[1]) -1*v_9(k, M, y[4], y[2], y[7], y[1]) +1*v_13(k, M, y[4], y[8], y[2], y[7]) +1*v_19(k, M, y[4], y[11], y[2], y[9])
    dydt[3] =  +1*v_1(k, M, y[3], y[2], y[0], y[1]) +1*v_3(k, M, y[5], y[4], y[3], y[0]) +1*v_5(k, M, y[6], y[3], y[0], y[1]) -1*v_7(k, M, y[3], y[1]) -1*v_7(k, M, y[3], y[1])
    dydt[4] =  -1*v_3(k, M, y[5], y[4], y[3], y[0]) +1*v_9(k, M, y[4], y[2], y[7], y[1]) -1*v_13(k, M, y[4], y[8], y[2], y[7]) -1*v_19(k, M, y[4], y[11], y[2], y[9])
    dydt[5] =  +1*v_3(k, M, y[5], y[4], y[3], y[0])
    dydt[6] =  -1*v_5(k, M, y[6], y[3], y[0], y[1])
    dydt[7] =  -1*v_9(k, M, y[4], y[2], y[7], y[1]) -1*v_11(k, M, y[8], y[7], y[1]) -1*v_11(k, M, y[8], y[7], y[1]) -1*v_13(k, M, y[4], y[8], y[2], y[7]) -1*v_17(k, M, y[11], y[1], y[7], y[9])
    dydt[8] =  +1*v_11(k, M, y[8], y[7], y[1]) +1*v_13(k, M, y[4], y[8], y[2], y[7])
    dydt[9] =  -1*v_15(k, M, y[10], y[1], y[9]) -1*v_15(k, M, y[10], y[1], y[9]) -1*v_17(k, M, y[11], y[1], y[7], y[9]) -1*v_19(k, M, y[4], y[11], y[2], y[9])
    dydt[10] =  +1*v_15(k, M, y[10], y[1], y[9])
    dydt[11] =  +1*v_17(k, M, y[11], y[1], y[7], y[9]) +1*v_19(k, M, y[4], y[11], y[2], y[9])
    dydt = np.transpose(dydt)
    return dydt

def df(y, M, k):
    df_list = []
    df_list.append(  -(k[1]*y[0]*y[1] - k[2]*y[2]*y[3]) -(k[3]*y[0]*y[4] - k[4]*y[3]*y[5]) +(k[5]*y[6]*y[1] - k[6]*y[0]*y[3]) )
    df_list.append(  -(k[1]*y[0]*y[1] - k[2]*y[2]*y[3]) -(k[5]*y[6]*y[1] - k[6]*y[0]*y[3]) +(k[7]*y[3]*y[3]*M - k[8]*y[1]*M) +(k[9]*y[7]*y[2] - k[10]*y[4]*y[1]*y[1]*y[1]) +(k[9]*y[7]*y[2] - k[10]*y[4]*y[1]*y[1]*y[1]) +(k[9]*y[7]*y[2] - k[10]*y[4]*y[1]*y[1]*y[1]) +(k[11]*y[7]*y[7] - k[12]*y[8]*y[1]*y[1]*y[1]) +(k[11]*y[7]*y[7] - k[12]*y[8]*y[1]*y[1]*y[1]) +(k[11]*y[7]*y[7] - k[12]*y[8]*y[1]*y[1]*y[1]) +(k[15]*y[9]*y[9] - k[16]*y[10]*y[1]*y[1]*y[1]) +(k[15]*y[9]*y[9] - k[16]*y[10]*y[1]*y[1]*y[1]) +(k[15]*y[9]*y[9] - k[16]*y[10]*y[1]*y[1]*y[1]) +(k[17]*y[7]*y[9] - k[18]*y[11]*y[1]*y[1]*y[1]) +(k[17]*y[7]*y[9] - k[18]*y[11]*y[1]*y[1]*y[1]) +(k[17]*y[7]*y[9] - k[18]*y[11]*y[1]*y[1]*y[1]) )
    df_list.append(  +(k[1]*y[0]*y[1] - k[2]*y[2]*y[3]) -(k[9]*y[7]*y[2] - k[10]*y[4]*y[1]*y[1]*y[1]) +(k[13]*y[4]*y[7] - k[14]*y[8]*y[2]) +(k[19]*y[4]*y[9] - k[20]*y[11]*y[2]) )
    df_list.append(  +(k[1]*y[0]*y[1] - k[2]*y[2]*y[3]) +(k[3]*y[0]*y[4] - k[4]*y[3]*y[5]) +(k[5]*y[6]*y[1] - k[6]*y[0]*y[3]) -(k[7]*y[3]*y[3]*M - k[8]*y[1]*M) -(k[7]*y[3]*y[3]*M - k[8]*y[1]*M) )
    df_list.append(  -(k[3]*y[0]*y[4] - k[4]*y[3]*y[5]) +(k[9]*y[7]*y[2] - k[10]*y[4]*y[1]*y[1]*y[1]) -(k[13]*y[4]*y[7] - k[14]*y[8]*y[2]) -(k[19]*y[4]*y[9] - k[20]*y[11]*y[2]) )
    df_list.append(  +(k[3]*y[0]*y[4] - k[4]*y[3]*y[5]) )
    df_list.append(  -(k[5]*y[6]*y[1] - k[6]*y[0]*y[3]) )
    df_list.append(  -(k[9]*y[7]*y[2] - k[10]*y[4]*y[1]*y[1]*y[1]) -(k[11]*y[7]*y[7] - k[12]*y[8]*y[1]*y[1]*y[1]) -(k[11]*y[7]*y[7] - k[12]*y[8]*y[1]*y[1]*y[1]) -(k[13]*y[4]*y[7] - k[14]*y[8]*y[2]) -(k[17]*y[7]*y[9] - k[18]*y[11]*y[1]*y[1]*y[1]) )
    df_list.append(  +(k[11]*y[7]*y[7] - k[12]*y[8]*y[1]*y[1]*y[1]) +(k[13]*y[4]*y[7] - k[14]*y[8]*y[2]) )
    df_list.append(  -(k[15]*y[9]*y[9] - k[16]*y[10]*y[1]*y[1]*y[1]) -(k[15]*y[9]*y[9] - k[16]*y[10]*y[1]*y[1]*y[1]) -(k[17]*y[7]*y[9] - k[18]*y[11]*y[1]*y[1]*y[1]) -(k[19]*y[4]*y[9] - k[20]*y[11]*y[2]) )
    df_list.append(  +(k[15]*y[9]*y[9] - k[16]*y[10]*y[1]*y[1]*y[1]) )
    df_list.append(  +(k[17]*y[7]*y[9] - k[18]*y[11]*y[1]*y[1]*y[1]) +(k[19]*y[4]*y[9] - k[20]*y[11]*y[2]) )
    return df_list

#OH + H2 -> H2O + H
v_1 = lambda k, M, H, H2O, OH, H2 : k[1]*OH*H2 - k[2]*H2O*H


#OH + CO -> H + CO2
v_3 = lambda k, M, CO2, CO, H, OH : k[3]*OH*CO - k[4]*H*CO2


#O + H2 -> OH + H
v_5 = lambda k, M, O, H, OH, H2 : k[5]*O*H2 - k[6]*OH*H


#H + H + M -> H2 + M
v_7 = lambda k, M, H, H2 : k[7]*H*H*M - k[8]*H2*M


#CH4 + H2O -> CO + H2 + H2 + H2
v_9 = lambda k, M, CO, H2O, CH4, H2 : k[9]*CH4*H2O - k[10]*CO*H2*H2*H2


#CH4 + CH4 -> C2H2 + H2 + H2 + H2
v_11 = lambda k, M, C2H2, CH4, H2 : k[11]*CH4*CH4 - k[12]*C2H2*H2*H2*H2


#CO + CH4 -> C2H2 + H2O
v_13 = lambda k, M, CO, C2H2, H2O, CH4 : k[13]*CO*CH4 - k[14]*C2H2*H2O


#NH3 + NH3 -> N2 + H2 + H2 + H2
v_15 = lambda k, M, N2, H2, NH3 : k[15]*NH3*NH3 - k[16]*N2*H2*H2*H2


#CH4 + NH3 -> HCN + H2 + H2 + H2
v_17 = lambda k, M, HCN, H2, CH4, NH3 : k[17]*CH4*NH3 - k[18]*HCN*H2*H2*H2


#CO + NH3 -> HCN + H2O
v_19 = lambda k, M, CO, HCN, H2O, NH3 : k[19]*CO*NH3 - k[20]*HCN*H2O


def rate_ans(sp):
    rate_str = {}
    re_sp_dic = {}
    rate_str["OH"] = [ -1*k[1]*y[0]*y[1], +1* k[2]*y[2]*y[3], -1*k[3]*y[0]*y[4], +1* k[4]*y[3]*y[5], +1*k[5]*y[6]*y[1], -1* k[6]*y[0]*y[3]]
    rate_str["H2"] = [ -1*k[1]*y[0]*y[1], +1* k[2]*y[2]*y[3], -1*k[5]*y[6]*y[1], +1* k[6]*y[0]*y[3], +1*k[7]*y[3]*y[3]*M, -1* k[8]*y[1]*M, +1*k[9]*y[7]*y[2], -1* k[10]*y[4]*y[1]*y[1]*y[1], +1*k[9]*y[7]*y[2], -1* k[10]*y[4]*y[1]*y[1]*y[1], +1*k[9]*y[7]*y[2], -1* k[10]*y[4]*y[1]*y[1]*y[1], +1*k[11]*y[7]*y[7], -1* k[12]*y[8]*y[1]*y[1]*y[1], +1*k[11]*y[7]*y[7], -1* k[12]*y[8]*y[1]*y[1]*y[1], +1*k[11]*y[7]*y[7], -1* k[12]*y[8]*y[1]*y[1]*y[1], +1*k[15]*y[9]*y[9], -1* k[16]*y[10]*y[1]*y[1]*y[1], +1*k[15]*y[9]*y[9], -1* k[16]*y[10]*y[1]*y[1]*y[1], +1*k[15]*y[9]*y[9], -1* k[16]*y[10]*y[1]*y[1]*y[1], +1*k[17]*y[7]*y[9], -1* k[18]*y[11]*y[1]*y[1]*y[1], +1*k[17]*y[7]*y[9], -1* k[18]*y[11]*y[1]*y[1]*y[1], +1*k[17]*y[7]*y[9], -1* k[18]*y[11]*y[1]*y[1]*y[1]]
    rate_str["H2O"] = [ +1*k[1]*y[0]*y[1], -1* k[2]*y[2]*y[3], -1*k[9]*y[7]*y[2], +1* k[10]*y[4]*y[1]*y[1]*y[1], +1*k[13]*y[4]*y[7], -1* k[14]*y[8]*y[2], +1*k[19]*y[4]*y[9], -1* k[20]*y[11]*y[2]]
    rate_str["H"] = [ +1*k[1]*y[0]*y[1], -1* k[2]*y[2]*y[3], +1*k[3]*y[0]*y[4], -1* k[4]*y[3]*y[5], +1*k[5]*y[6]*y[1], -1* k[6]*y[0]*y[3], -1*k[7]*y[3]*y[3]*M, +1* k[8]*y[1]*M, -1*k[7]*y[3]*y[3]*M, +1* k[8]*y[1]*M]
    rate_str["CO"] = [ -1*k[3]*y[0]*y[4], +1* k[4]*y[3]*y[5], +1*k[9]*y[7]*y[2], -1* k[10]*y[4]*y[1]*y[1]*y[1], -1*k[13]*y[4]*y[7], +1* k[14]*y[8]*y[2], -1*k[19]*y[4]*y[9], +1* k[20]*y[11]*y[2]]
    rate_str["CO2"] = [ +1*k[3]*y[0]*y[4], -1* k[4]*y[3]*y[5]]
    rate_str["O"] = [ -1*k[5]*y[6]*y[1], +1* k[6]*y[0]*y[3]]
    rate_str["CH4"] = [ -1*k[9]*y[7]*y[2], +1* k[10]*y[4]*y[1]*y[1]*y[1], -1*k[11]*y[7]*y[7], +1* k[12]*y[8]*y[1]*y[1]*y[1], -1*k[11]*y[7]*y[7], +1* k[12]*y[8]*y[1]*y[1]*y[1], -1*k[13]*y[4]*y[7], +1* k[14]*y[8]*y[2], -1*k[17]*y[7]*y[9], +1* k[18]*y[11]*y[1]*y[1]*y[1]]
    rate_str["C2H2"] = [ +1*k[11]*y[7]*y[7], -1* k[12]*y[8]*y[1]*y[1]*y[1], +1*k[13]*y[4]*y[7], -1* k[14]*y[8]*y[2]]
    rate_str["NH3"] = [ -1*k[15]*y[9]*y[9], +1* k[16]*y[10]*y[1]*y[1]*y[1], -1*k[15]*y[9]*y[9], +1* k[16]*y[10]*y[1]*y[1]*y[1], -1*k[17]*y[7]*y[9], +1* k[18]*y[11]*y[1]*y[1]*y[1], -1*k[19]*y[4]*y[9], +1* k[20]*y[11]*y[2]]
    rate_str["N2"] = [ +1*k[15]*y[9]*y[9], -1* k[16]*y[10]*y[1]*y[1]*y[1]]
    rate_str["HCN"] = [ +1*k[17]*y[7]*y[9], -1* k[18]*y[11]*y[1]*y[1]*y[1], +1*k[19]*y[4]*y[9], -1* k[20]*y[11]*y[2]]
    return np.array(rate_str[sp])

corr = kb/1.e6  #P0=1.e6

# the data of 'H2CO' is from Brucat's 2015
#C2H NASA 9 new from Brucat
nasa9 = {}
for i in [ _ for _ in spec_list]:
    nasa9[i] = np.loadtxt('thermo/NASA9/' + str(i) + '.txt')
    nasa9[i] = nasa9[i].flatten()
    nasa9[i,'low'] = nasa9[i][0:10]
    nasa9[i,'high'] = nasa9[i][10:20]

#H/RT
def h_RT(T,a):
    return -a[0]/T**2 + a[1]*np.log(T)/T + a[2] + a[3]*T/2. + a[4]*T**2/3. + a[5]*T**3/4. + a[6]*T**4/5. + a[8]/T

#s/R
def s_R(T,a):
    return -a[0]/T**2/2. -a[1]/T + a[2]*np.log(T) + a[3]*T + a[4]*T**2/2. + a[5]*T**3/3. + a[6]*T**4/4. + a[9]

#g/RT (non-dimensional)
def g_RT(T,a_low,a_high): # T has to be 200 < T < 6000 K
    gi = (T < 1000)*(h_RT(T, a_low)-s_R(T, a_low)) + (T >= 1000)*(h_RT(T, a_high)-s_R(T, a_high))

    return gi

def gibbs_sp(i,T):
    gi = g_RT(T,nasa9[i,'low'],nasa9[i,'high'])
    return gi

#cp/R
def cp_R(T,a):
    return a[0]/T**2 + a[1]/T + a[2] + a[3]*T + a[4]*T**2 + a[5]*T**3 + a[6]*T**4

#cp/R of species sp
def cp_R_sp(i,T):
    if np.any(np.logical_or(T < 200, T > 6000)):
        print ('T exceeds the valid range.')
    cp =  (T < 1000)*cp_R(T,nasa9[i,'low']) + (T >= 1000)*cp_R(T,nasa9[i,'high'])
    return cp

# Gibbs free energy:
def Gibbs(i,T):
    G={}
    G[1] = lambda T: np.exp( -(-1*gibbs_sp('OH',T)-1*gibbs_sp('H2',T)+1*gibbs_sp('H2O',T)+1*gibbs_sp('H',T) ) )
    G[3] = lambda T: np.exp( -(-1*gibbs_sp('OH',T)-1*gibbs_sp('CO',T)+1*gibbs_sp('H',T)+1*gibbs_sp('CO2',T) ) )
    G[5] = lambda T: np.exp( -(-1*gibbs_sp('O',T)-1*gibbs_sp('H2',T)+1*gibbs_sp('OH',T)+1*gibbs_sp('H',T) ) )
    G[7] = lambda T: np.exp( -(-1*gibbs_sp('H',T)-1*gibbs_sp('H',T)+1*gibbs_sp('H2',T) ) )*(corr*T)**1
    G[9] = lambda T: np.exp( -(-1*gibbs_sp('CH4',T)-1*gibbs_sp('H2O',T)+1*gibbs_sp('CO',T)+1*gibbs_sp('H2',T)+1*gibbs_sp('H2',T)+1*gibbs_sp('H2',T) ) )*(corr*T)**-2
    G[11] = lambda T: np.exp( -(-1*gibbs_sp('CH4',T)-1*gibbs_sp('CH4',T)+1*gibbs_sp('C2H2',T)+1*gibbs_sp('H2',T)+1*gibbs_sp('H2',T)+1*gibbs_sp('H2',T) ) )*(corr*T)**-2
    G[13] = lambda T: np.exp( -(-1*gibbs_sp('CO',T)-1*gibbs_sp('CH4',T)+1*gibbs_sp('C2H2',T)+1*gibbs_sp('H2O',T) ) )
    G[15] = lambda T: np.exp( -(-1*gibbs_sp('NH3',T)-1*gibbs_sp('NH3',T)+1*gibbs_sp('N2',T)+1*gibbs_sp('H2',T)+1*gibbs_sp('H2',T)+1*gibbs_sp('H2',T) ) )*(corr*T)**-2
    G[17] = lambda T: np.exp( -(-1*gibbs_sp('CH4',T)-1*gibbs_sp('NH3',T)+1*gibbs_sp('HCN',T)+1*gibbs_sp('H2',T)+1*gibbs_sp('H2',T)+1*gibbs_sp('H2',T) ) )*(corr*T)**-2
    G[19] = lambda T: np.exp( -(-1*gibbs_sp('CO',T)-1*gibbs_sp('NH3',T)+1*gibbs_sp('HCN',T)+1*gibbs_sp('H2O',T) ) )
    return G[i](T)


def symjac(y, M, k):
    nz = vulcan_cfg.nz
    dfdy = np.zeros(shape=[ni*nz, ni*nz])
    indx = []
    for j in range(ni):
        indx.append( np.arange(j,j+ni*nz,ni) )
    dfdy[indx[0], indx[0]] = -k[1]*y[:,1] - k[3]*y[:,4] - k[6]*y[:,3]
    dfdy[indx[0], indx[1]] = -k[1]*y[:,0] + k[5]*y[:,6]
    dfdy[indx[0], indx[2]] = k[2]*y[:,3]
    dfdy[indx[0], indx[3]] = k[2]*y[:,2] + k[4]*y[:,5] - k[6]*y[:,0]
    dfdy[indx[0], indx[4]] = -k[3]*y[:,0]
    dfdy[indx[0], indx[5]] = k[4]*y[:,3]
    dfdy[indx[0], indx[6]] = k[5]*y[:,1]
    dfdy[indx[0], indx[7]] = 0
    dfdy[indx[0], indx[8]] = 0
    dfdy[indx[0], indx[9]] = 0
    dfdy[indx[0], indx[10]] = 0
    dfdy[indx[0], indx[11]] = 0
    dfdy[indx[1], indx[0]] = -k[1]*y[:,1] + k[6]*y[:,3]
    dfdy[indx[1], indx[1]] = -M*k[8] - 9*k[10]*y[:,1]**2*y[:,4] - 9*k[12]*y[:,1]**2*y[:,8] - 9*k[16]*y[:,10]*y[:,1]**2 - 9*k[18]*y[:,11]*y[:,1]**2 - k[1]*y[:,0] - k[5]*y[:,6]
    dfdy[indx[1], indx[2]] = k[2]*y[:,3] + 3*k[9]*y[:,7]
    dfdy[indx[1], indx[3]] = 2*M*k[7]*y[:,3] + k[2]*y[:,2] + k[6]*y[:,0]
    dfdy[indx[1], indx[4]] = -3*k[10]*y[:,1]**3
    dfdy[indx[1], indx[5]] = 0
    dfdy[indx[1], indx[6]] = -k[5]*y[:,1]
    dfdy[indx[1], indx[7]] = 6*k[11]*y[:,7] + 3*k[17]*y[:,9] + 3*k[9]*y[:,2]
    dfdy[indx[1], indx[8]] = -3*k[12]*y[:,1]**3
    dfdy[indx[1], indx[9]] = 6*k[15]*y[:,9] + 3*k[17]*y[:,7]
    dfdy[indx[1], indx[10]] = -3*k[16]*y[:,1]**3
    dfdy[indx[1], indx[11]] = -3*k[18]*y[:,1]**3
    dfdy[indx[2], indx[0]] = k[1]*y[:,1]
    dfdy[indx[2], indx[1]] = 3*k[10]*y[:,1]**2*y[:,4] + k[1]*y[:,0]
    dfdy[indx[2], indx[2]] = -k[14]*y[:,8] - k[20]*y[:,11] - k[2]*y[:,3] - k[9]*y[:,7]
    dfdy[indx[2], indx[3]] = -k[2]*y[:,2]
    dfdy[indx[2], indx[4]] = k[10]*y[:,1]**3 + k[13]*y[:,7] + k[19]*y[:,9]
    dfdy[indx[2], indx[5]] = 0
    dfdy[indx[2], indx[6]] = 0
    dfdy[indx[2], indx[7]] = k[13]*y[:,4] - k[9]*y[:,2]
    dfdy[indx[2], indx[8]] = -k[14]*y[:,2]
    dfdy[indx[2], indx[9]] = k[19]*y[:,4]
    dfdy[indx[2], indx[10]] = 0
    dfdy[indx[2], indx[11]] = -k[20]*y[:,2]
    dfdy[indx[3], indx[0]] = k[1]*y[:,1] + k[3]*y[:,4] - k[6]*y[:,3]
    dfdy[indx[3], indx[1]] = 2*M*k[8] + k[1]*y[:,0] + k[5]*y[:,6]
    dfdy[indx[3], indx[2]] = -k[2]*y[:,3]
    dfdy[indx[3], indx[3]] = -4*M*k[7]*y[:,3] - k[2]*y[:,2] - k[4]*y[:,5] - k[6]*y[:,0]
    dfdy[indx[3], indx[4]] = k[3]*y[:,0]
    dfdy[indx[3], indx[5]] = -k[4]*y[:,3]
    dfdy[indx[3], indx[6]] = k[5]*y[:,1]
    dfdy[indx[3], indx[7]] = 0
    dfdy[indx[3], indx[8]] = 0
    dfdy[indx[3], indx[9]] = 0
    dfdy[indx[3], indx[10]] = 0
    dfdy[indx[3], indx[11]] = 0
    dfdy[indx[4], indx[0]] = -k[3]*y[:,4]
    dfdy[indx[4], indx[1]] = -3*k[10]*y[:,1]**2*y[:,4]
    dfdy[indx[4], indx[2]] = k[14]*y[:,8] + k[20]*y[:,11] + k[9]*y[:,7]
    dfdy[indx[4], indx[3]] = k[4]*y[:,5]
    dfdy[indx[4], indx[4]] = -k[10]*y[:,1]**3 - k[13]*y[:,7] - k[19]*y[:,9] - k[3]*y[:,0]
    dfdy[indx[4], indx[5]] = k[4]*y[:,3]
    dfdy[indx[4], indx[6]] = 0
    dfdy[indx[4], indx[7]] = -k[13]*y[:,4] + k[9]*y[:,2]
    dfdy[indx[4], indx[8]] = k[14]*y[:,2]
    dfdy[indx[4], indx[9]] = -k[19]*y[:,4]
    dfdy[indx[4], indx[10]] = 0
    dfdy[indx[4], indx[11]] = k[20]*y[:,2]
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
    dfdy[indx[5], indx[11]] = 0
    dfdy[indx[6], indx[0]] = k[6]*y[:,3]
    dfdy[indx[6], indx[1]] = -k[5]*y[:,6]
    dfdy[indx[6], indx[2]] = 0
    dfdy[indx[6], indx[3]] = k[6]*y[:,0]
    dfdy[indx[6], indx[4]] = 0
    dfdy[indx[6], indx[5]] = 0
    dfdy[indx[6], indx[6]] = -k[5]*y[:,1]
    dfdy[indx[6], indx[7]] = 0
    dfdy[indx[6], indx[8]] = 0
    dfdy[indx[6], indx[9]] = 0
    dfdy[indx[6], indx[10]] = 0
    dfdy[indx[6], indx[11]] = 0
    dfdy[indx[7], indx[0]] = 0
    dfdy[indx[7], indx[1]] = 3*k[10]*y[:,1]**2*y[:,4] + 6*k[12]*y[:,1]**2*y[:,8] + 3*k[18]*y[:,11]*y[:,1]**2
    dfdy[indx[7], indx[2]] = k[14]*y[:,8] - k[9]*y[:,7]
    dfdy[indx[7], indx[3]] = 0
    dfdy[indx[7], indx[4]] = k[10]*y[:,1]**3 - k[13]*y[:,7]
    dfdy[indx[7], indx[5]] = 0
    dfdy[indx[7], indx[6]] = 0
    dfdy[indx[7], indx[7]] = -4*k[11]*y[:,7] - k[13]*y[:,4] - k[17]*y[:,9] - k[9]*y[:,2]
    dfdy[indx[7], indx[8]] = 2*k[12]*y[:,1]**3 + k[14]*y[:,2]
    dfdy[indx[7], indx[9]] = -k[17]*y[:,7]
    dfdy[indx[7], indx[10]] = 0
    dfdy[indx[7], indx[11]] = k[18]*y[:,1]**3
    dfdy[indx[8], indx[0]] = 0
    dfdy[indx[8], indx[1]] = -3*k[12]*y[:,1]**2*y[:,8]
    dfdy[indx[8], indx[2]] = -k[14]*y[:,8]
    dfdy[indx[8], indx[3]] = 0
    dfdy[indx[8], indx[4]] = k[13]*y[:,7]
    dfdy[indx[8], indx[5]] = 0
    dfdy[indx[8], indx[6]] = 0
    dfdy[indx[8], indx[7]] = 2*k[11]*y[:,7] + k[13]*y[:,4]
    dfdy[indx[8], indx[8]] = -k[12]*y[:,1]**3 - k[14]*y[:,2]
    dfdy[indx[8], indx[9]] = 0
    dfdy[indx[8], indx[10]] = 0
    dfdy[indx[8], indx[11]] = 0
    dfdy[indx[9], indx[0]] = 0
    dfdy[indx[9], indx[1]] = 6*k[16]*y[:,10]*y[:,1]**2 + 3*k[18]*y[:,11]*y[:,1]**2
    dfdy[indx[9], indx[2]] = k[20]*y[:,11]
    dfdy[indx[9], indx[3]] = 0
    dfdy[indx[9], indx[4]] = -k[19]*y[:,9]
    dfdy[indx[9], indx[5]] = 0
    dfdy[indx[9], indx[6]] = 0
    dfdy[indx[9], indx[7]] = -k[17]*y[:,9]
    dfdy[indx[9], indx[8]] = 0
    dfdy[indx[9], indx[9]] = -4*k[15]*y[:,9] - k[17]*y[:,7] - k[19]*y[:,4]
    dfdy[indx[9], indx[10]] = 2*k[16]*y[:,1]**3
    dfdy[indx[9], indx[11]] = k[18]*y[:,1]**3 + k[20]*y[:,2]
    dfdy[indx[10], indx[0]] = 0
    dfdy[indx[10], indx[1]] = -3*k[16]*y[:,10]*y[:,1]**2
    dfdy[indx[10], indx[2]] = 0
    dfdy[indx[10], indx[3]] = 0
    dfdy[indx[10], indx[4]] = 0
    dfdy[indx[10], indx[5]] = 0
    dfdy[indx[10], indx[6]] = 0
    dfdy[indx[10], indx[7]] = 0
    dfdy[indx[10], indx[8]] = 0
    dfdy[indx[10], indx[9]] = 2*k[15]*y[:,9]
    dfdy[indx[10], indx[10]] = -k[16]*y[:,1]**3
    dfdy[indx[10], indx[11]] = 0
    dfdy[indx[11], indx[0]] = 0
    dfdy[indx[11], indx[1]] = -3*k[18]*y[:,11]*y[:,1]**2
    dfdy[indx[11], indx[2]] = -k[20]*y[:,11]
    dfdy[indx[11], indx[3]] = 0
    dfdy[indx[11], indx[4]] = k[19]*y[:,9]
    dfdy[indx[11], indx[5]] = 0
    dfdy[indx[11], indx[6]] = 0
    dfdy[indx[11], indx[7]] = k[17]*y[:,9]
    dfdy[indx[11], indx[8]] = 0
    dfdy[indx[11], indx[9]] = k[17]*y[:,7] + k[19]*y[:,4]
    dfdy[indx[11], indx[10]] = 0
    dfdy[indx[11], indx[11]] = -k[18]*y[:,1]**3 - k[20]*y[:,2]
    return dfdy


def neg_symjac(y, M, k):
    nz = vulcan_cfg.nz
    dfdy = np.zeros(shape=[ni*nz, ni*nz])
    indx = []
    for j in range(ni):
        indx.append( np.arange(j,j+ni*nz,ni) )
    dfdy[indx[0], indx[0]] = -(-k[1]*y[:,1] - k[3]*y[:,4] - k[6]*y[:,3])
    dfdy[indx[0], indx[1]] = -(-k[1]*y[:,0] + k[5]*y[:,6])
    dfdy[indx[0], indx[2]] = -(k[2]*y[:,3])
    dfdy[indx[0], indx[3]] = -(k[2]*y[:,2] + k[4]*y[:,5] - k[6]*y[:,0])
    dfdy[indx[0], indx[4]] = -(-k[3]*y[:,0])
    dfdy[indx[0], indx[5]] = -(k[4]*y[:,3])
    dfdy[indx[0], indx[6]] = -(k[5]*y[:,1])
    dfdy[indx[0], indx[7]] = -(0)
    dfdy[indx[0], indx[8]] = -(0)
    dfdy[indx[0], indx[9]] = -(0)
    dfdy[indx[0], indx[10]] = -(0)
    dfdy[indx[0], indx[11]] = -(0)
    dfdy[indx[1], indx[0]] = -(-k[1]*y[:,1] + k[6]*y[:,3])
    dfdy[indx[1], indx[1]] = -(-M*k[8] - 9*k[10]*y[:,1]**2*y[:,4] - 9*k[12]*y[:,1]**2*y[:,8] - 9*k[16]*y[:,10]*y[:,1]**2 - 9*k[18]*y[:,11]*y[:,1]**2 - k[1]*y[:,0] - k[5]*y[:,6])
    dfdy[indx[1], indx[2]] = -(k[2]*y[:,3] + 3*k[9]*y[:,7])
    dfdy[indx[1], indx[3]] = -(2*M*k[7]*y[:,3] + k[2]*y[:,2] + k[6]*y[:,0])
    dfdy[indx[1], indx[4]] = -(-3*k[10]*y[:,1]**3)
    dfdy[indx[1], indx[5]] = -(0)
    dfdy[indx[1], indx[6]] = -(-k[5]*y[:,1])
    dfdy[indx[1], indx[7]] = -(6*k[11]*y[:,7] + 3*k[17]*y[:,9] + 3*k[9]*y[:,2])
    dfdy[indx[1], indx[8]] = -(-3*k[12]*y[:,1]**3)
    dfdy[indx[1], indx[9]] = -(6*k[15]*y[:,9] + 3*k[17]*y[:,7])
    dfdy[indx[1], indx[10]] = -(-3*k[16]*y[:,1]**3)
    dfdy[indx[1], indx[11]] = -(-3*k[18]*y[:,1]**3)
    dfdy[indx[2], indx[0]] = -(k[1]*y[:,1])
    dfdy[indx[2], indx[1]] = -(3*k[10]*y[:,1]**2*y[:,4] + k[1]*y[:,0])
    dfdy[indx[2], indx[2]] = -(-k[14]*y[:,8] - k[20]*y[:,11] - k[2]*y[:,3] - k[9]*y[:,7])
    dfdy[indx[2], indx[3]] = -(-k[2]*y[:,2])
    dfdy[indx[2], indx[4]] = -(k[10]*y[:,1]**3 + k[13]*y[:,7] + k[19]*y[:,9])
    dfdy[indx[2], indx[5]] = -(0)
    dfdy[indx[2], indx[6]] = -(0)
    dfdy[indx[2], indx[7]] = -(k[13]*y[:,4] - k[9]*y[:,2])
    dfdy[indx[2], indx[8]] = -(-k[14]*y[:,2])
    dfdy[indx[2], indx[9]] = -(k[19]*y[:,4])
    dfdy[indx[2], indx[10]] = -(0)
    dfdy[indx[2], indx[11]] = -(-k[20]*y[:,2])
    dfdy[indx[3], indx[0]] = -(k[1]*y[:,1] + k[3]*y[:,4] - k[6]*y[:,3])
    dfdy[indx[3], indx[1]] = -(2*M*k[8] + k[1]*y[:,0] + k[5]*y[:,6])
    dfdy[indx[3], indx[2]] = -(-k[2]*y[:,3])
    dfdy[indx[3], indx[3]] = -(-4*M*k[7]*y[:,3] - k[2]*y[:,2] - k[4]*y[:,5] - k[6]*y[:,0])
    dfdy[indx[3], indx[4]] = -(k[3]*y[:,0])
    dfdy[indx[3], indx[5]] = -(-k[4]*y[:,3])
    dfdy[indx[3], indx[6]] = -(k[5]*y[:,1])
    dfdy[indx[3], indx[7]] = -(0)
    dfdy[indx[3], indx[8]] = -(0)
    dfdy[indx[3], indx[9]] = -(0)
    dfdy[indx[3], indx[10]] = -(0)
    dfdy[indx[3], indx[11]] = -(0)
    dfdy[indx[4], indx[0]] = -(-k[3]*y[:,4])
    dfdy[indx[4], indx[1]] = -(-3*k[10]*y[:,1]**2*y[:,4])
    dfdy[indx[4], indx[2]] = -(k[14]*y[:,8] + k[20]*y[:,11] + k[9]*y[:,7])
    dfdy[indx[4], indx[3]] = -(k[4]*y[:,5])
    dfdy[indx[4], indx[4]] = -(-k[10]*y[:,1]**3 - k[13]*y[:,7] - k[19]*y[:,9] - k[3]*y[:,0])
    dfdy[indx[4], indx[5]] = -(k[4]*y[:,3])
    dfdy[indx[4], indx[6]] = -(0)
    dfdy[indx[4], indx[7]] = -(-k[13]*y[:,4] + k[9]*y[:,2])
    dfdy[indx[4], indx[8]] = -(k[14]*y[:,2])
    dfdy[indx[4], indx[9]] = -(-k[19]*y[:,4])
    dfdy[indx[4], indx[10]] = -(0)
    dfdy[indx[4], indx[11]] = -(k[20]*y[:,2])
    dfdy[indx[5], indx[0]] = -(k[3]*y[:,4])
    dfdy[indx[5], indx[1]] = -(0)
    dfdy[indx[5], indx[2]] = -(0)
    dfdy[indx[5], indx[3]] = -(-k[4]*y[:,5])
    dfdy[indx[5], indx[4]] = -(k[3]*y[:,0])
    dfdy[indx[5], indx[5]] = -(-k[4]*y[:,3])
    dfdy[indx[5], indx[6]] = -(0)
    dfdy[indx[5], indx[7]] = -(0)
    dfdy[indx[5], indx[8]] = -(0)
    dfdy[indx[5], indx[9]] = -(0)
    dfdy[indx[5], indx[10]] = -(0)
    dfdy[indx[5], indx[11]] = -(0)
    dfdy[indx[6], indx[0]] = -(k[6]*y[:,3])
    dfdy[indx[6], indx[1]] = -(-k[5]*y[:,6])
    dfdy[indx[6], indx[2]] = -(0)
    dfdy[indx[6], indx[3]] = -(k[6]*y[:,0])
    dfdy[indx[6], indx[4]] = -(0)
    dfdy[indx[6], indx[5]] = -(0)
    dfdy[indx[6], indx[6]] = -(-k[5]*y[:,1])
    dfdy[indx[6], indx[7]] = -(0)
    dfdy[indx[6], indx[8]] = -(0)
    dfdy[indx[6], indx[9]] = -(0)
    dfdy[indx[6], indx[10]] = -(0)
    dfdy[indx[6], indx[11]] = -(0)
    dfdy[indx[7], indx[0]] = -(0)
    dfdy[indx[7], indx[1]] = -(3*k[10]*y[:,1]**2*y[:,4] + 6*k[12]*y[:,1]**2*y[:,8] + 3*k[18]*y[:,11]*y[:,1]**2)
    dfdy[indx[7], indx[2]] = -(k[14]*y[:,8] - k[9]*y[:,7])
    dfdy[indx[7], indx[3]] = -(0)
    dfdy[indx[7], indx[4]] = -(k[10]*y[:,1]**3 - k[13]*y[:,7])
    dfdy[indx[7], indx[5]] = -(0)
    dfdy[indx[7], indx[6]] = -(0)
    dfdy[indx[7], indx[7]] = -(-4*k[11]*y[:,7] - k[13]*y[:,4] - k[17]*y[:,9] - k[9]*y[:,2])
    dfdy[indx[7], indx[8]] = -(2*k[12]*y[:,1]**3 + k[14]*y[:,2])
    dfdy[indx[7], indx[9]] = -(-k[17]*y[:,7])
    dfdy[indx[7], indx[10]] = -(0)
    dfdy[indx[7], indx[11]] = -(k[18]*y[:,1]**3)
    dfdy[indx[8], indx[0]] = -(0)
    dfdy[indx[8], indx[1]] = -(-3*k[12]*y[:,1]**2*y[:,8])
    dfdy[indx[8], indx[2]] = -(-k[14]*y[:,8])
    dfdy[indx[8], indx[3]] = -(0)
    dfdy[indx[8], indx[4]] = -(k[13]*y[:,7])
    dfdy[indx[8], indx[5]] = -(0)
    dfdy[indx[8], indx[6]] = -(0)
    dfdy[indx[8], indx[7]] = -(2*k[11]*y[:,7] + k[13]*y[:,4])
    dfdy[indx[8], indx[8]] = -(-k[12]*y[:,1]**3 - k[14]*y[:,2])
    dfdy[indx[8], indx[9]] = -(0)
    dfdy[indx[8], indx[10]] = -(0)
    dfdy[indx[8], indx[11]] = -(0)
    dfdy[indx[9], indx[0]] = -(0)
    dfdy[indx[9], indx[1]] = -(6*k[16]*y[:,10]*y[:,1]**2 + 3*k[18]*y[:,11]*y[:,1]**2)
    dfdy[indx[9], indx[2]] = -(k[20]*y[:,11])
    dfdy[indx[9], indx[3]] = -(0)
    dfdy[indx[9], indx[4]] = -(-k[19]*y[:,9])
    dfdy[indx[9], indx[5]] = -(0)
    dfdy[indx[9], indx[6]] = -(0)
    dfdy[indx[9], indx[7]] = -(-k[17]*y[:,9])
    dfdy[indx[9], indx[8]] = -(0)
    dfdy[indx[9], indx[9]] = -(-4*k[15]*y[:,9] - k[17]*y[:,7] - k[19]*y[:,4])
    dfdy[indx[9], indx[10]] = -(2*k[16]*y[:,1]**3)
    dfdy[indx[9], indx[11]] = -(k[18]*y[:,1]**3 + k[20]*y[:,2])
    dfdy[indx[10], indx[0]] = -(0)
    dfdy[indx[10], indx[1]] = -(-3*k[16]*y[:,10]*y[:,1]**2)
    dfdy[indx[10], indx[2]] = -(0)
    dfdy[indx[10], indx[3]] = -(0)
    dfdy[indx[10], indx[4]] = -(0)
    dfdy[indx[10], indx[5]] = -(0)
    dfdy[indx[10], indx[6]] = -(0)
    dfdy[indx[10], indx[7]] = -(0)
    dfdy[indx[10], indx[8]] = -(0)
    dfdy[indx[10], indx[9]] = -(2*k[15]*y[:,9])
    dfdy[indx[10], indx[10]] = -(-k[16]*y[:,1]**3)
    dfdy[indx[10], indx[11]] = -(0)
    dfdy[indx[11], indx[0]] = -(0)
    dfdy[indx[11], indx[1]] = -(-3*k[18]*y[:,11]*y[:,1]**2)
    dfdy[indx[11], indx[2]] = -(-k[20]*y[:,11])
    dfdy[indx[11], indx[3]] = -(0)
    dfdy[indx[11], indx[4]] = -(k[19]*y[:,9])
    dfdy[indx[11], indx[5]] = -(0)
    dfdy[indx[11], indx[6]] = -(0)
    dfdy[indx[11], indx[7]] = -(k[17]*y[:,9])
    dfdy[indx[11], indx[8]] = -(0)
    dfdy[indx[11], indx[9]] = -(k[17]*y[:,7] + k[19]*y[:,4])
    dfdy[indx[11], indx[10]] = -(0)
    dfdy[indx[11], indx[11]] = -(-k[18]*y[:,1]**3 - k[20]*y[:,2])
    return dfdy
