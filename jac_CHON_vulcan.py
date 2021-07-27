'''
## Reaction ##

# Chemical Network and Photolysis Reactions

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
CH4 + CO -> C2H2 + H + OH
#S19
NH3 + NH3 -> N2 + H2 + H2 + H2
#S21
CH4 + NH3 -> HCN + H2 + H2 + H2
#S23
CO + NH3 -> HCN + H2O
#S25
C2H2 + H2 + H2 -> C2H6


## Mapping ##

OH: y[0], H2: y[1], H2O: y[2], H: y[3], CO: y[4], CO2: y[5], O: y[6], C: y[7], CH4: y[8], C2H2: y[9], C2H6: y[10], NH3: y[11], N2: y[12], HCN: y[13],

OH	0	 -1*v_1(k, M, y[2], y[1], y[3], y[0]) -1*v_3(k, M, y[4], y[5], y[3], y[0]) +1*v_5(k, M, y[3], y[1], y[6], y[0]) -1*v_7(k, M, y[4], y[7], y[3], y[0]) +1*v_17(k, M, y[9], y[8], y[4], y[3], y[0])
H2	1	 -1*v_1(k, M, y[2], y[1], y[3], y[0]) -1*v_5(k, M, y[3], y[1], y[6], y[0]) +1*v_9(k, M, y[1], y[3]) +1*v_11(k, M, y[4], y[2], y[1], y[8]) +1*v_11(k, M, y[4], y[2], y[1], y[8]) +1*v_11(k, M, y[4], y[2], y[1], y[8]) +1*v_13(k, M, y[9], y[1], y[8]) +1*v_13(k, M, y[9], y[1], y[8]) +1*v_13(k, M, y[9], y[1], y[8]) +1*v_15(k, M, y[10], y[1], y[8]) +1*v_19(k, M, y[1], y[12], y[11]) +1*v_19(k, M, y[1], y[12], y[11]) +1*v_19(k, M, y[1], y[12], y[11]) +1*v_21(k, M, y[13], y[1], y[11], y[8]) +1*v_21(k, M, y[13], y[1], y[11], y[8]) +1*v_21(k, M, y[13], y[1], y[11], y[8]) -1*v_25(k, M, y[10], y[9], y[1]) -1*v_25(k, M, y[10], y[9], y[1])
H2O	2	 +1*v_1(k, M, y[2], y[1], y[3], y[0]) -1*v_11(k, M, y[4], y[2], y[1], y[8]) +1*v_23(k, M, y[4], y[2], y[11], y[13])
H	3	 +1*v_1(k, M, y[2], y[1], y[3], y[0]) +1*v_3(k, M, y[4], y[5], y[3], y[0]) +1*v_5(k, M, y[3], y[1], y[6], y[0]) +1*v_7(k, M, y[4], y[7], y[3], y[0]) -1*v_9(k, M, y[1], y[3]) -1*v_9(k, M, y[1], y[3]) +1*v_17(k, M, y[9], y[8], y[4], y[3], y[0])
CO	4	 -1*v_3(k, M, y[4], y[5], y[3], y[0]) +1*v_7(k, M, y[4], y[7], y[3], y[0]) +1*v_11(k, M, y[4], y[2], y[1], y[8]) -1*v_17(k, M, y[9], y[8], y[4], y[3], y[0]) -1*v_23(k, M, y[4], y[2], y[11], y[13])
CO2	5	 +1*v_3(k, M, y[4], y[5], y[3], y[0])
O	6	 -1*v_5(k, M, y[3], y[1], y[6], y[0])
C	7	 -1*v_7(k, M, y[4], y[7], y[3], y[0])
CH4	8	 -1*v_11(k, M, y[4], y[2], y[1], y[8]) -1*v_13(k, M, y[9], y[1], y[8]) -1*v_13(k, M, y[9], y[1], y[8]) -1*v_15(k, M, y[10], y[1], y[8]) -1*v_15(k, M, y[10], y[1], y[8]) -1*v_17(k, M, y[9], y[8], y[4], y[3], y[0]) -1*v_21(k, M, y[13], y[1], y[11], y[8])
C2H2	9	 +1*v_13(k, M, y[9], y[1], y[8]) +1*v_17(k, M, y[9], y[8], y[4], y[3], y[0]) -1*v_25(k, M, y[10], y[9], y[1])
C2H6	10	 +1*v_15(k, M, y[10], y[1], y[8]) +1*v_25(k, M, y[10], y[9], y[1])
NH3	11	 -1*v_19(k, M, y[1], y[12], y[11]) -1*v_19(k, M, y[1], y[12], y[11]) -1*v_21(k, M, y[13], y[1], y[11], y[8]) -1*v_23(k, M, y[4], y[2], y[11], y[13])
N2	12	 +1*v_19(k, M, y[1], y[12], y[11])
HCN	13	 +1*v_21(k, M, y[13], y[1], y[11], y[8]) +1*v_23(k, M, y[4], y[2], y[11], y[13])
'''

#species list
spec_list = ['OH', 'H2', 'H2O', 'H', 'CO', 'CO2', 'O', 'C', 'CH4', 'C2H2', 'C2H6', 'NH3', 'N2', 'HCN']
# the total number of species
ni = 14
# the total number of reactions (forward and reverse)
nr = 26

# store the products and the reactants in the 1st and the 2nd element for reaction j (without M)
re_dict = {1:[['OH', 'H2'], ['H2O', 'H']], 2:[['H2O', 'H'], ['OH', 'H2']], 3:[['OH', 'CO'], ['H', 'CO2']], 4:[['H', 'CO2'], ['OH', 'CO']], 5:[['O', 'H2'], ['OH', 'H']], 6:[['OH', 'H'], ['O', 'H2']], 7:[['OH', 'C'], ['CO', 'H']], 8:[['CO', 'H'], ['OH', 'C']], 9:[['H', 'H'], ['H2']], 10:[['H2'], ['H', 'H']], 11:[['CH4', 'H2O'], ['CO', 'H2', 'H2', 'H2']], 12:[['CO', 'H2', 'H2', 'H2'], ['CH4', 'H2O']], 13:[['CH4', 'CH4'], ['C2H2', 'H2', 'H2', 'H2']], 14:[['C2H2', 'H2', 'H2', 'H2'], ['CH4', 'CH4']], 15:[['CH4', 'CH4'], ['C2H6', 'H2']], 16:[['C2H6', 'H2'], ['CH4', 'CH4']], 17:[['CH4', 'CO'], ['C2H2', 'H', 'OH']], 18:[['C2H2', 'H', 'OH'], ['CH4', 'CO']], 19:[['NH3', 'NH3'], ['N2', 'H2', 'H2', 'H2']], 20:[['N2', 'H2', 'H2', 'H2'], ['NH3', 'NH3']], 21:[['CH4', 'NH3'], ['HCN', 'H2', 'H2', 'H2']], 22:[['HCN', 'H2', 'H2', 'H2'], ['CH4', 'NH3']], 23:[['CO', 'NH3'], ['HCN', 'H2O']], 24:[['HCN', 'H2O'], ['CO', 'NH3']], 25:[['C2H2', 'H2', 'H2'], ['C2H6']], 26:[['C2H6'], ['C2H2', 'H2', 'H2']]}


# store the products and the reactants in the 1st and the 2nd element for reaction j (with M)
re_wM_dict = {1:[['OH', 'H2'], ['H2O', 'H']], 2:[['H2O', 'H'], ['OH', 'H2']], 3:[['OH', 'CO'], ['H', 'CO2']], 4:[['H', 'CO2'], ['OH', 'CO']], 5:[['O', 'H2'], ['OH', 'H']], 6:[['OH', 'H'], ['O', 'H2']], 7:[['OH', 'C'], ['CO', 'H']], 8:[['CO', 'H'], ['OH', 'C']], 9:[['H', 'H', 'M'], ['H2', 'M']], 10:[['H2', 'M'], ['H', 'H', 'M']], 11:[['CH4', 'H2O'], ['CO', 'H2', 'H2', 'H2']], 12:[['CO', 'H2', 'H2', 'H2'], ['CH4', 'H2O']], 13:[['CH4', 'CH4'], ['C2H2', 'H2', 'H2', 'H2']], 14:[['C2H2', 'H2', 'H2', 'H2'], ['CH4', 'CH4']], 15:[['CH4', 'CH4'], ['C2H6', 'H2']], 16:[['C2H6', 'H2'], ['CH4', 'CH4']], 17:[['CH4', 'CO'], ['C2H2', 'H', 'OH']], 18:[['C2H2', 'H', 'OH'], ['CH4', 'CO']], 19:[['NH3', 'NH3'], ['N2', 'H2', 'H2', 'H2']], 20:[['N2', 'H2', 'H2', 'H2'], ['NH3', 'NH3']], 21:[['CH4', 'NH3'], ['HCN', 'H2', 'H2', 'H2']], 22:[['HCN', 'H2', 'H2', 'H2'], ['CH4', 'NH3']], 23:[['CO', 'NH3'], ['HCN', 'H2O']], 24:[['HCN', 'H2O'], ['CO', 'NH3']], 25:[['C2H2', 'H2', 'H2'], ['C2H6']], 26:[['C2H6'], ['C2H2', 'H2', 'H2']]}


def chemdf(y, M, k):
    y = np.transpose(y)
    dydt = np.zeros(shape=y.shape)
    dydt[0] =  -1*v_1(k, M, y[2], y[1], y[3], y[0]) -1*v_3(k, M, y[4], y[5], y[3], y[0]) +1*v_5(k, M, y[3], y[1], y[6], y[0]) -1*v_7(k, M, y[4], y[7], y[3], y[0]) +1*v_17(k, M, y[9], y[8], y[4], y[3], y[0])
    dydt[1] =  -1*v_1(k, M, y[2], y[1], y[3], y[0]) -1*v_5(k, M, y[3], y[1], y[6], y[0]) +1*v_9(k, M, y[1], y[3]) +1*v_11(k, M, y[4], y[2], y[1], y[8]) +1*v_11(k, M, y[4], y[2], y[1], y[8]) +1*v_11(k, M, y[4], y[2], y[1], y[8]) +1*v_13(k, M, y[9], y[1], y[8]) +1*v_13(k, M, y[9], y[1], y[8]) +1*v_13(k, M, y[9], y[1], y[8]) +1*v_15(k, M, y[10], y[1], y[8]) +1*v_19(k, M, y[1], y[12], y[11]) +1*v_19(k, M, y[1], y[12], y[11]) +1*v_19(k, M, y[1], y[12], y[11]) +1*v_21(k, M, y[13], y[1], y[11], y[8]) +1*v_21(k, M, y[13], y[1], y[11], y[8]) +1*v_21(k, M, y[13], y[1], y[11], y[8]) -1*v_25(k, M, y[10], y[9], y[1]) -1*v_25(k, M, y[10], y[9], y[1])
    dydt[2] =  +1*v_1(k, M, y[2], y[1], y[3], y[0]) -1*v_11(k, M, y[4], y[2], y[1], y[8]) +1*v_23(k, M, y[4], y[2], y[11], y[13])
    dydt[3] =  +1*v_1(k, M, y[2], y[1], y[3], y[0]) +1*v_3(k, M, y[4], y[5], y[3], y[0]) +1*v_5(k, M, y[3], y[1], y[6], y[0]) +1*v_7(k, M, y[4], y[7], y[3], y[0]) -1*v_9(k, M, y[1], y[3]) -1*v_9(k, M, y[1], y[3]) +1*v_17(k, M, y[9], y[8], y[4], y[3], y[0])
    dydt[4] =  -1*v_3(k, M, y[4], y[5], y[3], y[0]) +1*v_7(k, M, y[4], y[7], y[3], y[0]) +1*v_11(k, M, y[4], y[2], y[1], y[8]) -1*v_17(k, M, y[9], y[8], y[4], y[3], y[0]) -1*v_23(k, M, y[4], y[2], y[11], y[13])
    dydt[5] =  +1*v_3(k, M, y[4], y[5], y[3], y[0])
    dydt[6] =  -1*v_5(k, M, y[3], y[1], y[6], y[0])
    dydt[7] =  -1*v_7(k, M, y[4], y[7], y[3], y[0])
    dydt[8] =  -1*v_11(k, M, y[4], y[2], y[1], y[8]) -1*v_13(k, M, y[9], y[1], y[8]) -1*v_13(k, M, y[9], y[1], y[8]) -1*v_15(k, M, y[10], y[1], y[8]) -1*v_15(k, M, y[10], y[1], y[8]) -1*v_17(k, M, y[9], y[8], y[4], y[3], y[0]) -1*v_21(k, M, y[13], y[1], y[11], y[8])
    dydt[9] =  +1*v_13(k, M, y[9], y[1], y[8]) +1*v_17(k, M, y[9], y[8], y[4], y[3], y[0]) -1*v_25(k, M, y[10], y[9], y[1])
    dydt[10] =  +1*v_15(k, M, y[10], y[1], y[8]) +1*v_25(k, M, y[10], y[9], y[1])
    dydt[11] =  -1*v_19(k, M, y[1], y[12], y[11]) -1*v_19(k, M, y[1], y[12], y[11]) -1*v_21(k, M, y[13], y[1], y[11], y[8]) -1*v_23(k, M, y[4], y[2], y[11], y[13])
    dydt[12] =  +1*v_19(k, M, y[1], y[12], y[11])
    dydt[13] =  +1*v_21(k, M, y[13], y[1], y[11], y[8]) +1*v_23(k, M, y[4], y[2], y[11], y[13])
    dydt = np.transpose(dydt)
    return dydt

def df(y, M, k):
    df_list = []
    df_list.append(  -(k[1]*y[0]*y[1] - k[2]*y[2]*y[3]) -(k[3]*y[0]*y[4] - k[4]*y[3]*y[5]) +(k[5]*y[6]*y[1] - k[6]*y[0]*y[3]) -(k[7]*y[0]*y[7] - k[8]*y[4]*y[3]) +(k[17]*y[8]*y[4] - k[18]*y[9]*y[3]*y[0]) )
    df_list.append(  -(k[1]*y[0]*y[1] - k[2]*y[2]*y[3]) -(k[5]*y[6]*y[1] - k[6]*y[0]*y[3]) +(k[9]*y[3]*y[3]*M - k[10]*y[1]*M) +(k[11]*y[8]*y[2] - k[12]*y[4]*y[1]*y[1]*y[1]) +(k[11]*y[8]*y[2] - k[12]*y[4]*y[1]*y[1]*y[1]) +(k[11]*y[8]*y[2] - k[12]*y[4]*y[1]*y[1]*y[1]) +(k[13]*y[8]*y[8] - k[14]*y[9]*y[1]*y[1]*y[1]) +(k[13]*y[8]*y[8] - k[14]*y[9]*y[1]*y[1]*y[1]) +(k[13]*y[8]*y[8] - k[14]*y[9]*y[1]*y[1]*y[1]) +(k[15]*y[8]*y[8] - k[16]*y[10]*y[1]) +(k[19]*y[11]*y[11] - k[20]*y[12]*y[1]*y[1]*y[1]) +(k[19]*y[11]*y[11] - k[20]*y[12]*y[1]*y[1]*y[1]) +(k[19]*y[11]*y[11] - k[20]*y[12]*y[1]*y[1]*y[1]) +(k[21]*y[8]*y[11] - k[22]*y[13]*y[1]*y[1]*y[1]) +(k[21]*y[8]*y[11] - k[22]*y[13]*y[1]*y[1]*y[1]) +(k[21]*y[8]*y[11] - k[22]*y[13]*y[1]*y[1]*y[1]) -(k[25]*y[9]*y[1]*y[1] - k[26]*y[10]) -(k[25]*y[9]*y[1]*y[1] - k[26]*y[10]) )
    df_list.append(  +(k[1]*y[0]*y[1] - k[2]*y[2]*y[3]) -(k[11]*y[8]*y[2] - k[12]*y[4]*y[1]*y[1]*y[1]) +(k[23]*y[4]*y[11] - k[24]*y[13]*y[2]) )
    df_list.append(  +(k[1]*y[0]*y[1] - k[2]*y[2]*y[3]) +(k[3]*y[0]*y[4] - k[4]*y[3]*y[5]) +(k[5]*y[6]*y[1] - k[6]*y[0]*y[3]) +(k[7]*y[0]*y[7] - k[8]*y[4]*y[3]) -(k[9]*y[3]*y[3]*M - k[10]*y[1]*M) -(k[9]*y[3]*y[3]*M - k[10]*y[1]*M) +(k[17]*y[8]*y[4] - k[18]*y[9]*y[3]*y[0]) )
    df_list.append(  -(k[3]*y[0]*y[4] - k[4]*y[3]*y[5]) +(k[7]*y[0]*y[7] - k[8]*y[4]*y[3]) +(k[11]*y[8]*y[2] - k[12]*y[4]*y[1]*y[1]*y[1]) -(k[17]*y[8]*y[4] - k[18]*y[9]*y[3]*y[0]) -(k[23]*y[4]*y[11] - k[24]*y[13]*y[2]) )
    df_list.append(  +(k[3]*y[0]*y[4] - k[4]*y[3]*y[5]) )
    df_list.append(  -(k[5]*y[6]*y[1] - k[6]*y[0]*y[3]) )
    df_list.append(  -(k[7]*y[0]*y[7] - k[8]*y[4]*y[3]) )
    df_list.append(  -(k[11]*y[8]*y[2] - k[12]*y[4]*y[1]*y[1]*y[1]) -(k[13]*y[8]*y[8] - k[14]*y[9]*y[1]*y[1]*y[1]) -(k[13]*y[8]*y[8] - k[14]*y[9]*y[1]*y[1]*y[1]) -(k[15]*y[8]*y[8] - k[16]*y[10]*y[1]) -(k[15]*y[8]*y[8] - k[16]*y[10]*y[1]) -(k[17]*y[8]*y[4] - k[18]*y[9]*y[3]*y[0]) -(k[21]*y[8]*y[11] - k[22]*y[13]*y[1]*y[1]*y[1]) )
    df_list.append(  +(k[13]*y[8]*y[8] - k[14]*y[9]*y[1]*y[1]*y[1]) +(k[17]*y[8]*y[4] - k[18]*y[9]*y[3]*y[0]) -(k[25]*y[9]*y[1]*y[1] - k[26]*y[10]) )
    df_list.append(  +(k[15]*y[8]*y[8] - k[16]*y[10]*y[1]) +(k[25]*y[9]*y[1]*y[1] - k[26]*y[10]) )
    df_list.append(  -(k[19]*y[11]*y[11] - k[20]*y[12]*y[1]*y[1]*y[1]) -(k[19]*y[11]*y[11] - k[20]*y[12]*y[1]*y[1]*y[1]) -(k[21]*y[8]*y[11] - k[22]*y[13]*y[1]*y[1]*y[1]) -(k[23]*y[4]*y[11] - k[24]*y[13]*y[2]) )
    df_list.append(  +(k[19]*y[11]*y[11] - k[20]*y[12]*y[1]*y[1]*y[1]) )
    df_list.append(  +(k[21]*y[8]*y[11] - k[22]*y[13]*y[1]*y[1]*y[1]) +(k[23]*y[4]*y[11] - k[24]*y[13]*y[2]) )
    return df_list

#OH + H2 -> H2O + H
v_1 = lambda k, M, H2O, H2, H, OH : k[1]*OH*H2 - k[2]*H2O*H


#OH + CO -> H + CO2
v_3 = lambda k, M, CO, CO2, H, OH : k[3]*OH*CO - k[4]*H*CO2


#O + H2 -> OH + H
v_5 = lambda k, M, H, H2, O, OH : k[5]*O*H2 - k[6]*OH*H


#OH + C -> CO + H
v_7 = lambda k, M, CO, C, H, OH : k[7]*OH*C - k[8]*CO*H


#H + H + M -> H2 + M
v_9 = lambda k, M, H2, H : k[9]*H*H*M - k[10]*H2*M


#CH4 + H2O -> CO + H2 + H2 + H2
v_11 = lambda k, M, CO, H2O, H2, CH4 : k[11]*CH4*H2O - k[12]*CO*H2*H2*H2


#CH4 + CH4 -> C2H2 + H2 + H2 + H2
v_13 = lambda k, M, C2H2, H2, CH4 : k[13]*CH4*CH4 - k[14]*C2H2*H2*H2*H2


#CH4 + CH4 -> C2H6 + H2
v_15 = lambda k, M, C2H6, H2, CH4 : k[15]*CH4*CH4 - k[16]*C2H6*H2


#CH4 + CO -> C2H2 + H + OH
v_17 = lambda k, M, C2H2, CH4, CO, H, OH : k[17]*CH4*CO - k[18]*C2H2*H*OH


#NH3 + NH3 -> N2 + H2 + H2 + H2
v_19 = lambda k, M, H2, N2, NH3 : k[19]*NH3*NH3 - k[20]*N2*H2*H2*H2


#CH4 + NH3 -> HCN + H2 + H2 + H2
v_21 = lambda k, M, HCN, H2, NH3, CH4 : k[21]*CH4*NH3 - k[22]*HCN*H2*H2*H2


#CO + NH3 -> HCN + H2O
v_23 = lambda k, M, CO, H2O, NH3, HCN : k[23]*CO*NH3 - k[24]*HCN*H2O


#C2H2 + H2 + H2 -> C2H6
v_25 = lambda k, M, C2H6, C2H2, H2 : k[25]*C2H2*H2*H2 - k[26]*C2H6


def rate_ans(sp):
    rate_str = {}
    re_sp_dic = {}
    rate_str["OH"] = [ -1*k[1]*y[0]*y[1], +1* k[2]*y[2]*y[3], -1*k[3]*y[0]*y[4], +1* k[4]*y[3]*y[5], +1*k[5]*y[6]*y[1], -1* k[6]*y[0]*y[3], -1*k[7]*y[0]*y[7], +1* k[8]*y[4]*y[3], +1*k[17]*y[8]*y[4], -1* k[18]*y[9]*y[3]*y[0]]
    rate_str["H2"] = [ -1*k[1]*y[0]*y[1], +1* k[2]*y[2]*y[3], -1*k[5]*y[6]*y[1], +1* k[6]*y[0]*y[3], +1*k[9]*y[3]*y[3]*M, -1* k[10]*y[1]*M, +1*k[11]*y[8]*y[2], -1* k[12]*y[4]*y[1]*y[1]*y[1], +1*k[11]*y[8]*y[2], -1* k[12]*y[4]*y[1]*y[1]*y[1], +1*k[11]*y[8]*y[2], -1* k[12]*y[4]*y[1]*y[1]*y[1], +1*k[13]*y[8]*y[8], -1* k[14]*y[9]*y[1]*y[1]*y[1], +1*k[13]*y[8]*y[8], -1* k[14]*y[9]*y[1]*y[1]*y[1], +1*k[13]*y[8]*y[8], -1* k[14]*y[9]*y[1]*y[1]*y[1], +1*k[15]*y[8]*y[8], -1* k[16]*y[10]*y[1], +1*k[19]*y[11]*y[11], -1* k[20]*y[12]*y[1]*y[1]*y[1], +1*k[19]*y[11]*y[11], -1* k[20]*y[12]*y[1]*y[1]*y[1], +1*k[19]*y[11]*y[11], -1* k[20]*y[12]*y[1]*y[1]*y[1], +1*k[21]*y[8]*y[11], -1* k[22]*y[13]*y[1]*y[1]*y[1], +1*k[21]*y[8]*y[11], -1* k[22]*y[13]*y[1]*y[1]*y[1], +1*k[21]*y[8]*y[11], -1* k[22]*y[13]*y[1]*y[1]*y[1], -1*k[25]*y[9]*y[1]*y[1], +1* k[26]*y[10], -1*k[25]*y[9]*y[1]*y[1], +1* k[26]*y[10]]
    rate_str["H2O"] = [ +1*k[1]*y[0]*y[1], -1* k[2]*y[2]*y[3], -1*k[11]*y[8]*y[2], +1* k[12]*y[4]*y[1]*y[1]*y[1], +1*k[23]*y[4]*y[11], -1* k[24]*y[13]*y[2]]
    rate_str["H"] = [ +1*k[1]*y[0]*y[1], -1* k[2]*y[2]*y[3], +1*k[3]*y[0]*y[4], -1* k[4]*y[3]*y[5], +1*k[5]*y[6]*y[1], -1* k[6]*y[0]*y[3], +1*k[7]*y[0]*y[7], -1* k[8]*y[4]*y[3], -1*k[9]*y[3]*y[3]*M, +1* k[10]*y[1]*M, -1*k[9]*y[3]*y[3]*M, +1* k[10]*y[1]*M, +1*k[17]*y[8]*y[4], -1* k[18]*y[9]*y[3]*y[0]]
    rate_str["CO"] = [ -1*k[3]*y[0]*y[4], +1* k[4]*y[3]*y[5], +1*k[7]*y[0]*y[7], -1* k[8]*y[4]*y[3], +1*k[11]*y[8]*y[2], -1* k[12]*y[4]*y[1]*y[1]*y[1], -1*k[17]*y[8]*y[4], +1* k[18]*y[9]*y[3]*y[0], -1*k[23]*y[4]*y[11], +1* k[24]*y[13]*y[2]]
    rate_str["CO2"] = [ +1*k[3]*y[0]*y[4], -1* k[4]*y[3]*y[5]]
    rate_str["O"] = [ -1*k[5]*y[6]*y[1], +1* k[6]*y[0]*y[3]]
    rate_str["C"] = [ -1*k[7]*y[0]*y[7], +1* k[8]*y[4]*y[3]]
    rate_str["CH4"] = [ -1*k[11]*y[8]*y[2], +1* k[12]*y[4]*y[1]*y[1]*y[1], -1*k[13]*y[8]*y[8], +1* k[14]*y[9]*y[1]*y[1]*y[1], -1*k[13]*y[8]*y[8], +1* k[14]*y[9]*y[1]*y[1]*y[1], -1*k[15]*y[8]*y[8], +1* k[16]*y[10]*y[1], -1*k[15]*y[8]*y[8], +1* k[16]*y[10]*y[1], -1*k[17]*y[8]*y[4], +1* k[18]*y[9]*y[3]*y[0], -1*k[21]*y[8]*y[11], +1* k[22]*y[13]*y[1]*y[1]*y[1]]
    rate_str["C2H2"] = [ +1*k[13]*y[8]*y[8], -1* k[14]*y[9]*y[1]*y[1]*y[1], +1*k[17]*y[8]*y[4], -1* k[18]*y[9]*y[3]*y[0], -1*k[25]*y[9]*y[1]*y[1], +1* k[26]*y[10]]
    rate_str["C2H6"] = [ +1*k[15]*y[8]*y[8], -1* k[16]*y[10]*y[1], +1*k[25]*y[9]*y[1]*y[1], -1* k[26]*y[10]]
    rate_str["NH3"] = [ -1*k[19]*y[11]*y[11], +1* k[20]*y[12]*y[1]*y[1]*y[1], -1*k[19]*y[11]*y[11], +1* k[20]*y[12]*y[1]*y[1]*y[1], -1*k[21]*y[8]*y[11], +1* k[22]*y[13]*y[1]*y[1]*y[1], -1*k[23]*y[4]*y[11], +1* k[24]*y[13]*y[2]]
    rate_str["N2"] = [ +1*k[19]*y[11]*y[11], -1* k[20]*y[12]*y[1]*y[1]*y[1]]
    rate_str["HCN"] = [ +1*k[21]*y[8]*y[11], -1* k[22]*y[13]*y[1]*y[1]*y[1], +1*k[23]*y[4]*y[11], -1* k[24]*y[13]*y[2]]
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
    G[7] = lambda T: np.exp( -(-1*gibbs_sp('OH',T)-1*gibbs_sp('C',T)+1*gibbs_sp('CO',T)+1*gibbs_sp('H',T) ) )
    G[9] = lambda T: np.exp( -(-1*gibbs_sp('H',T)-1*gibbs_sp('H',T)+1*gibbs_sp('H2',T) ) )*(corr*T)**1
    G[11] = lambda T: np.exp( -(-1*gibbs_sp('CH4',T)-1*gibbs_sp('H2O',T)+1*gibbs_sp('CO',T)+1*gibbs_sp('H2',T)+1*gibbs_sp('H2',T)+1*gibbs_sp('H2',T) ) )*(corr*T)**-2
    G[13] = lambda T: np.exp( -(-1*gibbs_sp('CH4',T)-1*gibbs_sp('CH4',T)+1*gibbs_sp('C2H2',T)+1*gibbs_sp('H2',T)+1*gibbs_sp('H2',T)+1*gibbs_sp('H2',T) ) )*(corr*T)**-2
    G[15] = lambda T: np.exp( -(-1*gibbs_sp('CH4',T)-1*gibbs_sp('CH4',T)+1*gibbs_sp('C2H6',T)+1*gibbs_sp('H2',T) ) )
    G[17] = lambda T: np.exp( -(-1*gibbs_sp('CH4',T)-1*gibbs_sp('CO',T)+1*gibbs_sp('C2H2',T)+1*gibbs_sp('H',T)+1*gibbs_sp('OH',T) ) )*(corr*T)**-1
    G[19] = lambda T: np.exp( -(-1*gibbs_sp('NH3',T)-1*gibbs_sp('NH3',T)+1*gibbs_sp('N2',T)+1*gibbs_sp('H2',T)+1*gibbs_sp('H2',T)+1*gibbs_sp('H2',T) ) )*(corr*T)**-2
    G[21] = lambda T: np.exp( -(-1*gibbs_sp('CH4',T)-1*gibbs_sp('NH3',T)+1*gibbs_sp('HCN',T)+1*gibbs_sp('H2',T)+1*gibbs_sp('H2',T)+1*gibbs_sp('H2',T) ) )*(corr*T)**-2
    G[23] = lambda T: np.exp( -(-1*gibbs_sp('CO',T)-1*gibbs_sp('NH3',T)+1*gibbs_sp('HCN',T)+1*gibbs_sp('H2O',T) ) )
    G[25] = lambda T: np.exp( -(-1*gibbs_sp('C2H2',T)-1*gibbs_sp('H2',T)-1*gibbs_sp('H2',T)+1*gibbs_sp('C2H6',T) ) )*(corr*T)**2
    return G[i](T)


def symjac(y, M, k):
    nz = vulcan_cfg.nz
    dfdy = np.zeros(shape=[ni*nz, ni*nz])
    indx = []
    for j in range(ni):
        indx.append( np.arange(j,j+ni*nz,ni) )
    dfdy[indx[0], indx[0]] = -k[18]*y[:,3]*y[:,9] - k[1]*y[:,1] - k[3]*y[:,4] - k[6]*y[:,3] - k[7]*y[:,7]
    dfdy[indx[0], indx[1]] = -k[1]*y[:,0] + k[5]*y[:,6]
    dfdy[indx[0], indx[2]] = k[2]*y[:,3]
    dfdy[indx[0], indx[3]] = -k[18]*y[:,0]*y[:,9] + k[2]*y[:,2] + k[4]*y[:,5] - k[6]*y[:,0] + k[8]*y[:,4]
    dfdy[indx[0], indx[4]] = k[17]*y[:,8] - k[3]*y[:,0] + k[8]*y[:,3]
    dfdy[indx[0], indx[5]] = k[4]*y[:,3]
    dfdy[indx[0], indx[6]] = k[5]*y[:,1]
    dfdy[indx[0], indx[7]] = -k[7]*y[:,0]
    dfdy[indx[0], indx[8]] = k[17]*y[:,4]
    dfdy[indx[0], indx[9]] = -k[18]*y[:,0]*y[:,3]
    dfdy[indx[0], indx[10]] = 0
    dfdy[indx[0], indx[11]] = 0
    dfdy[indx[0], indx[12]] = 0
    dfdy[indx[0], indx[13]] = 0
    dfdy[indx[1], indx[0]] = -k[1]*y[:,1] + k[6]*y[:,3]
    dfdy[indx[1], indx[1]] = -M*k[10] - 9*k[12]*y[:,1]**2*y[:,4] - 9*k[14]*y[:,1]**2*y[:,9] - k[16]*y[:,10] - k[1]*y[:,0] - 9*k[20]*y[:,12]*y[:,1]**2 - 9*k[22]*y[:,13]*y[:,1]**2 - 4*k[25]*y[:,1]*y[:,9] - k[5]*y[:,6]
    dfdy[indx[1], indx[2]] = 3*k[11]*y[:,8] + k[2]*y[:,3]
    dfdy[indx[1], indx[3]] = 2*M*k[9]*y[:,3] + k[2]*y[:,2] + k[6]*y[:,0]
    dfdy[indx[1], indx[4]] = -3*k[12]*y[:,1]**3
    dfdy[indx[1], indx[5]] = 0
    dfdy[indx[1], indx[6]] = -k[5]*y[:,1]
    dfdy[indx[1], indx[7]] = 0
    dfdy[indx[1], indx[8]] = 3*k[11]*y[:,2] + 6*k[13]*y[:,8] + 2*k[15]*y[:,8] + 3*k[21]*y[:,11]
    dfdy[indx[1], indx[9]] = -3*k[14]*y[:,1]**3 - 2*k[25]*y[:,1]**2
    dfdy[indx[1], indx[10]] = -k[16]*y[:,1] + 2*k[26]
    dfdy[indx[1], indx[11]] = 6*k[19]*y[:,11] + 3*k[21]*y[:,8]
    dfdy[indx[1], indx[12]] = -3*k[20]*y[:,1]**3
    dfdy[indx[1], indx[13]] = -3*k[22]*y[:,1]**3
    dfdy[indx[2], indx[0]] = k[1]*y[:,1]
    dfdy[indx[2], indx[1]] = 3*k[12]*y[:,1]**2*y[:,4] + k[1]*y[:,0]
    dfdy[indx[2], indx[2]] = -k[11]*y[:,8] - k[24]*y[:,13] - k[2]*y[:,3]
    dfdy[indx[2], indx[3]] = -k[2]*y[:,2]
    dfdy[indx[2], indx[4]] = k[12]*y[:,1]**3 + k[23]*y[:,11]
    dfdy[indx[2], indx[5]] = 0
    dfdy[indx[2], indx[6]] = 0
    dfdy[indx[2], indx[7]] = 0
    dfdy[indx[2], indx[8]] = -k[11]*y[:,2]
    dfdy[indx[2], indx[9]] = 0
    dfdy[indx[2], indx[10]] = 0
    dfdy[indx[2], indx[11]] = k[23]*y[:,4]
    dfdy[indx[2], indx[12]] = 0
    dfdy[indx[2], indx[13]] = -k[24]*y[:,2]
    dfdy[indx[3], indx[0]] = -k[18]*y[:,3]*y[:,9] + k[1]*y[:,1] + k[3]*y[:,4] - k[6]*y[:,3] + k[7]*y[:,7]
    dfdy[indx[3], indx[1]] = 2*M*k[10] + k[1]*y[:,0] + k[5]*y[:,6]
    dfdy[indx[3], indx[2]] = -k[2]*y[:,3]
    dfdy[indx[3], indx[3]] = -4*M*k[9]*y[:,3] - k[18]*y[:,0]*y[:,9] - k[2]*y[:,2] - k[4]*y[:,5] - k[6]*y[:,0] - k[8]*y[:,4]
    dfdy[indx[3], indx[4]] = k[17]*y[:,8] + k[3]*y[:,0] - k[8]*y[:,3]
    dfdy[indx[3], indx[5]] = -k[4]*y[:,3]
    dfdy[indx[3], indx[6]] = k[5]*y[:,1]
    dfdy[indx[3], indx[7]] = k[7]*y[:,0]
    dfdy[indx[3], indx[8]] = k[17]*y[:,4]
    dfdy[indx[3], indx[9]] = -k[18]*y[:,0]*y[:,3]
    dfdy[indx[3], indx[10]] = 0
    dfdy[indx[3], indx[11]] = 0
    dfdy[indx[3], indx[12]] = 0
    dfdy[indx[3], indx[13]] = 0
    dfdy[indx[4], indx[0]] = k[18]*y[:,3]*y[:,9] - k[3]*y[:,4] + k[7]*y[:,7]
    dfdy[indx[4], indx[1]] = -3*k[12]*y[:,1]**2*y[:,4]
    dfdy[indx[4], indx[2]] = k[11]*y[:,8] + k[24]*y[:,13]
    dfdy[indx[4], indx[3]] = k[18]*y[:,0]*y[:,9] + k[4]*y[:,5] - k[8]*y[:,4]
    dfdy[indx[4], indx[4]] = -k[12]*y[:,1]**3 - k[17]*y[:,8] - k[23]*y[:,11] - k[3]*y[:,0] - k[8]*y[:,3]
    dfdy[indx[4], indx[5]] = k[4]*y[:,3]
    dfdy[indx[4], indx[6]] = 0
    dfdy[indx[4], indx[7]] = k[7]*y[:,0]
    dfdy[indx[4], indx[8]] = k[11]*y[:,2] - k[17]*y[:,4]
    dfdy[indx[4], indx[9]] = k[18]*y[:,0]*y[:,3]
    dfdy[indx[4], indx[10]] = 0
    dfdy[indx[4], indx[11]] = -k[23]*y[:,4]
    dfdy[indx[4], indx[12]] = 0
    dfdy[indx[4], indx[13]] = k[24]*y[:,2]
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
    dfdy[indx[5], indx[12]] = 0
    dfdy[indx[5], indx[13]] = 0
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
    dfdy[indx[6], indx[12]] = 0
    dfdy[indx[6], indx[13]] = 0
    dfdy[indx[7], indx[0]] = -k[7]*y[:,7]
    dfdy[indx[7], indx[1]] = 0
    dfdy[indx[7], indx[2]] = 0
    dfdy[indx[7], indx[3]] = k[8]*y[:,4]
    dfdy[indx[7], indx[4]] = k[8]*y[:,3]
    dfdy[indx[7], indx[5]] = 0
    dfdy[indx[7], indx[6]] = 0
    dfdy[indx[7], indx[7]] = -k[7]*y[:,0]
    dfdy[indx[7], indx[8]] = 0
    dfdy[indx[7], indx[9]] = 0
    dfdy[indx[7], indx[10]] = 0
    dfdy[indx[7], indx[11]] = 0
    dfdy[indx[7], indx[12]] = 0
    dfdy[indx[7], indx[13]] = 0
    dfdy[indx[8], indx[0]] = k[18]*y[:,3]*y[:,9]
    dfdy[indx[8], indx[1]] = 3*k[12]*y[:,1]**2*y[:,4] + 6*k[14]*y[:,1]**2*y[:,9] + 2*k[16]*y[:,10] + 3*k[22]*y[:,13]*y[:,1]**2
    dfdy[indx[8], indx[2]] = -k[11]*y[:,8]
    dfdy[indx[8], indx[3]] = k[18]*y[:,0]*y[:,9]
    dfdy[indx[8], indx[4]] = k[12]*y[:,1]**3 - k[17]*y[:,8]
    dfdy[indx[8], indx[5]] = 0
    dfdy[indx[8], indx[6]] = 0
    dfdy[indx[8], indx[7]] = 0
    dfdy[indx[8], indx[8]] = -k[11]*y[:,2] - 4*k[13]*y[:,8] - 4*k[15]*y[:,8] - k[17]*y[:,4] - k[21]*y[:,11]
    dfdy[indx[8], indx[9]] = 2*k[14]*y[:,1]**3 + k[18]*y[:,0]*y[:,3]
    dfdy[indx[8], indx[10]] = 2*k[16]*y[:,1]
    dfdy[indx[8], indx[11]] = -k[21]*y[:,8]
    dfdy[indx[8], indx[12]] = 0
    dfdy[indx[8], indx[13]] = k[22]*y[:,1]**3
    dfdy[indx[9], indx[0]] = -k[18]*y[:,3]*y[:,9]
    dfdy[indx[9], indx[1]] = -3*k[14]*y[:,1]**2*y[:,9] - 2*k[25]*y[:,1]*y[:,9]
    dfdy[indx[9], indx[2]] = 0
    dfdy[indx[9], indx[3]] = -k[18]*y[:,0]*y[:,9]
    dfdy[indx[9], indx[4]] = k[17]*y[:,8]
    dfdy[indx[9], indx[5]] = 0
    dfdy[indx[9], indx[6]] = 0
    dfdy[indx[9], indx[7]] = 0
    dfdy[indx[9], indx[8]] = 2*k[13]*y[:,8] + k[17]*y[:,4]
    dfdy[indx[9], indx[9]] = -k[14]*y[:,1]**3 - k[18]*y[:,0]*y[:,3] - k[25]*y[:,1]**2
    dfdy[indx[9], indx[10]] = k[26]
    dfdy[indx[9], indx[11]] = 0
    dfdy[indx[9], indx[12]] = 0
    dfdy[indx[9], indx[13]] = 0
    dfdy[indx[10], indx[0]] = 0
    dfdy[indx[10], indx[1]] = -k[16]*y[:,10] + 2*k[25]*y[:,1]*y[:,9]
    dfdy[indx[10], indx[2]] = 0
    dfdy[indx[10], indx[3]] = 0
    dfdy[indx[10], indx[4]] = 0
    dfdy[indx[10], indx[5]] = 0
    dfdy[indx[10], indx[6]] = 0
    dfdy[indx[10], indx[7]] = 0
    dfdy[indx[10], indx[8]] = 2*k[15]*y[:,8]
    dfdy[indx[10], indx[9]] = k[25]*y[:,1]**2
    dfdy[indx[10], indx[10]] = -k[16]*y[:,1] - k[26]
    dfdy[indx[10], indx[11]] = 0
    dfdy[indx[10], indx[12]] = 0
    dfdy[indx[10], indx[13]] = 0
    dfdy[indx[11], indx[0]] = 0
    dfdy[indx[11], indx[1]] = 6*k[20]*y[:,12]*y[:,1]**2 + 3*k[22]*y[:,13]*y[:,1]**2
    dfdy[indx[11], indx[2]] = k[24]*y[:,13]
    dfdy[indx[11], indx[3]] = 0
    dfdy[indx[11], indx[4]] = -k[23]*y[:,11]
    dfdy[indx[11], indx[5]] = 0
    dfdy[indx[11], indx[6]] = 0
    dfdy[indx[11], indx[7]] = 0
    dfdy[indx[11], indx[8]] = -k[21]*y[:,11]
    dfdy[indx[11], indx[9]] = 0
    dfdy[indx[11], indx[10]] = 0
    dfdy[indx[11], indx[11]] = -4*k[19]*y[:,11] - k[21]*y[:,8] - k[23]*y[:,4]
    dfdy[indx[11], indx[12]] = 2*k[20]*y[:,1]**3
    dfdy[indx[11], indx[13]] = k[22]*y[:,1]**3 + k[24]*y[:,2]
    dfdy[indx[12], indx[0]] = 0
    dfdy[indx[12], indx[1]] = -3*k[20]*y[:,12]*y[:,1]**2
    dfdy[indx[12], indx[2]] = 0
    dfdy[indx[12], indx[3]] = 0
    dfdy[indx[12], indx[4]] = 0
    dfdy[indx[12], indx[5]] = 0
    dfdy[indx[12], indx[6]] = 0
    dfdy[indx[12], indx[7]] = 0
    dfdy[indx[12], indx[8]] = 0
    dfdy[indx[12], indx[9]] = 0
    dfdy[indx[12], indx[10]] = 0
    dfdy[indx[12], indx[11]] = 2*k[19]*y[:,11]
    dfdy[indx[12], indx[12]] = -k[20]*y[:,1]**3
    dfdy[indx[12], indx[13]] = 0
    dfdy[indx[13], indx[0]] = 0
    dfdy[indx[13], indx[1]] = -3*k[22]*y[:,13]*y[:,1]**2
    dfdy[indx[13], indx[2]] = -k[24]*y[:,13]
    dfdy[indx[13], indx[3]] = 0
    dfdy[indx[13], indx[4]] = k[23]*y[:,11]
    dfdy[indx[13], indx[5]] = 0
    dfdy[indx[13], indx[6]] = 0
    dfdy[indx[13], indx[7]] = 0
    dfdy[indx[13], indx[8]] = k[21]*y[:,11]
    dfdy[indx[13], indx[9]] = 0
    dfdy[indx[13], indx[10]] = 0
    dfdy[indx[13], indx[11]] = k[21]*y[:,8] + k[23]*y[:,4]
    dfdy[indx[13], indx[12]] = 0
    dfdy[indx[13], indx[13]] = -k[22]*y[:,1]**3 - k[24]*y[:,2]
    return dfdy


def neg_symjac(y, M, k):
    nz = vulcan_cfg.nz
    dfdy = np.zeros(shape=[ni*nz, ni*nz])
    indx = []
    for j in range(ni):
        indx.append( np.arange(j,j+ni*nz,ni) )
    dfdy[indx[0], indx[0]] = -(-k[18]*y[:,3]*y[:,9] - k[1]*y[:,1] - k[3]*y[:,4] - k[6]*y[:,3] - k[7]*y[:,7])
    dfdy[indx[0], indx[1]] = -(-k[1]*y[:,0] + k[5]*y[:,6])
    dfdy[indx[0], indx[2]] = -(k[2]*y[:,3])
    dfdy[indx[0], indx[3]] = -(-k[18]*y[:,0]*y[:,9] + k[2]*y[:,2] + k[4]*y[:,5] - k[6]*y[:,0] + k[8]*y[:,4])
    dfdy[indx[0], indx[4]] = -(k[17]*y[:,8] - k[3]*y[:,0] + k[8]*y[:,3])
    dfdy[indx[0], indx[5]] = -(k[4]*y[:,3])
    dfdy[indx[0], indx[6]] = -(k[5]*y[:,1])
    dfdy[indx[0], indx[7]] = -(-k[7]*y[:,0])
    dfdy[indx[0], indx[8]] = -(k[17]*y[:,4])
    dfdy[indx[0], indx[9]] = -(-k[18]*y[:,0]*y[:,3])
    dfdy[indx[0], indx[10]] = -(0)
    dfdy[indx[0], indx[11]] = -(0)
    dfdy[indx[0], indx[12]] = -(0)
    dfdy[indx[0], indx[13]] = -(0)
    dfdy[indx[1], indx[0]] = -(-k[1]*y[:,1] + k[6]*y[:,3])
    dfdy[indx[1], indx[1]] = -(-M*k[10] - 9*k[12]*y[:,1]**2*y[:,4] - 9*k[14]*y[:,1]**2*y[:,9] - k[16]*y[:,10] - k[1]*y[:,0] - 9*k[20]*y[:,12]*y[:,1]**2 - 9*k[22]*y[:,13]*y[:,1]**2 - 4*k[25]*y[:,1]*y[:,9] - k[5]*y[:,6])
    dfdy[indx[1], indx[2]] = -(3*k[11]*y[:,8] + k[2]*y[:,3])
    dfdy[indx[1], indx[3]] = -(2*M*k[9]*y[:,3] + k[2]*y[:,2] + k[6]*y[:,0])
    dfdy[indx[1], indx[4]] = -(-3*k[12]*y[:,1]**3)
    dfdy[indx[1], indx[5]] = -(0)
    dfdy[indx[1], indx[6]] = -(-k[5]*y[:,1])
    dfdy[indx[1], indx[7]] = -(0)
    dfdy[indx[1], indx[8]] = -(3*k[11]*y[:,2] + 6*k[13]*y[:,8] + 2*k[15]*y[:,8] + 3*k[21]*y[:,11])
    dfdy[indx[1], indx[9]] = -(-3*k[14]*y[:,1]**3 - 2*k[25]*y[:,1]**2)
    dfdy[indx[1], indx[10]] = -(-k[16]*y[:,1] + 2*k[26])
    dfdy[indx[1], indx[11]] = -(6*k[19]*y[:,11] + 3*k[21]*y[:,8])
    dfdy[indx[1], indx[12]] = -(-3*k[20]*y[:,1]**3)
    dfdy[indx[1], indx[13]] = -(-3*k[22]*y[:,1]**3)
    dfdy[indx[2], indx[0]] = -(k[1]*y[:,1])
    dfdy[indx[2], indx[1]] = -(3*k[12]*y[:,1]**2*y[:,4] + k[1]*y[:,0])
    dfdy[indx[2], indx[2]] = -(-k[11]*y[:,8] - k[24]*y[:,13] - k[2]*y[:,3])
    dfdy[indx[2], indx[3]] = -(-k[2]*y[:,2])
    dfdy[indx[2], indx[4]] = -(k[12]*y[:,1]**3 + k[23]*y[:,11])
    dfdy[indx[2], indx[5]] = -(0)
    dfdy[indx[2], indx[6]] = -(0)
    dfdy[indx[2], indx[7]] = -(0)
    dfdy[indx[2], indx[8]] = -(-k[11]*y[:,2])
    dfdy[indx[2], indx[9]] = -(0)
    dfdy[indx[2], indx[10]] = -(0)
    dfdy[indx[2], indx[11]] = -(k[23]*y[:,4])
    dfdy[indx[2], indx[12]] = -(0)
    dfdy[indx[2], indx[13]] = -(-k[24]*y[:,2])
    dfdy[indx[3], indx[0]] = -(-k[18]*y[:,3]*y[:,9] + k[1]*y[:,1] + k[3]*y[:,4] - k[6]*y[:,3] + k[7]*y[:,7])
    dfdy[indx[3], indx[1]] = -(2*M*k[10] + k[1]*y[:,0] + k[5]*y[:,6])
    dfdy[indx[3], indx[2]] = -(-k[2]*y[:,3])
    dfdy[indx[3], indx[3]] = -(-4*M*k[9]*y[:,3] - k[18]*y[:,0]*y[:,9] - k[2]*y[:,2] - k[4]*y[:,5] - k[6]*y[:,0] - k[8]*y[:,4])
    dfdy[indx[3], indx[4]] = -(k[17]*y[:,8] + k[3]*y[:,0] - k[8]*y[:,3])
    dfdy[indx[3], indx[5]] = -(-k[4]*y[:,3])
    dfdy[indx[3], indx[6]] = -(k[5]*y[:,1])
    dfdy[indx[3], indx[7]] = -(k[7]*y[:,0])
    dfdy[indx[3], indx[8]] = -(k[17]*y[:,4])
    dfdy[indx[3], indx[9]] = -(-k[18]*y[:,0]*y[:,3])
    dfdy[indx[3], indx[10]] = -(0)
    dfdy[indx[3], indx[11]] = -(0)
    dfdy[indx[3], indx[12]] = -(0)
    dfdy[indx[3], indx[13]] = -(0)
    dfdy[indx[4], indx[0]] = -(k[18]*y[:,3]*y[:,9] - k[3]*y[:,4] + k[7]*y[:,7])
    dfdy[indx[4], indx[1]] = -(-3*k[12]*y[:,1]**2*y[:,4])
    dfdy[indx[4], indx[2]] = -(k[11]*y[:,8] + k[24]*y[:,13])
    dfdy[indx[4], indx[3]] = -(k[18]*y[:,0]*y[:,9] + k[4]*y[:,5] - k[8]*y[:,4])
    dfdy[indx[4], indx[4]] = -(-k[12]*y[:,1]**3 - k[17]*y[:,8] - k[23]*y[:,11] - k[3]*y[:,0] - k[8]*y[:,3])
    dfdy[indx[4], indx[5]] = -(k[4]*y[:,3])
    dfdy[indx[4], indx[6]] = -(0)
    dfdy[indx[4], indx[7]] = -(k[7]*y[:,0])
    dfdy[indx[4], indx[8]] = -(k[11]*y[:,2] - k[17]*y[:,4])
    dfdy[indx[4], indx[9]] = -(k[18]*y[:,0]*y[:,3])
    dfdy[indx[4], indx[10]] = -(0)
    dfdy[indx[4], indx[11]] = -(-k[23]*y[:,4])
    dfdy[indx[4], indx[12]] = -(0)
    dfdy[indx[4], indx[13]] = -(k[24]*y[:,2])
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
    dfdy[indx[5], indx[12]] = -(0)
    dfdy[indx[5], indx[13]] = -(0)
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
    dfdy[indx[6], indx[12]] = -(0)
    dfdy[indx[6], indx[13]] = -(0)
    dfdy[indx[7], indx[0]] = -(-k[7]*y[:,7])
    dfdy[indx[7], indx[1]] = -(0)
    dfdy[indx[7], indx[2]] = -(0)
    dfdy[indx[7], indx[3]] = -(k[8]*y[:,4])
    dfdy[indx[7], indx[4]] = -(k[8]*y[:,3])
    dfdy[indx[7], indx[5]] = -(0)
    dfdy[indx[7], indx[6]] = -(0)
    dfdy[indx[7], indx[7]] = -(-k[7]*y[:,0])
    dfdy[indx[7], indx[8]] = -(0)
    dfdy[indx[7], indx[9]] = -(0)
    dfdy[indx[7], indx[10]] = -(0)
    dfdy[indx[7], indx[11]] = -(0)
    dfdy[indx[7], indx[12]] = -(0)
    dfdy[indx[7], indx[13]] = -(0)
    dfdy[indx[8], indx[0]] = -(k[18]*y[:,3]*y[:,9])
    dfdy[indx[8], indx[1]] = -(3*k[12]*y[:,1]**2*y[:,4] + 6*k[14]*y[:,1]**2*y[:,9] + 2*k[16]*y[:,10] + 3*k[22]*y[:,13]*y[:,1]**2)
    dfdy[indx[8], indx[2]] = -(-k[11]*y[:,8])
    dfdy[indx[8], indx[3]] = -(k[18]*y[:,0]*y[:,9])
    dfdy[indx[8], indx[4]] = -(k[12]*y[:,1]**3 - k[17]*y[:,8])
    dfdy[indx[8], indx[5]] = -(0)
    dfdy[indx[8], indx[6]] = -(0)
    dfdy[indx[8], indx[7]] = -(0)
    dfdy[indx[8], indx[8]] = -(-k[11]*y[:,2] - 4*k[13]*y[:,8] - 4*k[15]*y[:,8] - k[17]*y[:,4] - k[21]*y[:,11])
    dfdy[indx[8], indx[9]] = -(2*k[14]*y[:,1]**3 + k[18]*y[:,0]*y[:,3])
    dfdy[indx[8], indx[10]] = -(2*k[16]*y[:,1])
    dfdy[indx[8], indx[11]] = -(-k[21]*y[:,8])
    dfdy[indx[8], indx[12]] = -(0)
    dfdy[indx[8], indx[13]] = -(k[22]*y[:,1]**3)
    dfdy[indx[9], indx[0]] = -(-k[18]*y[:,3]*y[:,9])
    dfdy[indx[9], indx[1]] = -(-3*k[14]*y[:,1]**2*y[:,9] - 2*k[25]*y[:,1]*y[:,9])
    dfdy[indx[9], indx[2]] = -(0)
    dfdy[indx[9], indx[3]] = -(-k[18]*y[:,0]*y[:,9])
    dfdy[indx[9], indx[4]] = -(k[17]*y[:,8])
    dfdy[indx[9], indx[5]] = -(0)
    dfdy[indx[9], indx[6]] = -(0)
    dfdy[indx[9], indx[7]] = -(0)
    dfdy[indx[9], indx[8]] = -(2*k[13]*y[:,8] + k[17]*y[:,4])
    dfdy[indx[9], indx[9]] = -(-k[14]*y[:,1]**3 - k[18]*y[:,0]*y[:,3] - k[25]*y[:,1]**2)
    dfdy[indx[9], indx[10]] = -(k[26])
    dfdy[indx[9], indx[11]] = -(0)
    dfdy[indx[9], indx[12]] = -(0)
    dfdy[indx[9], indx[13]] = -(0)
    dfdy[indx[10], indx[0]] = -(0)
    dfdy[indx[10], indx[1]] = -(-k[16]*y[:,10] + 2*k[25]*y[:,1]*y[:,9])
    dfdy[indx[10], indx[2]] = -(0)
    dfdy[indx[10], indx[3]] = -(0)
    dfdy[indx[10], indx[4]] = -(0)
    dfdy[indx[10], indx[5]] = -(0)
    dfdy[indx[10], indx[6]] = -(0)
    dfdy[indx[10], indx[7]] = -(0)
    dfdy[indx[10], indx[8]] = -(2*k[15]*y[:,8])
    dfdy[indx[10], indx[9]] = -(k[25]*y[:,1]**2)
    dfdy[indx[10], indx[10]] = -(-k[16]*y[:,1] - k[26])
    dfdy[indx[10], indx[11]] = -(0)
    dfdy[indx[10], indx[12]] = -(0)
    dfdy[indx[10], indx[13]] = -(0)
    dfdy[indx[11], indx[0]] = -(0)
    dfdy[indx[11], indx[1]] = -(6*k[20]*y[:,12]*y[:,1]**2 + 3*k[22]*y[:,13]*y[:,1]**2)
    dfdy[indx[11], indx[2]] = -(k[24]*y[:,13])
    dfdy[indx[11], indx[3]] = -(0)
    dfdy[indx[11], indx[4]] = -(-k[23]*y[:,11])
    dfdy[indx[11], indx[5]] = -(0)
    dfdy[indx[11], indx[6]] = -(0)
    dfdy[indx[11], indx[7]] = -(0)
    dfdy[indx[11], indx[8]] = -(-k[21]*y[:,11])
    dfdy[indx[11], indx[9]] = -(0)
    dfdy[indx[11], indx[10]] = -(0)
    dfdy[indx[11], indx[11]] = -(-4*k[19]*y[:,11] - k[21]*y[:,8] - k[23]*y[:,4])
    dfdy[indx[11], indx[12]] = -(2*k[20]*y[:,1]**3)
    dfdy[indx[11], indx[13]] = -(k[22]*y[:,1]**3 + k[24]*y[:,2])
    dfdy[indx[12], indx[0]] = -(0)
    dfdy[indx[12], indx[1]] = -(-3*k[20]*y[:,12]*y[:,1]**2)
    dfdy[indx[12], indx[2]] = -(0)
    dfdy[indx[12], indx[3]] = -(0)
    dfdy[indx[12], indx[4]] = -(0)
    dfdy[indx[12], indx[5]] = -(0)
    dfdy[indx[12], indx[6]] = -(0)
    dfdy[indx[12], indx[7]] = -(0)
    dfdy[indx[12], indx[8]] = -(0)
    dfdy[indx[12], indx[9]] = -(0)
    dfdy[indx[12], indx[10]] = -(0)
    dfdy[indx[12], indx[11]] = -(2*k[19]*y[:,11])
    dfdy[indx[12], indx[12]] = -(-k[20]*y[:,1]**3)
    dfdy[indx[12], indx[13]] = -(0)
    dfdy[indx[13], indx[0]] = -(0)
    dfdy[indx[13], indx[1]] = -(-3*k[22]*y[:,13]*y[:,1]**2)
    dfdy[indx[13], indx[2]] = -(-k[24]*y[:,13])
    dfdy[indx[13], indx[3]] = -(0)
    dfdy[indx[13], indx[4]] = -(k[23]*y[:,11])
    dfdy[indx[13], indx[5]] = -(0)
    dfdy[indx[13], indx[6]] = -(0)
    dfdy[indx[13], indx[7]] = -(0)
    dfdy[indx[13], indx[8]] = -(k[21]*y[:,11])
    dfdy[indx[13], indx[9]] = -(0)
    dfdy[indx[13], indx[10]] = -(0)
    dfdy[indx[13], indx[11]] = -(k[21]*y[:,8] + k[23]*y[:,4])
    dfdy[indx[13], indx[12]] = -(0)
    dfdy[indx[13], indx[13]] = -(-k[22]*y[:,1]**3 - k[24]*y[:,2])
