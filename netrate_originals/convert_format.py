import numpy as np

dire = 'Tsai2022_net_rates/solar_lodders2019/'
sol = ''

netrat = ['CH4-C2H2','CH4-CO','CH4-HCN','CO-CH4-C2H2','CO-HCN','NH3-N2']
lnet = len(netrat)

for n in range(lnet):

  fname = dire+'NCHO-'+sol+'solar-'+netrat[n]+'-netrate.txt'
  data = np.loadtxt(fname,skiprows=2,usecols=(0,1,2))

  T = data[:,0]
  P =  data[:,1]
  kf = data[:,2]

  T = np.unique(T)
  P = np.unique(P)

  nT = len(T)
  nP = len(P)
  nf = len(kf)

  print(nT, T)
  print(nP, P)

  fout = '1x_NHCO-solar-'+netrat[n]+'.txt'

  f = open(fout,'w')

  f.write('# For the net reaction'  + '\n')
  f.write(str(nT) + ' ' + str(nP) + ' ' + str(nf) + '\n')
  f.write('T [K]' + '\n')
  f.write(" ".join(str(g) for g in T[:]) + '\n')
  f.write('p [bar]' + '\n')
  f.write(" ".join(str(g) for g in P[:]) + '\n')
  f.write('kf [cm6 s-1]' + '\n')
  for i in range(nf):
    f.write(str(kf[i]) + '\n')
