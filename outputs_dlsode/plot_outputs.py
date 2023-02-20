import numpy as np
import matplotlib.pylab as plt

integrator = ['dlsode','../outputs_dvode/dvode']
sp = ['OH','H2','H2O','H','CO','CO2','O','CH4','C2H2','NH3','N2','HCN']

nint = len(integrator)
nsp = len(sp)

fig = plt.figure()

for i in range(nint):
  data = np.loadtxt(integrator[i]+'.txt',skiprows=1)
  time = data[:,1]
  VMR = data[:,2:]
  ivmr = 2
  plt.plot(time,VMR[:,ivmr],label=integrator[i])


plt.legend()

plt.yscale('log')
plt.xscale('log')

plt.show()
