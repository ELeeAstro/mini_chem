import numpy as np
import matplotlib.pylab as plt
import pickle


filename = '../vulcan_benchmark_data/net-T1500KlogP7.0-NCHO-solar_hot_noHe.vul'

ivmr = 0
with open(filename, 'rb') as handle:
  data = pickle.load(handle)
  time_V = data['variable']['t_time']
  mix_time_V = np.array(data['variable']['y_time'])[:,0,ivmr] /float(data['atm']['n_0'])
  print(mix_time_V)
# sp the name of species, e.g., ‘H2O’


integrator = ['dlsodes']
sp = ['OH','H2','H2O','H','CO','CO2','O','CH4','C2H2','NH3','N2','HCN']

nint = len(integrator)
nsp = len(sp)

fig = plt.figure()

plt.plot(time_V,mix_time_V,label='Vulcan',ls='dotted',lw=3,c='black')

for i in range(nint):
  data = np.loadtxt(integrator[i]+'.txt',skiprows=1)
  time = data[:,1]
  VMR = data[:,2:]
  plt.plot(time,VMR[:,ivmr],label=integrator[i])


plt.legend()

plt.yscale('log')
plt.xscale('log')

plt.show()
