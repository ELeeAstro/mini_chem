import numpy as np
import matplotlib.pylab as plt
import seaborn as sns

fname = 'dlsode.txt'

nlay = 100

sp = ['OH','H2','H2O','H','CO','CO2','O','CH4','C2H2','NH3','N2','HCN','He']
nsp = len(sp)

data = np.loadtxt(fname,skiprows=1)

nt = int(len(data[:,0])/nlay)

print(nlay, nsp, nt)

time = np.zeros(nt)
VMR = np.zeros((nt,nlay,nsp))

nj = 0
for i in range(nt):
  time[i] = data[nj,1]
  for j in range(nlay):
    VMR[i,j,:] = data[nj,2:]
    nj = nj + 1
  print(i,time[i])

# Read T-p-Kzz profile
fname = 'T_p_Kzz.txt'
data = np.loadtxt(fname)
pl = data[:,1]

col = sns.color_palette("hls", nsp)

for i in range(nsp):
  plt.plot(VMR[0,:,i],pl[:],ls='dashed',c=col[i])
  plt.plot(VMR[-1,:,i],pl[:],ls='solid',label=sp[i],c=col[i])


plt.title(str(time[-1]) + ' s',fontsize=14)
plt.legend() 

plt.yscale('log')
plt.xscale('log')

plt.gca().invert_yaxis()

plt.ylabel('p [bar]',fontsize=14)
plt.xlabel('VMR',fontsize=14)

plt.xlim(1e-20,1.0)

plt.show()



