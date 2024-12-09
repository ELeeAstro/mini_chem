import numpy as np
import matplotlib.pylab as plt

fname = 'dlsode.txt'

nlay = 100
nsp = 13

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


ni = 50
vi = 0

for i in range(nlay):
  plt.plot(time[:],VMR[:,i,vi])

plt.yscale('log')

plt.show()



