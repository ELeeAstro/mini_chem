import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
import seaborn as sns

# Nice Fonts
plt.rc('text', usetex=False)
plt.rc('font', family='serif')

path = '../netrate_tables/'

fname = 'solar-CH4-C2H2.txt'
title = r'CH$_{4}$ + CH$_{4}$ $\rightarrow$ C$_{2}$H$_{2}$ + 3H$_{2}$'

# fname = path+'solar-CH4-CO.txt'
# title = r'CH$_{4}$ + H$_{2}$O $\rightarrow$ CO + 3H$_{2}$'

# fname = path+'solar-CH4-HCN.txt'
# title = r'CH$_{4}$ + NH$_{3}$ $\rightarrow$ HCN + 3H$_{2}$'

# fname = path+'solar-CO-CH4-C2H2.txt'
# title = r'CO + CH$_{4}$  $\rightarrow$ C$_{2}$H$_{2}$ + H$_{2}$O'

# fname = path+'solar-CO-HCN.txt'
# title = r'CO + NH$_{3}$ $\rightarrow$ HCN + H$_{2}$O'

# fname = path+'solar-NH3-N2.txt'
# title = r'NH$_{3}$ + NH$_{3}$ $\rightarrow$ N$_{2}$ + 3H$_{2}$'

T = np.loadtxt(path+fname,max_rows=1,skiprows=3)
p = np.loadtxt(path+fname,max_rows=1,skiprows=5)

nT = len(T)
nP = len(p)

kf_1D = np.loadtxt(path+fname,skiprows=7)
kf1 = np.zeros((nT,nP))
kf1 = np.log10(kf_1D.reshape((nT,nP)))

kf = np.zeros((nP,nT))
for i in range(nP):
    kf[i,:] = kf1[:,i]

fig = plt.figure()

zmin = np.amin(kf)
zmax = np.amax(kf)
lev = np.linspace(zmin,zmax,20)

cmap = sns.color_palette("rocket", as_cmap=True)
#cmap = 'RdYlBu_r'
fmt = '%d'

CS = plt.contourf(T,p,kf,levels=lev,extend='both',cmap=cmap)
# Colour bar and formatting
for c in CS.collections:
    c.set_edgecolor("face")
CB = plt.colorbar(CS, format=fmt)
CB.solids.set_rasterized(True)

xticks = [500,1000,1500,2000,2500,3000,3500]
xticks_lab = ['500','1000','1500','2000','2500','3000','3500']
plt.xlim(300,3500)
plt.xticks(xticks,xticks_lab)

plt.yscale('log')
yticks = [1000,100,10,1,0.1,0.01,1e-3,1e-4,1e-5,1e-6]#,1e-7,1e-8]
yticks_lab = ['1000','100','10','1','0.1','0.01','10$^{-3}$','10$^{-4}$','10$^{-5}$','10$^{-6}$']#,'10$^{-7}$','10$^{-8}$']
plt.ylim(1000,1e-6)
plt.yticks(yticks,yticks_lab)

#plt.gca().invert_yaxis()
plt.title(title,fontsize=10)
plt.ylabel(r'p$_{\mathrm{gas}}$ [bar]',fontsize=12)
plt.xlabel(r'T$_{\mathrm{gas}}$ [K]',fontsize=12)
CB.set_label(r'$\log_{10}$ k$_{\mathrm{f}}$ [cm$^{6}$s$^{-1}$]',fontsize=12)


plt.tick_params(axis='both',which='major',labelsize=11)

# Save figure
#plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

plt.savefig(fname+'.pdf',dpi=300,bbox_inches='tight')

plt.show()
