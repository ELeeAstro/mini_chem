import numpy as np 
import pickle

plot_EQ = False
output = open('NEQ-HD189-vulcan.txt', "w")


vul = 'HD189.vul'
with open(vul, 'rb') as handle:
  vul = pickle.load(handle)
  
species = vul['variable']['species'] 
out_species = ['OH','H2','H2O','H','CO','CO2','O','CH4','C2H2','NH3','N2','HCN','He']


ost = '{:<8s}'.format('(dyn/cm2)')  + '{:>9s}'.format('(K)') + '{:>9s}'.format('(cm)') + '\n'
ost += '{:<8s}'.format('Pressure')  + '{:>9s}'.format('Temp')+ '{:>9s}'.format('Hight')
for sp in out_species: ost += '{:>10s}'.format(sp) 
ost +='\n'
 
for n, p in enumerate(vul['atm']['pco']):
    ost += '{:<8.3E}'.format(p)  + '{:>8.1f}'.format(vul['atm']['Tco'][n])  + '{:>10.2E}'.format(vul['atm']['zco'][n])
    for sp in out_species:
        if plot_EQ == True:
            ost += '{:>10.2E}'.format(vul['variable']['y_ini'][n,species.index(sp)]/vul['atm']['n_0'][n])
        else: 
            ost += '{:>10.2E}'.format(vul['variable']['ymix'][n,species.index(sp)])
    ost += '\n'

ost = ost[:-1]
output.write(ost)
output.close()

