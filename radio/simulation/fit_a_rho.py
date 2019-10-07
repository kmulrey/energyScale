import numpy as np
import os
import cPickle
import matplotlib.pyplot as plt
import sys
import tools as tools

import fluence as flu
import radiation_energy as rad
import sim_functions as sim
import helper as helper
from scipy.optimize import minimize
import matplotlib.ticker as ticker

plt.ion()

b_scale=2.03
b=1.8
mag=2.03
average_density = 6.5e-4 #in g/cm^3#atmc.get_density(average_zenith, average_xmax) * 1e-3  # in kg/m^3
p0=0.250524463912
p1=-2.95290494

p0_=0.239
p1_=-3.13

data_dir='/Users/kmulrey/LOFAR/energy/energyScale/data/radio/sim/'

coreas_file=data_dir+'compiled_sims_proton.dat'

em_energy_P,energy_P,zenith_P,azimuth_P,xmax_P,alpha_P,S_basic_P,Srd_1_P,Srd_2_P,Erad_P,a_P,density_P,Erad_vxB_P,Erad_vxvxB_P=tools.read_file_plain(data_dir+'compiled_sims_proton.dat')
em_energy_Fe,energy_Fe,zenith_Fe,azimuth_Fe,xmax_Fe,alpha_Fe,S_basic_Fe,Srd_1_Fe,Srd_2_Fe,Erad_Fe,a_Fe,density_Fe,Erad_vxB_Fe,Erad_vxvxB_Fe=tools.read_file_plain(data_dir+'compiled_sims_iron.dat')

def get_Srd_alpha(Erad,alpha):
    return Erad/(np.sin(alpha)**2*mag**b)

def get_Srd(Erad,alpha,a):
    a1=a/(mag**0.9)
    return Erad/(a1**2+(1-a1**2)*(np.sin(alpha)**2*mag**b))

def get_Srd(Erad,alpha,a,density):
    p0=0.25
    p1=-3.08
    avg_rho=0.65 #km/m3
    a1=a/(mag**0.9)

    return Erad/((a1**2+(1-a1**2)*(np.sin(alpha)**2*mag**b))*(a-p0+p0*np.exp(p1*(density*1e3-avg_rho))**2))


#x = np.linspace(10**4,10**8,100)

y_plt=np.sin(alpha_P)*np.sqrt(Erad_vxvxB_P/Erad_vxB_P)*mag**0.9
x_plt=density_P*1e3

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)


ax1.plot(x_plt,y_plt,'.')



ax1.grid()
plt.show()
raw_input()
plt.close()


'''
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1,aspect=1)


ax1.plot(em_energy_P,Erad_P,'.',color='blue',alpha=0.2,label='proton')
#ax1.plot(em_energy_P,get_Srd_alpha(Erad_P,alpha_P),'.',color='red',alpha=0.2,label='proton')
ax1.plot(em_energy_P,Erad_P/(a_P**2+(1-a_P**2)*(np.sin(alpha_P)**2*mag**b)),'.',color='green',alpha=0.2,label='proton')

ax1.plot(em_energy_P,get_Srd(Erad_P,alpha_P,a_P,density_P),'.',color='red',alpha=0.2,label='proton')

#ax1.plot(em_energy_P,Srd_1_P,'.',color='red',alpha=0.2,label='proton')

#ax1.plot(em_energy_P,Srd_2_P,'.',color='purple',alpha=0.2,label='proton')



ax1.set_yscale('log')
ax1.set_xscale('log')



ax1.grid()
plt.show()
raw_input()
plt.close()
'''

'''

#ax1.plot(x,x,color='black')

ax1.plot(em_energy_P,S_P_use,'.',color='blue',alpha=0.2,label='proton')
ax1.plot(x,rad_em_fit_p,color='blue')

ax1.plot(em_energy_Fe,S_Fe_use,'.',color='red',alpha=0.2,label='iron')
ax1.plot(x,rad_em_fit_fe,color='red')

ax1.set_yscale('log')
ax1.set_xscale('log')

ax1.legend(numpoints=1,loc='upper left',framealpha=1,facecolor='white')

ax1.set_xlabel('EM Energy (eV)',fontsize=12)
ax1.set_ylabel('S$_{RD} \ (eV)$',fontsize=12)

ax1.axis([2e16,3e18,7e3,1e8])


###
ax2.plot(energy_P,S_P_use,'.',color='blue',alpha=0.2,label='proton')
ax2.plot(x,rad_tot_fit_p,color='blue')


ax2.plot(energy_Fe,S_Fe_use,'.',color='red',alpha=0.2,label='iron')
ax2.plot(x,rad_tot_fit_fe,color='red')

ax2.set_yscale('log')
ax2.set_xscale('log')

ax2.legend(numpoints=1,loc='upper left',framealpha=1,facecolor='white')

ax2.set_xlabel('Total CR Energy (eV)',fontsize=12)
ax2.set_ylabel('S$_{RD} (eV)$',fontsize=12)
ax2.axis([2e16,3e18,7e3,1e8])

ax1.grid()
ax2.grid()
plt.show()
raw_input()
plt.close()
'''
