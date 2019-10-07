import numpy as np
import os
import cPickle
import matplotlib.pyplot as plt
import sys

sys.path.insert(0, '../simulation')
import sim_functions as sim


plt.ion()
#
data_dir='../../data/radio/sim_local/'
#data_dir=''

mag=2.03
average_density = 6.5e-4 #in g/cm^3#atmc.get_density(average_zenith, average_xmax) * 1e-3  # in kg/m^3
p0=0.250524463912
p1=-2.95290494

p0_=0.239
p1_=-3.13

def read_file_old(filename):
    
    infile=open(filename,'r')
    info=cPickle.load(infile)
    infile.close()
    
    
    em_energy=info['em_energy']
    energy=info['energy']
    zenith=info['zenith']
    azimuth=info['azimuth']
    xmax=info['xmax']
    alpha=info['alpha']
    Erad=info['Erad']
    Erad_ce=info['Erad_ce']
    Erad_gm=info['Erad_gm']
    density=info['density']
    
    #Erad_ce=Erad_ce+Erad_ce*.51
    
    #Erad=Erad+Erad*.11-Erad*0.0336
    #Erad_ce=Erad_ce+Erad_ce*.11-Erad_ce*0.0336
    #Erad_gm=Erad_gm+Erad_gm*.01#-Erad_gm*0.0336
    
    return em_energy,energy,zenith,azimuth,xmax,alpha,Erad,Erad_gm,Erad_ce,density


def read_file(filename):
    
    infile=open(filename,'r')
    info=cPickle.load(infile)
    infile.close()
    
    
    em_energy=info['em_energy']
    energy=info['energy']
    zenith=info['zenith']
    azimuth=info['azimuth']
    xmax=info['xmax']
    alpha=info['alpha']
    Erad=info['Erad']
    Erad_ce=info['Erad_ce']
    Erad_gm=info['Erad_gm']
    density=info['density']
    density2=info['density2']

    #Erad_ce=Erad_ce+Erad_ce*.51

    #Erad=Erad+Erad*.11-Erad*0.0336
    #Erad_ce=Erad_ce+Erad_ce*.11-Erad_ce*0.0336
    #Erad_gm=Erad_gm+Erad_gm*.01#-Erad_gm*0.0336

    return em_energy,energy,zenith,azimuth,xmax,alpha,Erad,Erad_gm,Erad_ce,density,density2


em_energy_P,energy_P,zenith_P,azimuth_P,xmax_P,alpha_P,Erad_P,Erad_gm_P,Erad_ce_P,density_P,density2_P=read_file('compiled_sims_proton_debug2.dat')
em_energy_Fe,energy_Fe,zenith_Fe,azimuth_Fe,xmax_Fe,alpha_Fe,Erad_Fe,Erad_gm_Fe,Erad_ce_Fe,density_Fe=read_file_old(data_dir+'compiled_sims_proton.dat')

x_test=np.arange(0,1.2,0.01)
y_test=0.43*(np.exp(1.11*(x_test-0.65)))-.24
y_plt_P=np.sin(alpha_Fe)*np.sqrt(Erad_ce_Fe/Erad_gm_Fe)*mag**0.9#*mag**0.9=np.sin(alpha_P)*
x_plt_P=density_Fe*1e3
cl_P=np.cos(zenith_P)
y_plt_Fe=np.sin(alpha_P)*np.sqrt(Erad_ce_P/Erad_gm_P)*mag**0.9#*mag**0.9=np.sin(alpha_Fe)*
x_plt_Fe=density2_P*1e3#*np.cos(zenith_P)
cl_Fe=(zenith_Fe)

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)

ax1.plot(x_plt_P,y_plt_P,'.')
ax1.plot(x_plt_Fe,y_plt_Fe,'x')

#ax1.scatter(x_plt_P,y_plt_P,s=10,c=cl_P,vmin=np.sqrt(2)/2,vmax=1)
#ax1.scatter(x_plt_Fe,y_plt_Fe,s=10,c=cl_Fe,vmin=np.sqrt(2)/2,vmax=1)
ax1.plot(x_test,y_test,color='black')



ax1.grid()
plt.show()
raw_input()
plt.close()
