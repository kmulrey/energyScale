import numpy as np
import os
import cPickle
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '../simulation')
import sim_functions as sim
from scipy.optimize import curve_fit

plt.ion()

data_dir='../../data/radio/sim/'
mag=2.03
average_density = 6.5e-4 # in kg/m^3
bad_event=73889001

def read_file(filename):
    
    infile=open(filename,'r')
    info=cPickle.load(infile)
    infile.close()
    
    event=info['event']
    em_energy=info['em_energy'][event!=bad_event]
    energy=info['energy'][event!=bad_event]
    zenith=info['zenith'][event!=bad_event]
    azimuth=info['azimuth'][event!=bad_event]
    xmax=info['xmax'][event!=bad_event]
    alpha=info['alpha'][event!=bad_event]
    Erad=info['Erad'][event!=bad_event]
    Erad_ce=info['Erad_ce'][event!=bad_event]
    Erad_gm=info['Erad_gm'][event!=bad_event]
    density=info['density'][event!=bad_event]
    density2=info['density2'][event!=bad_event]
    event=event[event!=bad_event]
    
    
    Erad=Erad+Erad*.11-Erad*0.0336
    Erad_ce=Erad_ce+Erad_ce*.11-Erad_ce*0.0336
    Erad_gm=Erad_gm+Erad_gm*.11-Erad_gm*0.0336
    
    return em_energy,energy,zenith,azimuth,xmax,alpha,Erad,Erad_gm,Erad_ce,density,density2,event


em_energy_P,energy_P,zenith_P,azimuth_P,xmax_P,alpha_P,Erad_P,Erad_gm_P,Erad_ce_P,density_P,density2_P,event_P=read_file(data_dir+'compiled_sims_proton_debug2.dat')
em_energy_Fe,energy_Fe,zenith_Fe,azimuth_Fe,xmax_Fe,alpha_Fe,Erad_Fe,Erad_gm_Fe,Erad_ce_Fe,density_Fe,density2_Fe,event_Fe=read_file(data_dir+'compiled_sims_iron_debug2.dat')


em_energy_P_US,energy_P_US,zenith_P_US,azimuth_P_US,xmax_P_US,alpha_P_US,Erad_P_US,Erad_gm_P_US,Erad_ce_P_US,density_P_US,density2_P_US,event_P_US=read_file(data_dir+'compiled_sims_proton_US_standard.dat')
em_energy_Fe_US,energy_Fe_US,zenith_Fe_US,azimuth_Fe_US,xmax_Fe_US,alpha_Fe_US,Erad_Fe_US,Erad_gm_Fe_US,Erad_ce_Fe_US,density_Fe_US,density2_Fe_US,event_Fe_US=read_file(data_dir+'compiled_sims_iron_US_standard.dat')



rho_avg=np.average(np.concatenate([density_P*1e3*np.cos(zenith_P),density_Fe*1e3*np.cos(zenith_Fe)]))
rho_avg_US=np.average(np.concatenate([density_P_US*1e3*np.cos(zenith_P_US),density_Fe_US*1e3*np.cos(zenith_Fe_US)]))

print rho_avg
def a_function(x,p0,p1,p2):
    return p2+p0*np.exp(p1*(x-rho_avg))

y_plt_P_a=np.sin(alpha_P)*np.sqrt(Erad_ce_P/Erad_gm_P)#*mag**0.9
x_plt_P_a=density_P*1e3*np.cos(zenith_P)
y_plt_Fe_a=np.sin(alpha_Fe)*np.sqrt(Erad_ce_Fe/Erad_gm_Fe)#*mag**0.9
x_plt_Fe_a=density_Fe*1e3*np.cos(zenith_Fe)

y_plt_P_a_US=np.sin(alpha_P_US)*np.sqrt(Erad_ce_P_US/Erad_gm_P_US)#*mag**0.9
x_plt_P_a_US=density_P_US*1e3*np.cos(zenith_P_US)
y_plt_Fe_a_US=np.sin(alpha_Fe_US)*np.sqrt(Erad_ce_Fe_US/Erad_gm_Fe_US)#*mag**0.9
x_plt_Fe_a_US=density_Fe_US*1e3*np.cos(zenith_Fe_US)

xdata_a=np.concatenate([x_plt_P_a,x_plt_Fe_a])
ydata_a=np.concatenate([y_plt_P_a,y_plt_Fe_a])

xdata_a_US=np.concatenate([x_plt_P_a_US,x_plt_Fe_a_US])
ydata_a_US=np.concatenate([y_plt_P_a_US,y_plt_Fe_a_US])

popt, pcov = curve_fit(a_function, xdata_a, ydata_a,p0=[0.43,1.11,-0.24], bounds=([0.2,0.8,-0.3], [0.6, 1.5,-0.0]))
p0=popt[0]
p1=popt[1]
p2=popt[2]
print popt

popt_US, pcov_US = curve_fit(a_function, xdata_a_US, ydata_a_US,p0=[0.43,1.11,-0.24], bounds=([0.2,0.8,-0.3], [0.6, 1.5,-0.0]))
p0_US=popt_US[0]
p1_US=popt_US[1]
p2_US=popt_US[2]


# fit from PRL
x_test=np.arange(0,1.2,0.01)
y_test=0.43*(np.exp(1.11*(x_test-0.64)))-.24

y_fit=a_function(x_test,popt[0],popt[1],popt[2])

# plot charge excess ratio



comb_plot_x=np.concatenate([x_plt_P_a,x_plt_Fe_a])
comb_plot_y=np.concatenate([y_plt_P_a,y_plt_Fe_a])
comb_plot_c=np.concatenate([np.sin(alpha_P),np.sin(alpha_Fe)])
#comb_plot_c=np.concatenate([np.degrees(zenith_P),np.degrees(zenith_Fe)])

#comb_plot_c=comb_plot_c[comb_plot_y>0.05]
#comb_plot_x=comb_plot_x[comb_plot_y>0.05]
#comb_plot_y=comb_plot_y[comb_plot_y>0.05]



min=0.4#np.degrees(0)#0.4
max=0.9#np.degrees(1)#0.9#

fig = plt.figure(figsize=(6.7,5))
ax1 = fig.add_subplot(1,1,1)


#sc = ax1.scatter(comb_plot_x,comb_plot_y, alpha=0.6,s=1,vmin=min,vmax=max,c=comb_plot_c)
#cbar = plt.colorbar(sc)
#cbar.ax.set_ylabel(r'sin $\alpha$', rotation=90)

ax1.plot(x_plt_P_a,y_plt_P_a,'.',alpha=0.2,color='blue',label='P')
ax1.plot(x_plt_Fe_a,y_plt_Fe_a,'.',alpha=0.2,color='red',label='Fe')

#sc = ax1.scatter(x_plt_Fe,y_plt_Fe, alpha=0.6,s=1,vmin=np.radians(0),vmax=np.radians(60),c=zenith_Fe)
#ax1.plot(x_test,y_fit,color='black',label=r'$p_2+p_0 exp(p_1 \cdot(\rho - \left< \rho \right>))$, $p_0$={0:.2f}, $p_1$={1:.2f}, $p_2$={2:.2f}'.format(p0,p1,p2))
#ax2.plot(x_test,y_fit,color='black',label='fit')
#ax1.plot(x_test,y_test,color='red',label='AERA PRL, $p_0$={0:.2f}, $p_1$={1:.2f}, $p_2$={2:.2f}'.format(0.43,1.11,-0.24))

#cmap='gnuplot2_r'#

#ax2.hist2d(x_plt_P_a,y_plt_P_a, bins=(100, 50),cmap=cmap)

ax1.set_xlim([0.3,1.3])
ax1.set_ylim([0.03,.4])
#ax1.set_ylim([0.03,.7])

ax1.grid()

#ax2.set_xlim([0.2,1.3])
#ax2.set_ylim([0.0,.3])
#ax2.grid()

ax1.legend(numpoints=1,loc='upper left',framealpha=1,facecolor='white')
#ax2.legend(numpoints=1,loc='upper left',framealpha=1,facecolor='white')

ax1.set_xlabel('density at shower maximum (kg/m$^3$)',fontsize=12)
ax1.set_ylabel(r'a=sin $\alpha\sqrt{E^{ce}_{rad}/E^{geo}_{rad}}$',fontsize=12)# (B_{LOFAR}/B_{Auger})^{0.9}


#ax2.set_xlabel('density at shower maximum (kg/m$^3$)',fontsize=12)
#ax2.set_ylabel(r'a=sin $\alpha\sqrt{E^{ce}_{rad}/E^{geo}_{rad}}$',fontsize=12)

fig.tight_layout()
plt.show()
raw_input()

plt.close()

def return_a(rho,avg,p0,p1,p2):
    a=np.zeros(len(rho))
    for i in np.arange(len(rho)):
        a[i]= p2+p0*np.exp(p1*(rho[i]-avg))
    return a


a_P=return_a(density_P*1e3*np.cos(zenith_P),rho_avg,p0,p1,p2)/mag**0.9
a_Fe=return_a(density_Fe*1e3*np.cos(zenith_Fe),rho_avg,p0,p1,p2)/mag**0.9

a_P_US=return_a(density_P_US*1e3*np.cos(zenith_P_US),rho_avg_US,p0_US,p1_US,p2_US)/mag**0.9
a_Fe_US=return_a(density_Fe_US*1e3*np.cos(zenith_Fe_US),rho_avg_US,p0_US,p1_US,p2_US)/mag**0.9


Srd_P=Erad_P/np.sin(alpha_P)**2/mag**1.8
Srd_Fe=Erad_Fe/np.sin(alpha_Fe)**2/mag**1.8

Srd_P_1=Erad_P/(a_P**2+(1-a_P**2)*np.sin(alpha_P)**2)/mag**1.8
Srd_Fe_1=Erad_Fe/(a_Fe**2+(1-a_Fe**2)*np.sin(alpha_Fe)**2)/mag**1.8

Srd_P_US=Erad_P_US/np.sin(alpha_P_US)**2/mag**1.8
Srd_Fe_US=Erad_Fe_US/np.sin(alpha_Fe_US)**2/mag**1.8

Srd_P_1_US=Erad_P_US/(a_P_US**2+(1-a_P_US**2)*np.sin(alpha_P_US)**2)/mag**1.8
Srd_Fe_1_US=Erad_Fe_US/(a_Fe_US**2+(1-a_Fe_US**2)*np.sin(alpha_Fe_US)**2)/mag**1.8



a_em_p,b_em_p,a_sig_em_p,b_sig_em_p=sim.getFit(em_energy_P,Srd_P_1)
a_em_fe,b_em_fe,a_sig_em_fe,b_sig_em_fe=sim.getFit(em_energy_Fe,Srd_Fe_1)

a_em_p_US,b_em_p_US,a_sig_em_p_US,b_sig_em_p_US=sim.getFit(em_energy_P_US,Srd_P_1_US)
a_em_fe_US,b_em_fe_US,a_sig_em_fe_US,b_sig_em_fe_US=sim.getFit(em_energy_Fe_US,Srd_Fe_1_US)

print a_em_p,b_em_p
print a_em_fe,b_em_fe

EM_sr_P=1e18*np.power((Srd_P_1/(a_em_p*1e7)),(1/b_em_p))
EM_sr_Fe=1e18*np.power((Srd_Fe_1/(a_em_fe*1e7)),(1/b_em_fe))

EM_sr_P_US=1e18*np.power((Srd_P_1_US/(a_em_p_US*1e7)),(1/b_em_p_US))
EM_sr_Fe_US=1e18*np.power((Srd_Fe_1_US/(a_em_fe_US*1e7)),(1/b_em_fe_US))


def e_function(x,p0,p1,p2):
    return p2+p0*np.exp(p1*(x-rho_avg))


x_plt_P=density_P*1e3*np.cos(zenith_P)
y_plt_P=EM_sr_P/em_energy_P
x_plt_Fe=density_Fe*1e3*np.cos(zenith_Fe)
y_plt_Fe=EM_sr_Fe/em_energy_Fe


x_plt_P_US=density_P_US*1e3*np.cos(zenith_P_US)
y_plt_P_US=EM_sr_P_US/em_energy_P_US
x_plt_Fe_US=density_Fe_US*1e3*np.cos(zenith_Fe_US)
y_plt_Fe_US=EM_sr_Fe_US/em_energy_Fe_US


xdata=np.concatenate([x_plt_P,x_plt_Fe])
ydata=np.concatenate([y_plt_P,y_plt_Fe])
cdata=np.concatenate([np.sin(alpha_P),np.sin(alpha_Fe)])


xdata_US=np.concatenate([x_plt_P_US,x_plt_Fe_US])
ydata_US=np.concatenate([y_plt_P_US,y_plt_Fe_US])
cdata_US=np.concatenate([np.sin(alpha_P_US),np.sin(alpha_Fe_US)])
#cdata=np.concatenate([zenith_P,zenith_Fe])

popt, pcov = curve_fit(a_function, xdata, ydata,p0=[0.25,-3.08,0.77], bounds=([0.1,-5.0,0.2], [1.8, 0.0,1.2]))
print popt
p0e=popt[0]
p1e=popt[1]
p2e=popt[2]

popt_p, pcov_p = curve_fit(a_function, x_plt_P,y_plt_P,p0=[0.25,-3.08,0.77], bounds=([0.1,-5.0,0.0], [1.8, 0.0,1.2]))
popt_fe, pcov = curve_fit(a_function, x_plt_Fe,y_plt_Fe,p0=[0.25,-3.08,0.77], bounds=([0.1,-5.0,0.0], [1.8, 0.0,1.2]))
print popt_fe

p0e_p=popt_p[0]
p1e_p=popt_p[1]
p2e_p=popt_p[2]

p0e_fe=popt_fe[0]
p1e_fe=popt_fe[1]
p2e_fe=popt_fe[2]


popt_US, pcov_US = curve_fit(a_function, xdata_US, ydata_US,p0=[0.25,-3.08,0.77], bounds=([0.1,-5.0,0.2], [1.8, 0.0,1.2]))


popt_p_US, pcov_p_US = curve_fit(a_function, x_plt_P_US,y_plt_P_US,p0=[0.25,-3.08,0.77], bounds=([0.1,-5.0,0.0], [1.8, 0.0,1.2]))
popt_fe_US, pcov_US = curve_fit(a_function, x_plt_Fe_US,y_plt_Fe_US,p0=[0.25,-3.08,0.77], bounds=([0.1,-5.0,0.0], [1.8, 0.0,1.2]))

p0e_US=popt_US[0]
p1e_US=popt_US[1]
p2e_US=popt_US[2]

p0e_p_US=popt_p_US[0]
p1e_p_US=popt_p_US[1]
p2e_p_US=popt_p_US[2]

p0e_fe_US=popt_fe_US[0]
p1e_fe_US=popt_fe_US[1]
p2e_fe_US=popt_fe_US[2]

x_test=np.arange(0,1.2,0.01)
y_test=0.25*(np.exp(-3.08*(x_test-0.65)))+0.77
y_fit=p0e*(np.exp(p1e*(x_test-rho_avg)))+p2e
y_fit_p=p0e_p*(np.exp(p1e_p*(x_test-rho_avg)))+p2e_p
y_fit_fe=p0e_fe*(np.exp(p1e_fe*(x_test-rho_avg)))+p2e_fe

y_test_US=0.25*(np.exp(-3.08*(x_test-0.65)))+0.77
y_fit_US=p0e_US*(np.exp(p1e_US*(x_test-rho_avg_US)))+p2e_US
y_fit_p_US=p0e_p_US*(np.exp(p1e_p_US*(x_test-rho_avg_US)))+p2e_p_US
y_fit_fe_US=p0e_fe_US*(np.exp(p1e_fe_US*(x_test-rho_avg_US)))+p2e_fe_US


#  plot charge excess ratio vs density

fig = plt.figure()
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)

#sc1 = ax1.scatter(x_plt_P,y_plt_P, alpha=0.6,s=1,vmin=np.radians(0),vmax=np.radians(60),c=zenith_P)
#sc2 = ax1.scatter(x_plt_Fe,y_plt_Fe, alpha=0.6,s=1,vmin=np.radians(0),vmax=np.radians(60),c=zenith_Fe)
#sc = ax1.scatter(x_plt_P,y_plt_P, alpha=0.6,s=1,vmin=min,vmax=max,c=np.sin(alpha_P))

#sc = ax1.scatter(xdata,ydata, alpha=0.6,s=1,vmin=min,vmax=max,c=cdata)
#cbar = plt.colorbar(sc)
#cbar.ax.set_ylabel(r'sin $\alpha$', rotation=90)
#cbar.ax.set_ylabel(r'zenith $\theta$', rotation=90)
ax1.plot(x_plt_P,y_plt_P,'.',alpha=0.2,color='blue',label='P')
ax1.plot(x_plt_Fe,y_plt_Fe,'.',alpha=0.2,color='red',label='Fe')
ax1.legend(numpoints=1,loc='upper left',framealpha=1,facecolor='white')

ax1.set_xlabel('density at shower maximum (kg/m$^3$)',fontsize=12)
ax1.set_ylabel(r'$E_{em}(S_{RD})/E_{em}$',fontsize=12)# (B_{LOFAR}/B_{Auger})^{0.9}


ax1.set_xlim([0.2,1.3])
ax1.set_ylim([0.6,2.0])



ax2.plot(x_plt_P_US,y_plt_P_US,'.',alpha=0.2,color='blue',label='P')
ax2.plot(x_plt_Fe_US,y_plt_Fe_US,'.',alpha=0.2,color='red',label='Fe')
ax2.legend(numpoints=1,loc='upper left',framealpha=1,facecolor='white')

ax2.set_xlabel('density at shower maximum (kg/m$^3$)',fontsize=12)
ax2.set_ylabel(r'$E_{em}(S_{RD})/E_{em}$',fontsize=12)# (B_{LOFAR}/B_{Auger})^{0.9}

ax2.set_xlim([0.2,1.3])
ax2.set_ylim([0.6,2.0])




ax1.grid()
ax2.grid()

fig.tight_layout()
plt.show()
raw_input()
plt.close()


fig = plt.figure(figsize=(10,5))
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)

cmap='gnuplot2_r'#

ax1.hist2d(xdata,ydata, bins=(100, 50),cmap=cmap,label='Fe')

#ax1.hist2d(x_plt_P,y_plt_P, bins=(100, 50),cmap=cmap,label='P')
##ax1.hist2d(x_plt_Fe,y_plt_Fe, bins=(100, 50),cmap=cmap,label='Fe')

ax1.plot(x_test,y_fit_p,color='blue',label=r'proton')#: $p_2+p_0 exp(p_1 \cdot(\rho - \left< \rho \right>))$, $p_0$={0:.2f}, $p_1$={1:.2f}, $p_2$={2:.2f}'.format(p0e_p,p1e_p,p2e_p))
ax1.plot(x_test,y_fit_fe,color='red',label=r'iron')#: $p_2+p_0 exp(p_1 \cdot(\rho - \left< \rho \right>))$, $p_0$={0:.2f}, $p_1$={1:.2f}, $p_2$={2:.2f}'.format(p0e_fe,p1e_fe,p2e_fe))
ax1.plot(x_test,y_test,color='black',label=r'PRL')#: $p_2+p_0 exp(p_1 \cdot(\rho - \left< \rho \right>))$, $p_0$={0:.2f}, $p_1$={1:.2f}, $p_2$={2:.2f}'.format(p0e_p,p1e_p,p2e_p))


ax2.hist2d(xdata_US,ydata_US, bins=(100, 50),cmap=cmap,label='Fe')

ax2.plot(x_test,y_fit_p_US,color='blue',label=r'proton')#: $p_2+p_0 exp(p_1 \cdot(\rho - \left< \rho \right>))$, $p_0$={0:.2f}, $p_1$={1:.2f}, $p_2$={2:.2f}'.format(p0e_p,p1e_p,p2e_p))
ax2.plot(x_test,y_fit_fe_US,color='red',label=r'iron')#: $p_2+p_0 exp(p_1 \cdot(\rho - \left< \rho \right>))$, $p_0$={0:.2f}, $p_1$={1:.2f}, $p_2$={2:.2f}'.format(p0e_fe,p1e_fe,p2e_fe))
ax2.plot(x_test,y_test,color='black',label=r'PRL')
#ax1.plot(x_test,y_test,color='red',label='AERA PRL, $p_0$={0:.2f}, $p_1$={1:.2f}, $p_2$={2:.2f}'.format(0.25,-3.08,0.77))

#ax1.plot(x_test,y_fit,color='black',label=r'$p_2+p_0 exp(p_1 \cdot(\rho - \left< \rho \right>))$, $p_0$={0:.2f}, $p_1$={1:.2f}, $p_2$={2:.2f}'.format(p0e,p1e,p2e))
#ax1.plot(x_test,y_test,color='red',label='AERA PRL, $p_0$={0:.2f}, $p_1$={1:.2f}, $p_2$={2:.2f}'.format(0.25,-3.08,0.77))


ax1.legend(numpoints=1,loc='upper left',framealpha=1,facecolor='white')

ax1.set_xlabel('density at shower maximum (kg/m$^3$)',fontsize=12)
ax1.set_ylabel(r'$E_{em}(S_{RD})/E_{em}$',fontsize=12)# (B_{LOFAR}/B_{Auger})^{0.9}


ax1.set_xlim([0.2,1.3])
ax1.set_ylim([0.6,2.0])


ax2.legend(numpoints=1,loc='upper left',framealpha=1,facecolor='white')

ax2.set_xlabel('density at shower maximum (kg/m$^3$)',fontsize=12)
ax2.set_ylabel(r'$E_{em}(S_{RD})/E_{em}$',fontsize=12)# (B_{LOFAR}/B_{Auger})^{0.9}


ax2.set_xlim([0.2,1.3])
ax2.set_ylim([0.6,2.0])





ax1.grid()
ax2.grid()

#plt.colorbar(sc1)
fig.tight_layout()
plt.show()
raw_input()
plt.close()


p0a=0.25
p1a=-3.08
p2a=0.77


Srd_P_2=Erad_P/(a_P**2+(1-a_P**2)*np.sin(alpha_P)**2)/(1-p0e_p+p0e_p*np.exp(p1e_p*(density_P*1e3*np.cos(zenith_P)-rho_avg)))**2/mag**1.8
Srd_Fe_2=Erad_Fe/(a_Fe**2+(1-a_Fe**2)*np.sin(alpha_Fe)**2)/(1-p0e_fe+p0e_fe*np.exp(p1e_fe*(density_Fe*1e3*np.cos(zenith_Fe)-rho_avg)))**2/mag**1.8

#Srd_P_2_a=Erad_P/(a_P**2+(1-a_P**2)*np.sin(alpha_P)**2)/(1-p0a+p0a*np.exp(p1a*(density_P*1e3*np.cos(zenith_P)-rho_avg)))**2/mag**1.8
#Srd_Fe_2_a=Erad_Fe/(a_Fe**2+(1-a_Fe**2)*np.sin(alpha_Fe)**2)/(1-p0a+p0a*np.exp(p1a*(density_Fe*1e3*np.cos(zenith_Fe)-rho_avg)))**2/mag**1.8

a_em_p,b_em_p,a_sig_em_p,b_sig_em_p=sim.getFit(em_energy_P,Srd_P_2)
a_em_fe,b_em_fe,a_sig_em_fe,b_sig_em_fe=sim.getFit(em_energy_Fe,Srd_Fe_2)

a_tot_p,b_tot_p,a_sig_tot_p,b_sig_tot_p=sim.getFit(energy_P,Srd_P_2)
a_tot_fe,b_tot_fe,a_sig_tot_fe,b_sig_tot_fe=sim.getFit(energy_Fe,Srd_Fe_2)

print '\n\nEM energy fit'
print '___________________'
print '{0:.3f} +/- {1:.3f}       {2:.3f} +/- {3:.3f}'.format(a_em_p,a_sig_em_p,b_em_p,b_sig_em_p)
print '{0:.3f} +/- {1:.3f}       {2:.3f} +/- {3:.3f}'.format(a_em_fe,a_sig_em_fe,b_em_fe,b_sig_em_fe)


print '\n\nTotal energy fit'
print '___________________'
print '{0:.3f} +/- {1:.3f}       {2:.3f} +/- {3:.3f}'.format(a_tot_p,a_sig_tot_p,b_tot_p,b_sig_tot_p)
print '{0:.3f} +/- {1:.3f}       {2:.3f} +/- {3:.3f}'.format(a_tot_fe,a_sig_tot_fe,b_tot_fe,b_sig_tot_fe)

x = np.linspace(10**16,10**18.5,100)
x_data = np.linspace(10**17,10**18,100)

#EM_sr_P=1e18*np.power((Srd_P_2/(a_em_p*1e7)),(1/b_em_p))
#EM_sr_Fe=1e18*np.power((Srd_P_2/(a_em_fe*1e7)),(1/b_em_fe))

'''
rad_em_fit_p=sim.energy_to_rad(x,a_em_p,b_em_p)
rad_em_fit_fe=sim.energy_to_rad(x,a_em_fe,b_em_fe)

rad_tot_fit_p=sim.energy_to_rad(x,a_tot_p,b_tot_p)
rad_tot_fit_fe=sim.energy_to_rad(x,a_tot_fe,b_tot_fe)

em_recon_0=sim.rad_to_energy(Srd_P,a_em_p,b_em_p)
em_recon_1=sim.rad_to_energy(Srd_P_1,a_em_p,b_em_p)
em_recon_2=sim.rad_to_energy(Srd_P_2,a_em_p,b_em_p)


dist0=em_recon_0[xmax_P<870]/em_energy_P[xmax_P<870]
dist1=em_recon_1[xmax_P<870]/em_energy_P[xmax_P<870]
dist2=em_recon_2[xmax_P<870]/em_energy_P[xmax_P<870]

bins=np.arange(0.8,1.2,0.003)
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)

ax1.axvline(x=1,ymax=450,linestyle=':',color='black')


ax1.hist(dist0,bins=bins,histtype='step',color='black',label='no corrections')
ax1.hist(dist1,bins=bins,histtype='step',color='blue',label='charge excess correction')
ax1.hist(dist2,bins=bins,histtype='step',color='green',label='slant depth correction')
ax1.legend(loc='upper left',framealpha=1,facecolor='white')
ax1.set_xlabel(r'$E_{em}(S_{RD})/E_{em}$',fontsize=12)

fig.tight_layout()

plt.show()
raw_input()

plt.close()
'''



'''
fig = plt.figure()#figsize=(8,6))
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)

#ax1.plot(em_energy_P[xmax_P<870],Srd_P[xmax_P<870],'.',color='black',alpha=0.5,label='P, no correction')
#ax1.plot(em_energy_P[xmax_P<870],Srd_P_1[xmax_P<870],'.',color='blue',alpha=0.5,label='P, charge excess correction')
#ax1.plot(em_energy_P[xmax_P<870],Srd_P_2[xmax_P<870],'.',color='green',alpha=0.5,label='P, slant depth correction')


#sc1 = ax1.scatter(em_energy_P,Srd_P_2, alpha=0.6,s=1,vmin=500,vmax=900,c=xmax_P) #vmin=np.radians(0),vmax=np.radians(60)
#sc2 = ax1.scatter(em_energy_Fe,Srd_Fe_2, alpha=0.6,s=1,vmin=np.radians(0),vmax=np.radians(60),c=alpha_Fe) #vmin=np.radians(0),vmax=np.radians(60)

ax1.plot(em_energy_P[xmax_P<870],Srd_P_2[xmax_P<870],'.',color='blue',alpha=0.2,label='P')
ax1.plot(em_energy_Fe,Srd_Fe_2,'.',color='red',alpha=0.2,label='Fe')

ax2.plot(energy_P[xmax_P<870],Srd_P_2[xmax_P<870],'.',color='blue',alpha=0.2,label='P')
ax2.plot(energy_Fe,Srd_Fe_2,'.',color='red',alpha=0.2,label='Fe')

ax1.plot(x,rad_em_fit_p,':',color='black',label=r'P: A={0:.3f}, B={1:.3f}'.format(a_em_p,b_em_p))
ax1.plot(x,rad_em_fit_fe,color='black',label=r'Fe: A={0:.3f}, B={1:.3f}'.format(a_em_fe,b_em_fe))

ax2.plot(x,rad_tot_fit_p,':',color='black',label=r'P: A={0:.3f}, B={1:.3f}'.format(a_tot_p,b_tot_p))
ax2.plot(x,rad_tot_fit_fe,color='black',label=r'Fe: A={0:.3f}, B={1:.3f}'.format(a_tot_fe,b_tot_fe))

ax1.legend(numpoints=1,loc='upper left',framealpha=1,facecolor='white')
ax2.legend(numpoints=1,loc='upper left',framealpha=1,facecolor='white')

ax1.set_xlabel('energy in EM cascade (eV)',fontsize=12)
ax1.set_ylabel(r'corrected radiation energy $S^{\rho}_{RD}$ (eV)',fontsize=12)# (B_{LOFAR}/B_{Auger})^{0.9}
ax2.set_xlabel('energy in total cascade (eV)',fontsize=12)
ax2.set_ylabel(r'corrected radiation energy $S^{\rho}_{RD}$ (eV)',fontsize=12)# (B_{LOFAR}/B_{Auger})^{0.9}

#ax1.set_xlabel('energy in EM cascade (eV)',fontsize=12)
#ax1.set_ylabel(r'radiation energy (eV)',fontsize=12)# (B_{LOFAR}/B_{Auger})^{0.9}

#ax1.legend(numpoints=1,loc='upper left',framealpha=1,facecolor='white')

ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.axis([1e16,3e18,2e3,1e8])
ax1.grid()

ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.axis([1e16,3e18,2e3,1e8])
ax2.grid()
fig.tight_layout()

plt.show()
raw_input()

plt.close()
'''
