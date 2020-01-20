import numpy as np
import os
import cPickle
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, 'simulation/')
import tools as tools

import fluence as flu
import radiation_energy as rad
import sim_functions as sim
import helper as helper
from scipy.optimize import minimize
import matplotlib.ticker as ticker


plt.ion()


mag=2.03

data_dir='../data/radio/lofar_events/'

rho_avg=0.72  # average density, 0.72 kg/m3
def return_a(rho,avg,p0,p1,p2):
    a= p2+p0*np.exp(p1*(rho-avg))
    return a

# charge excess parameters
p0a=0.2
p1a=1.27
p2a=-0.08

# slant depth parameters proton
p0s_p=1.02
p1s_p=-0.49
p2s_p=0

# slant depth parameters iron
p0s_fe=0.98
p1s_fe=-0.47
p2s_fe=0

# energy parameters proton no correction
#A_p=1.292
#B_p=2.035
#A_em_p=1.669
#B_em_p=2.007

# energy parameters iron no correction
#A_fe=1.100
#B_fe=2.079
#A_em_fe=1.625
#B_em_fe=2.02

# energy parameters proton 11% correction
A_p=1.391
B_p=2.035
A_em_p=1.797
B_em_p=2.007

# energy parameters iron 11% correction
A_fe=1.184
B_fe=2.025
A_em_fe=1.749
B_em_fe=2.02


coreas_file=data_dir+'compiled_sims_all_final.dat'
ldf_file=data_dir+'parameterization_radiation_energy.p'
uncertainty_file=data_dir+'uncertainty_info.dat'

event,Erad,Erad_ce,Erad_gm,azimuth,zenith,xmax,em_energy,energy,density,p_ratio,d_ratio,type,alpha,sigma_e,sigma_e_radio,sigma_r=tools.combine_files(coreas_file,ldf_file,uncertainty_file)




### ldf already has div, sin2alpha

particle_scale=(6.5/6.7)#(1.1)*(6.1/6.7)



Erad=Erad*p_ratio
Erad=Erad+Erad*0.11-Erad*0.0336
Eparticle=energy*d_ratio*particle_scale
Eradio=energy*np.sqrt(p_ratio)
Eradio=Eradio+Eradio*0.05
em_energy=em_energy*np.sqrt(p_ratio)
cut=0.4
#d_ratio_use=d_ratio[sigma_e<cut]
#p_ratio_use=p_ratio[sigma_e<cut]

################################### transform sim into energy ###################################

Srd=np.zeros(len(event))
#print Srd_corr.shape


def return_Srd(Erad,zenith,density,aplpha,type):
    Srd_2=np.zeros(len(Erad))
    for i in np.arange(len(Erad)):
        a=return_a(density[i]*1e3*np.cos(zenith[i]),rho_avg,p0a,p1a,p2a)/mag**0.9
        Srd=Erad[i]/np.sin(alpha[i])**2/mag**1.8
        Srd_1=Erad[i]/(a**2+(1-a**2)*np.sin(alpha[i])**2)/mag**1.8
        if type[i]==0:
            Srd_2[i]=Erad[i]/(a**2+(1-a**2)*np.sin(alpha[i])**2)/(1-p0s_p+p0s_p*np.exp(p1s_p*(density[i]*1e3*np.cos(zenith[i])-rho_avg)))**2/mag**1.8
        else:
            Srd_2[i]=Erad[i]/(a**2+(1-a**2)*np.sin(alpha[i])**2)/(1-p0s_fe+p0s_fe*np.exp(p1s_fe*(density[i]*1e3*np.cos(zenith[i])-rho_avg)))**2/mag**1.8

    return Srd_2

Srd=return_Srd(Erad,zenith,density,alpha,type)


em_energy_Srd=np.zeros(len(Srd))
energy_Srd=np.zeros(len(Srd))

for i in np.arange(len(Srd)):
    if type[i]==0:
        em_energy_Srd[i]=tools.rad_to_energy(Srd[i],A_em_p,B_em_p)
        energy_Srd[i]=tools.rad_to_energy(Srd[i],A_p,B_p)

    if type[i]==1:
        em_energy_Srd[i]=tools.rad_to_energy(Srd[i],A_em_fe,B_em_fe)
        energy_Srd[i]=tools.rad_to_energy(Srd[i],A_fe,B_fe)



line=np.arange(1e16,1e19,1e16)
#line=np.arange(4e4,4e8,1e8)

x_fit = np.linspace(10**16,10**18.5,100)

rad_p_fit=sim.energy_to_rad(x_fit,A_p,B_p)
rad_fe_fit=sim.energy_to_rad(x_fit,A_fe,B_fe)

x=Eparticle
y=Srd


x_err=x*sigma_e#*d_ratio
y_err=y*sigma_e_radio#*np.sqrt(p_ratio)
#x=x[sigma_e<0.7]
#y=y[sigma_e<0.5]
#x_err=x_err[sigma_e<0.5]
#y_err=y_err[sigma_e<0.5]
#print x

########### lines from predictions
x_fit = np.linspace(10**16,10**18.5,1000)

rad_p_fit=sim.energy_to_rad(x_fit,A_p,B_p)
rad_fe_fit=sim.energy_to_rad(x_fit,A_fe,B_fe)


x=Eparticle
y=Eradio

x_err=x*sigma_e#*d_ratio
y_err=y*sigma_e_radio#*np.sqrt(p_ratio)
#x=x[sigma_e<0.8]
#y=y[sigma_e<0.8]
#x_err=x_err[sigma_e<0.8]
#y_err=y_err[sigma_e<0.8]


y_label='LOFAR radio-based energy,  (eV)'
#y_label='fluence-based energy, E$_{flu.}$ (eV)'
x_label='LORA particle-based energy,(eV)'


## plot radio energy comparison energy

'''
fig = plt.figure(figsize=(5,5))
ax1 = fig.add_subplot(1,1,1,aspect=1)
ax1.plot(x_fit,x_fit,color='black')
ax1.errorbar(x, y, xerr=x_err,yerr=y_err, marker='o',elinewidth=0.7,linestyle='none',markeredgecolor='black',ecolor='black',label='LOFAR events')

#ax1.plot(plt_x_line,plt_y_line,color='purple',linewidth=3,label='LOFAR prediction')

ax1.set_yscale('log')
ax1.set_xscale('log')

ax1.legend(numpoints=1,loc='upper left',framealpha=1,facecolor='white')

ax1.set_xlabel(x_label,fontsize=12)
ax1.set_ylabel(y_label,fontsize=12)

ax1.axis([2e16,3e18,2e16,3e18])

#ax1.axis([3e16,3e18,2e16,3e18])
#ax1.axis([2e16,3e18,7e3,1e8])

plt.tight_layout()
plt.grid()
plt.show()
raw_input()
plt.close()
'''









x=Eparticle
y=1.79*1e7*np.power(em_energy/1e18,2.006)#Srd
x_err=x*sigma_e#*d_ratio
y_err=y*sigma_e_radio#*np.sqrt(p_ratio)
x=x[sigma_e<0.8]
y=y[sigma_e<0.8]
x_err=x_err[sigma_e<0.8]
y_err=2*y_err[sigma_e<0.8]

err_use=np.sqrt(x_err*x_err+x_err*x_err)

print len(x)
print len(y)
print len(x_err)
print len(y_err)

################### FIT LOFAR

a0=1.3
b0=2.025
max_fit=7e19
min_fit=1.5e14
#res_data_comp=minimize(tools.e,[a0,b0],args=(x[(x>min_fit)*(x<max_fit)],y[(x>min_fit)*(x<max_fit)],x_err[(x>min_fit)*(x<max_fit)]),method='Nelder-Mead')
res_data_comp=minimize(tools.e_one,[a0,],args=(x[(x>min_fit)*(x<max_fit)],y[(x>min_fit)*(x<max_fit)],err_use[(x>min_fit)*(x<max_fit)]),method='Nelder-Mead')

a_lofar=res_data_comp['x'][0]
b_lofar=2#res_data_comp['x'][1]

print a_lofar
print b_lofar



###################



s_recon=np.zeros(len(y))

for i in np.arange(len(y)):
    s_recon[i]=tools.energy_to_rad(y[i],b_lofar,a_lofar)


y1=y
x1=x
y2=s_recon
sig1=y_err

chi2=tools.chi2(np.log10(y1),np.log10(y2),np.log10(sig1))/(len(y1))
print '____________________'
print 'chi2: {0}'.format(chi2)
print '____________________'

###################


#x_err=x*sigma_e#*d_ratio
#y_err=y*2*sigma_e_radio#*p_ratio


#x=x[sigma_e<0.8]
#y=y[sigma_e<0.8]
#x_err=x_err[sigma_e<0.8]
#y_err=y_err[sigma_e<0.8]

#x_label='direct scaling energy, $E_{dir. scal.}$ (eV)'

#y_label='direct scaling energy (eV)'
y_label='radiation energy, S$_{RD}$ (eV) / S$_{RD, proton}$'
y_label_log='LOFAR radiation energy, log(S$_{RD}$/eV)'


#y_label='fluence-based energy, E$_{flu.}$ (eV)'
x_label_log='LORA cosmic-ray energy, log( E$_{part.}$/eV)'
x_label='LORA cosmic-ray energy, E$_{part.}$ (eV)'

#hist_label='2 (E$_{flu.}$-E$_{dir. scal.}$)/(E$_{flu.}$+E$_{dir. scal.}$)'
hist_label='2 (E$_{part.}$-E$_{lofar\ fit}$)/(E$_{part.}$+E$_{lofar\ fit}$)'


fig = plt.figure(figsize=(5,5))
ax1 = fig.add_subplot(1,1,1)#,aspect=1)
#ax1.plot(x_fit,rad_p_fit,color='blue',linestyle=':',label='CORSIKA/CoREAS prediction, P')
#ax1.plot(x_fit,rad_fe_fit,color='red',linestyle=':',label='CORSIKA/CoREAS prediction, Fe')
ax1.fill_between(x_fit,tools.pred_sys_up(x_fit),tools.pred_sys_down(x_fit),alpha=0.3,color='g')
ax1.plot(x_fit,tools.prediction(x_fit),label="AERA PRL prediction",c='g',linewidth=3)
ax1.plot(x_fit,tools.lofar_prediction(x_fit,a_lofar,b_lofar),label="LOFAR prediction",c='purple',linewidth=3)
#ax1.fill_between(x_fit,tools.lofar_prediction(x_fit,a_lofar,b_lofar)*0.73,tools.lofar_prediction(x_fit,a_lofar,b_lofar)*1.27,facecolor='none',hatch="\\\\",edgecolor='purple')

ax1.fill_between(x_fit,tools.lofar_prediction(x_fit,a_lofar,b_lofar)*0.73,tools.lofar_prediction(x_fit,a_lofar,b_lofar)*1.27,alpha=0.3,color='purple')
#ax1.fill_between(x_fit,tools.pred_sys_up(x_fit),tools.pred_sys_down(x_fit),alpha=0.3,color='g')
#ax1.plot(x_fit,tools.prediction(x_fit),label="AERA PRL prediction",c='g',linewidth=3)


ax1.errorbar(x, y, xerr=x_err,yerr=y_err, marker='o',elinewidth=0.7,linestyle='none',markeredgecolor='black',ecolor='black',label='LOFAR events')
#ax1.plot(Eparticle[(xmax<800)*(sigma_r<7)],Srd[(xmax<800)*(sigma_r<7)],'x',color='red')
#ax1.plot(Eparticle[(xmax>550)*(xmax<800)],Srd[(xmax<800)*(xmax>550)],'x',color='red')



#ax1.plot(plt_x_line,plt_y_line,color='purple',linewidth=3,label='LOFAR prediction')

ax1.set_yscale('log')
ax1.set_xscale('log')

ax1.legend(numpoints=1,loc='upper left',framealpha=1,facecolor='white')

ax1.set_xlabel(x_label,fontsize=12)
ax1.set_ylabel(y_label,fontsize=12)

ax1.axis([2e16,3e18,1e4,2e8])

#ax1.axis([3e16,3e18,2e16,3e18])
#ax1.axis([2e16,3e18,7e3,1e8])

plt.tight_layout()
plt.grid()
plt.show()
raw_input()
plt.close()



fig = plt.figure(figsize=(5,5))
ax1 = fig.add_subplot(1,1,1)#,aspect=1)

ax1.fill_between(x_fit,tools.pred_sys_up(x_fit)/rad_p_fit,tools.pred_sys_down(x_fit)/rad_p_fit,alpha=0.3,color='g')
ax1.plot(x_fit,tools.prediction(x_fit)/rad_p_fit,label="AERA PRL prediction",c='g',linewidth=3)
ax1.plot(x_fit,tools.lofar_prediction(x_fit,a_lofar,b_lofar)/rad_p_fit,label="LOFAR prediction",c='purple',linewidth=3)
#ax1.fill_between(x_fit,tools.lofar_prediction(x_fit,a_lofar,b_lofar)*0.73,tools.lofar_prediction(x_fit,a_lofar,b_lofar)*1.27,facecolor='none',hatch="\\\\",edgecolor='purple')
ax1.fill_between(x_fit,tools.lofar_prediction(x_fit,a_lofar,b_lofar)*0.73/rad_p_fit,tools.lofar_prediction(x_fit,a_lofar,b_lofar)*1.27/rad_p_fit,alpha=0.3,color='purple')

ax1.plot(x_fit,rad_p_fit/rad_p_fit,color='blue',linestyle=':',linewidth=3,label='CORSIKA/CoREAS prediction, P')
ax1.plot(x_fit,rad_fe_fit/rad_p_fit,color='red',linestyle=':',linewidth=3,label='CORSIKA/CoREAS prediction, Fe')

#ax1.set_yscale('log')
ax1.set_xscale('log')

ax1.legend(numpoints=1,loc='upper right',framealpha=1,facecolor='white')

ax1.set_xlabel(x_label,fontsize=12)
ax1.set_ylabel(y_label,fontsize=12)

#ax1.axis([2e16,3e18,1e4,2e8])

#ax1.axis([3e16,3e18,2e16,3e18])
#ax1.axis([2e16,3e18,7e3,1e8])

plt.tight_layout()
plt.grid()
plt.show()
raw_input()
plt.close()








'''
fig = plt.figure(figsize=(5,5))
ax1 = fig.add_subplot(1,1,1)#,aspect=1)
#ax1.plot(x_fit,rad_p_fit,color='blue',linestyle=':',label='CORSIKA/CoREAS prediction, P')
#ax1.plot(x_fit,rad_fe_fit,color='red',linestyle=':',label='CORSIKA/CoREAS prediction, Fe')
#ax1.fill_between(x_fit,tools.pred_sys_up(x_fit),tools.pred_sys_down(x_fit),alpha=0.3,color='g')
#ax1.plot(x_fit,tools.prediction(x_fit),label="AERA PRL prediction",c='g',linewidth=3)
ax1.plot(x_fit,tools.lofar_prediction(x_fit,a_lofar,b_lofar),label="LOFAR prediction",c='purple',linewidth=3)


sc = ax1.scatter(Eparticle,Srd, alpha=1,s=5,vmin=np.min(sigma_e),vmax=0.5,c=sigma_e)
cbar = plt.colorbar(sc)
cbar.ax.set_ylabel(r'$\sigma_{E radio}$', rotation=90)

ax1.set_yscale('log')
ax1.set_xscale('log')

ax1.legend(numpoints=1,loc='upper left',framealpha=1,facecolor='white')

ax1.set_xlabel(x_label,fontsize=12)
ax1.set_ylabel(y_label,fontsize=12)

ax1.axis([2e16,3e18,1e4,2e8])

#ax1.axis([3e16,3e18,2e16,3e18])
#ax1.axis([2e16,3e18,7e3,1e8])

plt.tight_layout()
plt.grid()
plt.show()
raw_input()
plt.close()
'''




# contour
'''
fig = plt.figure(figsize=(5,5))
ax1 = fig.add_subplot(1,1,1)#,aspect=1)
ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))

ax1.axhline(y=4.0, color='grey', linestyle='-',linewidth=0.5)
ax1.axhline(y=5.0, color='grey', linestyle='-',linewidth=0.5)
ax1.axhline(y=6.0, color='grey', linestyle='-',linewidth=0.5)
ax1.axhline(y=7.0, color='grey', linestyle='-',linewidth=0.5)
ax1.axhline(y=8.0, color='grey', linestyle='-',linewidth=0.5)

ax1.axvline(x=16.0, color='grey', linestyle='-',linewidth=0.5)
ax1.axvline(x=17.0, color='grey', linestyle='-',linewidth=0.5)
ax1.axvline(x=18.0, color='grey', linestyle='-',linewidth=0.5)


cmap='gnuplot2_r'
ax1.hist2d(np.log10(x), np.log10(y), bins=(40, 40),cmap=cmap)#plt.cm.jet)

ax1.axis([np.log10(2e16),np.log10(3e18),np.log10(1e4),np.log10(2e8)])
ax1.fill_between(np.log10(x_fit),np.log10(tools.pred_sys_up(x_fit)),np.log10(tools.pred_sys_down(x_fit)),alpha=0.3,color='g')

#ax1.fill_between(np.log10(plt_x_line),np.log10(plt_y_line-plt_y_line*0.6),np.log10(plt_y_line+plt_y_line*0.6),alpha=0.3,color='purple')

ax1.plot(np.log10(x_fit),np.log10(tools.prediction(x_fit)),label="AERA PRL prediction",c='g',linewidth=3)
ax1.plot(np.log10(plt_x_line),np.log10(plt_y_line),color='purple',linewidth=3,label='LOFAR prediction')

ax1.plot(np.log10(x_fit),np.log10(rad_p_fit),color='blue',linestyle=':',label='CORSIKA/CoREAS prediction, P')
ax1.plot(np.log10(x_fit),np.log10(rad_fe_fit),color='red',linestyle=':',label='CORSIKA/CoREAS prediction, Fe')

ax1.legend(numpoints=1,loc='upper left',framealpha=1,facecolor='white')

ax1.set_xlabel(x_label_log,fontsize=12)
ax1.set_ylabel(y_label_log,fontsize=12)


plt.show()
raw_input()
plt.close()
'''



'''

######## hist
hist_label='2 (S$_{RD,AERA}$-S$_{RD,LOFAR}$)/(E$_{RD,AERA}$+E$_{RD,LOFAR}$)'
hist_label='2 (E$_{part.}$-E$_{radio}$)/(E$_{part.}$+E$_{radio}$)'

z=1.069*1e7*np.power(em_energy/1e18,1.994)#Eradio#1.79*1e7*np.power(em_energy/1e18,2.006)#Srd#Eradio
x=1.58*1e7*np.power(Eparticle/1e18,1.98)#Srd#EparticleEparticle#1.58*1e7*np.power(Eparticle/1e18,1.98)#Srd#Eparticle

diff=2*(x-z)/(x+z)
bin_start=-2.0
bin_stop=4.0
binWidth=0.05
#nBins=50.0
nBins=int((bin_stop-bin_start)/binWidth)
bins=np.arange(bin_start,bin_stop,binWidth)

p0,p1,p2, x_hist, y_hist=tools.fit_hist(diff,bin_start,bin_stop,int(nBins),-100)


fig = plt.figure(figsize=(4,4))
ax1 = fig.add_subplot(1,1,1)

ax1.hist(diff,bins=bins,alpha=0.5)
ax1.plot(x_hist,y_hist,color='green',label='mean {0:.2f}\n spread {1:.2f}'.format(p1,p2))

ax1.set_xlabel(hist_label,fontsize=12)

ax1.legend(numpoints=1,loc='upper left')
ax1.set_xlim([-1.0, 1.0])

plt.tight_layout()
plt.show()
raw_input()
plt.close()
'''
'''
bin_start=0.0
bin_stop=1.0

binWidthR=0.007
nBinsR=int((bin_stop-bin_start)/binWidthR)
binsR=np.arange(bin_start,bin_stop,binWidthR)

binWidthP=0.02
nBinsP=int((bin_stop-bin_start)/binWidthP)
binsP=np.arange(bin_start,bin_stop,binWidthP)

p0_r,p1_r,p2_r, x_hist_r, y_hist_r=tools.fit_hist(sigma_e_radio,bin_start,bin_stop,int(nBinsR),-0)
p0_p,p1_p,p2_p, x_hist_p, y_hist_p=tools.fit_hist_land(sigma_e,bin_start,bin_stop,int(nBinsP),0)
'''
#bins=np.arange(0,1,0.01)
'''
fig = plt.figure()
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)

ax1.hist(sigma_e_radio,alpha=0.7,bins=binsR)
ax1.plot(x_hist_r, y_hist_r,color='green',linewidth=2,label='mean {0:.2f}\n spread {1:.2f}'.format(p1_r,p2_r))

ax2.hist(sigma_e,alpha=0.7,bins=binsP)
ax2.plot(x_hist_p, y_hist_p,color='green',linewidth=2,label='mean {0:.2f}\n spread {1:.2f}'.format(p1_p,p2_p))

ax1.set_xlabel('std($f_r$)',fontsize=12)
ax2.set_xlabel('std($f_p$)',fontsize=12)


ax1.legend(numpoints=1,loc='upper right')
ax2.legend(numpoints=1,loc='upper right')

ax1.set_xlim([0.0,0.2])
ax2.set_xlim([0.0,0.5])

#ax1.set_yscale('log')
#ax2.set_yscale('log')


plt.show()
raw_input()
plt.close()
'''

'''
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
line=np.arange(0,5,1)
ax1.plot(line,line)
ax1.plot(np.sqrt(p_ratio),d_ratio*(particle_scale),'.')
plt.show()
raw_input()
plt.close()

'''


'''

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
line=np.arange(0,5,1)
#ax1.plot(np.degrees(zenith),sigma_e_radio,'.')
ax1.plot(sigma_e,sigma_e_radio,'.')

ax1.set_xlabel(r'$\sigma_{particle}$',fontsize=12)

ax1.set_ylabel(r'$\sigma_{radio}$',fontsize=12)


plt.show()
raw_input()
plt.close()

'''
