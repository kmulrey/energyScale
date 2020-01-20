import numpy as np
from optparse import OptionParser
import os
import cPickle
import os.path
import matplotlib.pyplot as plt
import scipy.interpolate as intp
from scipy.interpolate import interp1d
from matplotlib import cm
import re
from scipy import signal
from scipy.signal import hilbert
import sys

sys.path.insert(0, '../simulation')

import fluence as flu
import radiation_energy as rad
import sim_functions as sim
import helper as helper
import debug_sim as debug_sim


plt.ion()


#event_path='/Users/kmulrey/radio/events/'
#event_path='/vol/astro7/lofar/sim/pipeline/events/'
event_path='/vol/astro3/lofar/sim/pipeline/events/'

def integrate(r,flu0,flu1):
    n=len(r)
    dr=r[1]-r[0]
    integral=0
    for i in np.arange(n-1):
        r0=r[i]
        r1=r[i+1]
        val0=r0*(flu0[i]+flu1[i])
        val1=r1*(flu0[i+1]+flu1[i+1])
        integral=integral+(val0+val1)*0.5*dr
    
    
    
    return 2*np.pi*integral*6.2415e18 # to eV

def integrate_one_pol(r,flu):
    n=len(r)
    dr=r[1]-r[0]
    integral=0
    for i in np.arange(n-1):
        r0=r[i]
        r1=r[i+1]
        val0=r0*(flu[i])
        val1=r1*(flu[i+1])
        integral=integral+(val0+val1)*0.5*dr
    
    
    
    return 2*np.pi*integral*6.2415e18 # to eV


parser = OptionParser()
parser.add_option('-a', '--start',type='int',help='line number of event to start',default=0)
#parser.add_option('-b', '--stop',type='int',help='line number of event to stop',default=0)
parser.add_option('-t', '--type',type='string',help='type',default='proton')
#parser.add_option('-o', '--output',type='string',help='output',default='event_information.dat')

(options, args) = parser.parse_args()
a=options.start
#b=options.stop
#outfilename=options.output
type=options.type


#outfilename='/Users/kmulrey/LOFAR/energy/energyScale/data/radio/sim_local/output/'+type+'_'+str(a)+'_debug2.dat'
outfilename='output_US_standard/'+type+'_'+str(a)+'.dat'
# info to record

em_energy_list=[]
energy_list=[]
zenith_list=[]
az_list=[]
xmax_list=[]
alpha_list=[]
Erad_list=[]
Erad_gm_list=[]
Erad_ce_list=[]
rho_list=[]
rho2_list=[]

event_list=[]

#filename=type+'_event_list.txt'
filename=type+'_runs.txt'

sim_dir=[]
runnr=[]
count=0


with open(filename) as f:
    for line in f:
        sim_dir.append(line.split()[0])
        
        runnr.append(line.split()[1])
        count=count+1


for i in np.arange(a,a+1):
    #try:
    for r in np.arange(1):

        event= sim_dir[i].split('events/')[1].split('/')[0]
        print '{0}  {1}'.format(sim_dir[i], runnr[i])
        
        ant_pos,times,efield,zenith,az,energy,xmax=debug_sim.get_efield(sim_dir[i], runnr[i])
        em_energy,other_energy_hold,total_energy_hold=sim.getEM(sim_dir[i], runnr[i])
        alpha=sim.GetAlpha(zenith,az,1.1837)
        atm=sim.get_US_atm(event)
  
        hi=debug_sim.get_vertical_height(xmax,atm)   # height in cm
        at=debug_sim.get_atmosphere(hi,atm)
        rho=debug_sim.return_density(hi, atm)
        rho2=debug_sim.return_density(hi/np.cos(zenith), atm)

        dmax=helper.get_distance_xmax_geometric(zenith, xmax, atm)


        ant_pos_shower=sim.ground_to_shower(ant_pos,zenith,az,B=1.1837)
        e_filt,time_filt=debug_sim.lofar_filter(efield,times)
        fluence=debug_sim.calculate_energy_fluence_vector(e_filt,time_filt*1e-9, signal_window=100., remove_noise=True)


        pos_uvw_vxb=ant_pos_shower[0::8]
        pos_uvw_vxvxb=ant_pos_shower[2::8]
        neg_uvw_vxvxb=ant_pos_shower[6::8]

        fluence_vxvxb_0=np.concatenate([fluence[2::8].T[0],fluence[6::8].T[0]])
        fluence_vxvxb_1=np.concatenate([fluence[2::8].T[1],fluence[6::8].T[1]])

        pos_vxvxb_all=np.concatenate([neg_uvw_vxvxb.T[1],pos_uvw_vxvxb.T[1]])

        inds = pos_vxvxb_all.argsort()

        sorted_pos=pos_vxvxb_all[inds]
        sorted_fluence_gm=fluence_vxvxb_0[inds]
        sorted_fluence_ce=fluence_vxvxb_1[inds]
        xnew = np.linspace(0, 400, num=1000, endpoint=True)
        f0 = interp1d(sorted_pos, sorted_fluence_gm, kind='cubic')
        f1 = interp1d(sorted_pos, sorted_fluence_ce, kind='cubic')

        Erad=integrate(xnew,f0(xnew),f1(xnew))
        Erad_gm=integrate_one_pol(xnew,f0(xnew))
        Erad_ce=integrate_one_pol(xnew,f1(xnew))



        em_energy_list.append(em_energy*1e9)
        energy_list.append(energy*1e9)
        zenith_list.append(zenith)
        az_list.append(az)
        xmax_list.append(xmax)
        alpha_list.append(alpha)
        Erad_list.append(Erad)
        Erad_gm_list.append(Erad_gm)
        Erad_ce_list.append(Erad_ce)
        rho_list.append(rho)
        rho2_list.append(rho2)

        event_list.append(int(event))

info={'em_energy':np.asarray(em_energy_list),
    'total_energy':np.asarray(energy_list),
    'zenith':np.asarray(zenith_list),
    'azimuth':np.asarray(az_list),
    'energy':np.asarray(energy_list),
    'xmax':np.asarray(xmax_list),
    'alpha':np.asarray(alpha_list),
    'Erad':np.asarray(Erad_list),
    'Erad_gm':np.asarray(Erad_gm_list),
    'Erad_ce':np.asarray(Erad_ce_list),
    'density':np.asarray(rho_list),
    'density2':np.asarray(rho2_list),

    'event':np.asarray(event_list),

}

outfile=open(outfilename,'w')
cPickle.dump(info,outfile)
outfile.close()

