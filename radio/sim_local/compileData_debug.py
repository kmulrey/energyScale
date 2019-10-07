import numpy as np
import os
import cPickle


data_dir='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/radio/simulation/output/'
#data_dir='/Users/kmulrey/LOFAR/energy/energyScale/data/radio/sim_local/output/'
#out_dir='/Users/kmulrey/LOFAR/energy/energyScale/data/radio/sim_local/'
out_dir=''
type='proton'

file_list=[]

for file in os.listdir(data_dir):
    if file.startswith(type):
        if 'debug2' in file:
            file_list.append(data_dir+file)


outfilename='compiled_sims_'+type+'_debug2.dat'

em_energy=np.empty([0])
energy=np.empty([0])
zenith=np.empty([0])
azimuth=np.empty([0])
energy=np.empty([0])
xmax=np.empty([0])
alpha=np.empty([0])
Erad=np.empty([0])
Erad_gm=np.empty([0])
Erad_ce=np.empty([0])
density=np.empty([0])
density2=np.empty([0])
event=np.empty([0])


for i in np.arange(len(file_list)):
    infile=open(file_list[i],'r')
    info=cPickle.load(infile)
    infile.close()
 
    em_energy_hold=info['em_energy']
    energy_hold=info['energy']
    zenith_hold=info['zenith']
    azimuth_hold=info['azimuth']
    xmax_hold=info['xmax']
    alpha_hold=info['alpha']
    Erad_hold=info['Erad']
    Erad_gm_hold=info['Erad_gm']
    Erad_ce_hold=info['Erad_ce']
    density_hold=info['density']
    density2_hold=info['density2']
    event_hold=info['event']

    energy=np.concatenate((energy,energy_hold))
    em_energy=np.concatenate((em_energy,em_energy_hold))
    zenith=np.concatenate((zenith,zenith_hold))
    azimuth=np.concatenate((azimuth,azimuth_hold))
    xmax=np.concatenate((xmax,xmax_hold))
    alpha=np.concatenate((alpha,alpha_hold))
    Erad=np.concatenate((Erad,Erad_hold))
    Erad_gm=np.concatenate((Erad_gm,Erad_gm_hold))
    Erad_ce=np.concatenate((Erad_ce,Erad_ce_hold))
    density=np.concatenate((density,density_hold))
    density2=np.concatenate((density2,density2_hold))
    event=np.concatenate((event,event_hold))



info={'em_energy':np.asarray(em_energy),
    'energy':np.asarray(energy),
    'zenith':np.asarray(zenith),
    'azimuth':np.asarray(azimuth),
    'xmax':np.asarray(xmax),
    'alpha':np.asarray(alpha),
    'Erad':np.asarray(Erad),
    'Erad_gm':np.asarray(Erad_gm),
    'Erad_ce':np.asarray(Erad_ce),
    'density':np.asarray(density),
    'density2':np.asarray(density2),
    'event':np.asarray(event),

}



outfile=open(out_dir+outfilename,'w')
cPickle.dump(info,outfile)
outfile.close()










