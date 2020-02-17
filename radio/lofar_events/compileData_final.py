import numpy as np
import os
import cPickle


#data_dir='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/radio/lofar_events/output/'
data_dir='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/energyScale/radio/lofar_events/output/'
#data_dir='output/'

file_list=[]

for file in os.listdir(data_dir):
    file_list.append(data_dir+file)

outfilename='compiled_sims_all_final_feb17.dat'








em_energy=np.empty([0])
energy=np.empty([0])
zenith=np.empty([0])
azimuth=np.empty([0])
xmax=np.empty([0])
xmax_fit=np.empty([0])
core_x=np.empty([0])
core_y=np.empty([0])
x_off=np.empty([0])
y_off=np.empty([0])
p_ratio=np.empty([0])
d_ratio=np.empty([0])
combchi2=np.empty([0])
radiochi2=np.empty([0])
event=np.empty([0])
alpha=np.empty([0])
Erad=np.empty([0])
Erad_ce=np.empty([0])
Erad_gm=np.empty([0])
clip=np.empty([0])
density=np.empty([0])
dmax=np.empty([0])
type=np.empty([0])
converged=np.empty([0])
lora_count=np.empty([0])


for i in np.arange(len(file_list)):
    if '2020' in file_list[i]:
        print file_list[i]
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
        Erad_ce_hold=info['Erad_ce']
        Erad_gm_hold=info['Erad_gm']
        clip_hold=info['clip']
        xmax_fit_hold=info['xmax_fit']
        core_x_hold=info['core_x']
        core_y_hold=info['core_y']
        x_off_hold=info['x_off']
        y_off_hold=info['y_off']
        p_ratio_hold=info['p_ratio']
        d_ratio_hold=info['d_ratio']
        combchi2_hold=info['combchi2']
        radiochi2_hold=info['radiochi2']
        event_hold=info['event']
        density_hold=info['density']
        dmax_hold=info['dmax']
        type_hold=info['type']
        converged_hold=info['converged']
        lora_count_hold=info['lora_density']
        
        
        
        
        #try:
        for i in np.arange(1):
            energy=np.concatenate((energy,energy_hold))
            em_energy=np.concatenate((em_energy,em_energy_hold))
            zenith=np.concatenate((zenith,zenith_hold))
            azimuth=np.concatenate((azimuth,azimuth_hold))
            xmax=np.concatenate((xmax,xmax_hold))
            alpha=np.concatenate((alpha,alpha_hold))
            Erad=np.concatenate((Erad,Erad_hold))
            Erad_ce=np.concatenate((Erad_ce,Erad_ce_hold))
            Erad_gm=np.concatenate((Erad_gm,Erad_gm_hold))

            clip=np.concatenate((clip,clip_hold))

            xmax_fit=np.concatenate((xmax_fit,xmax_fit_hold))
            core_x=np.concatenate((core_x,core_x_hold))
            core_y=np.concatenate((core_y,core_y_hold))
            x_off=np.concatenate((x_off,x_off_hold))
            y_off=np.concatenate((y_off,y_off_hold))
            p_ratio=np.concatenate((p_ratio,p_ratio_hold))
            d_ratio=np.concatenate((d_ratio,d_ratio_hold))
            combchi2=np.concatenate((combchi2,combchi2_hold))
            radiochi2=np.concatenate((radiochi2,radiochi2_hold))
            event=np.concatenate((event,event_hold))
            density=np.concatenate((density,density_hold))
            dmax=np.concatenate((dmax,dmax_hold))
            type=np.concatenate((type,type_hold))
            converged=np.concatenate((converged,converged_hold))
            lora_count=np.concatenate((lora_count,lora_count_hold))


        # except:
        #print 'oops'



info={'em_energy':np.asarray(em_energy),
    'energy':np.asarray(energy),
    'zenith':np.asarray(zenith),
    'azimuth':np.asarray(azimuth),
    'energy':np.asarray(energy),
    'xmax':np.asarray(xmax),
    'alpha':np.asarray(alpha),
    'Erad':np.asarray(Erad),
    'Erad_ce':np.asarray(Erad_ce),
    'Erad_gm':np.asarray(Erad_gm),
    'clip':np.asarray(clip),
    'xmax_fit':np.asarray(xmax_fit),
    'core_x':np.asarray(core_x),
    'core_y':np.asarray(core_y),
    'x_off':np.asarray(x_off),
    'y_off':np.asarray(y_off),
    'p_ratio':np.asarray(p_ratio),
    'd_ratio':np.asarray(d_ratio),
    'combchi2':np.asarray(combchi2),
    'radiochi2':np.asarray(radiochi2),
    'event':np.asarray(event),
    'density':np.asarray(density),
    'dmax':np.asarray(dmax),
    'type':np.asarray(type),
    'converged':np.asarray(converged),
    'lora_count':np.asarray(lora_count)

}



outfile=open(outfilename,'w')
cPickle.dump(info,outfile)
outfile.close()









