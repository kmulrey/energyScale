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
import sys
sys.path.insert(0, '../simulation')


import fluence as flu
import radiation_energy as rad
import sim_functions as sim
import helper as helper
import sim_functions as sim
sys.path.insert(0, '../sim_local')

import debug_sim as debug_sim


#import slope_functions as slope
#import density

freq_a=30.0
freq_b=77.0
nAnt=160
go_ahead=1

n_spoke=8
n_rad=20

dist=np.arange(25,n_rad*25+25,25)




parser = OptionParser()
parser.add_option('-e', '--event',type='int',help='line number of event',default=0)

(options, args) = parser.parse_args()
e=int(options.event)


#reco_dir='/vol/astro7/lofar/sim/pipeline/production_analysis_oct2018/'
reco_dir='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/results/production_analysis_radio_only_FINAL/'
event_path='/vol/astro7/lofar/sim/pipeline/events/'

output_dir='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/energyScale/radio/lofar_events/output/'


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

def return_file_list(dirIron,dirProton,filelist,xmax_list):
    
    for file in os.listdir(dirIron):
        if file.endswith(".long"):
            longfile=dirIron+file
            runnr=file.split(".")[0]
            runnr=runnr.split("T")[1]
            sim_path=dirIron+'SIM'+runnr+'_coreas/'
            hillas = np.genfromtxt(re.findall("PARAMETERS.*",open(longfile,'r').read()))[2:]
            #print 'iron: {0}      {1}   {2}'.format(file,hillas[2],sim_path)
            filelist.append(sim_path)
            xmax_list.append(hillas[2])
    
    
    for file in os.listdir(dirProton):
        if file.endswith(".long"):
            longfile=dirProton+file
            runnr=file.split(".")[0]
            runnr=runnr.split("T")[1]
            sim_path=dirProton+'SIM'+runnr+'_coreas/'
            #print sim_path
            hillas = np.genfromtxt(re.findall("PARAMETERS.*",open(longfile,'r').read()))[2:]
            #print 'proton: {0}      {1}  {2}'.format(file,hillas[2],sim_path)
            filelist.append(sim_path)
            xmax_list.append(hillas[2])

    return filelist,xmax_list

def find_SIM_dir(event_no,xmax):
    
    event=str(int(event_no))
    filelist=[]
    xmax_list=[]
    
    try:
        dirIron=event_path+event+'/2/coreas/iron/'
        dirProton=event_path+event+'/2/coreas/proton/'
        filelist,xmax_list=return_file_list(dirIron,dirProton,filelist,xmax_list)
    except:
        #print 'no 2 directory'
        pass
    try:
        dirIron=event_path+event+'/1/coreas/iron/'
        dirProton=event_path+event+'/1/coreas/proton/'
        filelist,xmax_list=return_file_list(dirIron,dirProton,filelist,xmax_list)
    except:
        #print 'no 1 directory'
        pass
    
    try:
        dirIron=event_path+event+'/0/coreas/iron/'
        dirProton=event_path+event+'/0/coreas/proton/'
        #print dirProton
        filelist,xmax_list=return_file_list(dirIron,dirProton,filelist,xmax_list)
    except:
        #print 'no 0 directory'
        pass

    closest_xmax=100000
    closest_index=100000
    min_diff=10000
    

    for i in np.arange(len(xmax_list)):
        min_diff_here=np.abs(xmax_list[i]-xmax)
        if min_diff_here<min_diff:
            closest_xmax=xmax_list[i]
            closest_index=i
            min_diff=min_diff_here

    try:
        return_file=filelist[closest_index]
        return_xmax=xmax_list[closest_index]
    except:
        return_file='null'
        return_xmax=0

    return_dir=(return_file.split('_')[0]).split('SIM')[0]
    return_sim=(return_file.split('_')[0]).split('SIM')[1]

    return return_dir, return_sim,return_xmax



event_file=open('energy_events.txt','r')
event_list=np.genfromtxt(event_file,dtype='int')
event_file.close()

nEvents=len(event_list)

############# variables that i want to save ################



# fit parameters
core_x_list=[]
core_y_list=[]
x_off_list=[]
y_off_list=[]
events_list=[]
p_ratio_list=[]
d_ratio_list=[]
combchi2_list=[]
radiochi2_list=[]
xmax_fit_list=[]

# sim parameters
clip_ratio_list=[]
rho_list=[]
em_energy_list=[]
total_energy_list=[]
zenith_list=[]
azimuth_list=[]
energy_list=[]
xmax_list=[]
dmax_list=[]
alpha_list=[]
Erad_list=[]
Erad_gm_list=[]
Erad_ce_list=[]
type_list=[]

converge_list=[]
lora_dens_list=[]




#for i in np.arange(nEvents):
for i in np.arange(e,e+1):

    print '______________________'
    print event_list[i]

    my_file3 = reco_dir+'reco'+str(int(event_list[i]))+'_3.dat'
    my_file2 = reco_dir+'reco'+str(int(event_list[i]))+'_2.dat'
    my_file1 = reco_dir+'reco'+str(int(event_list[i]))+'_1.dat'
    my_file = reco_dir+'reco'+str(int(event_list[i]))+'.dat'



    if os.path.exists(my_file3)==True:
        use_file=my_file3
    elif os.path.exists(my_file2)==True:
        use_file=my_file2
    elif os.path.exists(my_file1)==True:
        use_file=my_file1
    elif os.path.exists(my_file)==True:
        use_file=my_file
    else:
        print 'no file for {0}'.format(int(event_list[i]))
        continue


    
    for g in np.arange(1):
        analysisinfo = cPickle.load(open(use_file,'r'))
        #print analysisinfo.keys()


        xmaxreco_temp=analysisinfo['xmaxreco']
        xmax_temp=analysisinfo['realxmax']
        xmaxbest_temp=analysisinfo['xmaxbest']

        core_x_temp=analysisinfo['core_x']+analysisinfo['xoff']
        core_y_temp=analysisinfo['core_y']+analysisinfo['yoff']
        x_off_temp=analysisinfo['xoff']
        y_off_temp=analysisinfo['yoff']
        p_ratio_temp=analysisinfo['p_ratio']
        d_ratio_temp=analysisinfo['d_ratio']
        combchi2_temp=analysisinfo['combchi2']
        radiochi2_temp=analysisinfo['radiochi2']
        energy_temp=analysisinfo['energy']
        zenith_temp=analysisinfo['zenith']
        azimuth_temp=analysisinfo['azimuth']
        converge_temp1=analysisinfo['fitconverged']
        lora_dens_temp=analysisinfo['lora_dens']
        converge_temp=0
        if converge_temp1==True:
            converge_temp=1

        print 'converged: {0}'.format(converge_temp)

        
        sim_dir,runnr,closest_xmax=find_SIM_dir(int(event_list[i]),xmaxreco_temp)
      #event= sim_dir.split('events/')[1].split('/')[0]
        
        if 'proton' in sim_dir:
            type=0
        elif 'iron' in sim_dir:
            type=1
        else:
            type=2

        print 'type: {0}'.format(type)


        event= sim_dir.split('events/')[1].split('/')[0]

        print sim_dir.split('events/')[1]
        print '{0}  {1}'.format(sim_dir, runnr)
    
        ant_pos,times,efield,zenith,az,energy,xmax=debug_sim.get_efield(sim_dir, runnr)
        em_energy,other_energy_hold,total_energy_hold=sim.getEM(sim_dir, runnr)
        alpha=sim.GetAlpha(zenith,az,1.1837)
        atm=sim.get_atm(event)
        
        hi=debug_sim.get_vertical_height(xmax,atm)   # height in cm
        at=debug_sim.get_atmosphere(hi,atm)
        rho=debug_sim.return_density(hi, atm)
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
        clip_ratio=rad.get_clipping(dmax)

        
        
        em_energy_list.append(em_energy*1e9)
        energy_list.append(energy*1e9)
        zenith_list.append(zenith)
        azimuth_list.append(az)
        xmax_list.append(xmax)
        alpha_list.append(alpha)
        Erad_list.append(Erad)
        Erad_gm_list.append(Erad_gm)
        Erad_ce_list.append(Erad_ce)
        rho_list.append(rho)
        events_list.append(int(event))
        clip_ratio_list.append(clip_ratio)
        type_list.append(type)
        dmax_list.append(dmax)


        xmax_fit_list.append(xmaxreco_temp)
        core_x_list.append(core_x_temp)
        core_y_list.append(core_y_temp)
        x_off_list.append(x_off_temp)
        y_off_list.append(y_off_temp)
        p_ratio_list.append(p_ratio_temp)
        d_ratio_list.append(d_ratio_temp)
        combchi2_list.append(combchi2_temp)
        radiochi2_list.append(radiochi2_temp)
        
        converge_list.append(converge_temp)
        lora_dens_list.append(len(lora_dens_temp))



    

#   except:
#       print'can\'t open file'


info={'em_energy':np.asarray(em_energy_list),
    'zenith':np.asarray(zenith_list),
    'azimuth':np.asarray(azimuth_list),
    'energy':np.asarray(energy_list),
    'xmax':np.asarray(xmax_list),
    'xmax_fit':np.asarray(xmax_fit_list),
    'core_x':np.asarray(core_x_list),
    'core_y':np.asarray(core_y_list),
    'x_off':np.asarray(x_off_list),
    'y_off':np.asarray(y_off_list),
    'p_ratio':np.asarray(p_ratio_list),
    'd_ratio':np.asarray(d_ratio_list),
    'combchi2':np.asarray(combchi2_list),
    'radiochi2':np.asarray(radiochi2_list),
    'event':np.asarray(events_list),
    'alpha':np.asarray(alpha_list),
    'Erad':np.asarray(Erad_list),
    'Erad_gm':np.asarray(Erad_gm_list),
    'Erad_ce':np.asarray(Erad_ce_list),

    'clip':np.asarray(clip_ratio_list),
    'density':np.asarray(rho_list),
    'dmax':np.asarray(dmax_list),
    'type':np.asarray(type_list),
    'converged':np.asarray(converge_list),
    'lora_density':np.asarray(lora_dens_list)

}

print 'event okay'

#outfile=open('output/event_'+str(int(event_list[e]))+'_information_oct8_2019.dat','w')
outfile=open(output_dir+'/event_'+str(int(event_list[e]))+'_information_feb17_2020.dat','w')

cPickle.dump(info,outfile)
outfile.close()


print 'file written'


