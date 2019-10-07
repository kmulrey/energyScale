import numpy as np
from optparse import OptionParser
import os
import cPickle
import os.path
import re
from glob import glob


event_path='/Users/kmulrey/radio/events/'
protonfile=open('proton_event_list.txt','w')
ironfile=open('iron_event_list.txt','w')



events = os.listdir(event_path)

def write_file_list(dirIron,dirProton,filelistProton,filelistIron):
    if os.path.isdir(dirIron):
        for file in os.listdir(dirIron):
            if file.endswith(".long"):
                longfile=dirIron+file
                runnr=file.split(".")[0]
                runnr=runnr.split("T")[1]
                filelistIron.append(runnr)
                ironfile.write('{0} {1}\n'.format(dirIron,runnr))
    
    if os.path.isdir(dirProton):

        for file in os.listdir(dirProton):
        
            if file.endswith(".long"):
                longfile=dirProton+file
                runnr=file.split(".")[0]
                runnr=runnr.split("T")[1]
                filelistProton.append(runnr)
                protonfile.write('{0} {1}\n'.format(dirProton,runnr))



for e in np.arange(len(events)):
    event=events[e]

    try:
        dirIron=event_path+event+'/2/coreas/iron/'
        dirProton=event_path+event+'/2/coreas/proton/'
        write_file_list(dirIron,dirProton,filelistProton,filelistIron)
    except:
        #print 'no 2 directory'
        pass
    try:
        dirIron=event_path+event+'/1/coreas/iron/'
        dirProton=event_path+event+'/1/coreas/proton/'
        write_file_list(dirIron,dirProton,filelistProton,filelistIron)
    except:
        #print 'no 1 directory'
        pass
    
    try:
        dirIron=event_path+event+'/0/coreas/iron/'
        dirProton=event_path+'/0/coreas/proton/'
        write_file_list(dirIron,dirProton,filelistProton,filelistIron)
    except:
        #print 'no 0 directory'
        pass



