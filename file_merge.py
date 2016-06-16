import numpy as np
import data
import binaries
import matplotlib.pyplot as plt


def file_merge(files, newfile):
    '''
    This takes multiple files with the format kic best_period orbital_period
    '''
    kic = []
    best_period = []
    strength = []
    p_orb = []
    for n in files:
        f = open(n)
        lines = f.read().split('\n')
        lines.pop()

        for i,line in enumerate(lines):
            lines[i]=line.split(' ')
        f.close()  

        for t in range(len(lines)):
            kic.append(lines[t][0])
            best_period.append(float(lines[t][1]))
            strength.append(float(lines[t][2]))
            p_orb.append(float(lines[t][3]))

            
    array = np.array([[kic[0], best_period[0], strength[0], p_orb[0]]])
    
    for n in range(1,len(kic)):
        array = np.append(array,[[kic[n], best_period[n], strength[n], p_orb[n]]],axis=0)
    np.savetxt(newfile, array, fmt ='%.18s')

  
