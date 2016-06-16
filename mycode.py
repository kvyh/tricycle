## eventually should be changed to be object-oriented
import numpy as np
import data
import binaries
import matplotlib.pyplot as plt
import lightcurve
from interpacf import interpolated_acf, dominant_period

def autocorrelate(allkics, write = False, append = ''):
    '''
    This uses the autocorrelation function interpacf
    '''
    period_dict = {}
    for n, kic in enumerate(allkics):
        binary = binaries.RealBinary(kic)
        time, flux, err, quarter = binary.curve_cut()
        flux -= np.mean(flux)
        lag, acf = interpolated_acf(time, flux)
        detected_period, strength = dominant_period(lag, acf)
        period_dict['kic'+str(kic)] = [detected_period, strength, binary.p_orb]
        
    if write:
        kics = []
        best_period = []
        strength = []
        p_orb = []
        for i,k in period_dict.iteritems():
            kics.append(i)
            best_period.append(k[0])
            strength.append(k[1])
            p_orb.append(k[2])
        array = np.array([[kics[0], best_period[0], strength[0], p_orb[0]]])
        for n in range(1,len(kics)):
            
            array = np.append(array,[[kics[n], best_period[n], strength[n], p_orb[n]]],axis=0)
        np.savetxt('periods'+append, array, fmt ='%.18s')
 
    print 'done with autocorrelation'        
    return period_dict

def period_power_dict(allkics, write = False, append = ''):
    period_power = {}
    best_periods = {}
    for n,kic in enumerate(allkics):
        binary=binaries.RealBinary(kic)

        period, power, best_period = binary.periodogram()
        period_power['kic'+str(kic)] = [period, power, binary.p_orb]
        where = np.abs(period - best_period).argmin()
        strength = power[where]
        best_periods['kic'+str(kic)] = [best_period, strength, binary.p_orb]
        
    #for printing
    if write:
        kics = []
        best_period = []
        strength = []
        p_orb = []
        for i,k in best_periods.iteritems():
            kics.append(i)
            best_period.append(k[0])
            strength.append(k[1])
            p_orb.append(k[2])
        array = np.array([[kics[0], best_period[0], strength[0], p_orb[0]]])
        for n in range(1,len(kics)):
            
            array = np.append(array,[[kics[n], best_period[n], strength[n], p_orb[n]]],axis=0)
        np.savetxt('periods'+append, array, fmt ='%.18s')
 
    print 'done with periodograms'
    return period_power, best_periods
    

def strongest_periods(period_power_dict):
    strongest_periods = {}
    for i,k in period_power_dict.iteritems():
        #this goes over each [period],[power]
        period = np.array(k[0])
        power = np.array(k[1])
        p_orb = np.array(k[2])
        period_mask = []
        for j in range(len(period)):
            if period[j]<45.0 and period[j]>.025:
                period_mask.append(j)
        index = np.where(power == max(power[period_mask]))
        strongest_periods[str(i)]=(float(period[index[0]]), float(power[index[0]]), float(p_orb))
    return strongest_periods

def pullfile(filename = 'best_periods'):
    f = open(filename)
    lines = f.read().split('\n')
    lines.pop()

    for i,line in enumerate(lines):
        lines[i]=line.split(' ')
    f.close()  

    kic = []
    best_period = []
    strength = []
    p_orb = []
    for t in range(len(lines)):
        kic.append(lines[t][0])
        best_period.append(float(lines[t][1]))
        strength.append(float(lines[t][2]))
        p_orb.append(float(lines[t][3]))
    return kic, best_period, strength, p_orb

def plot_EBs(xlim = (0,200), ylim = (0,10), dictionary = None, inputs = False ,filename = 'periodogram_all', minstrength = 0, maxstrength = 1):
    '''
    plots the ratio of p_orbital/p_rotational on the y-axis, and p_orbital on the x-axis
    '''

        
    if inputs:
        kic = []
        best_period = []
        strength = []
        p_orb = []
        for kics, k in dictionary.iteritems():
            kic.append(kics)
            best_period.append(k[0])
            strength.append(k[1])
            p_orb.append(k[2])
            
    else:
        kic, best_period, strength, p_orb = pullfile(filename)
    
    
    p_orb_rot = []
    for i in range(len(kic)):

        p_orb_rot.append(p_orb[i]/best_period[i])

    mask = np.where(np.array(strength) >= minstrength)
    p_orb_rot = np.array(p_orb_rot)[mask]
    p_orb = np.array(p_orb)[mask]
    strength = np.array(strength)[mask]
    
    mask2 = np.where(strength <= maxstrength)
    p_orb_rot = p_orb_rot[mask2]
    p_orb = p_orb[mask2]
    strength = strength[mask2]
    

    plt.plot(p_orb, p_orb_rot, linestyle = 'none', marker='o', markersize=1)
    if inputs:
        plt.title('Plotted from' + str(dictionary))
        #this may need modification
    else:
        plt.title('Plotted from ' + filename + '\nshowing EBs with strength between ' + str(minstrength) + ' and ' + str(maxstrength)) 
    plt.ylabel('orbital period/rotation period')
    plt.xlabel('orbital period')
    plt.ylim(ylim)
    plt.xlim(xlim)
    plt.show()

#allkics = data.select_kics()


'''
for n in range(len(allkics)/100):
    autocorrelate(allkics[100*n:100+100*n], write = True, append = '_autocor_'+ str(n))
    period_power_dict(allkics[100*n:100+100*n], write = True, append = '_periodogram_' + str(n))
autocorrelate(allkics[2800:], write = True, append = '_autocor_last')
period_power_dict(allkics[2800:], write = True, append = '_periodogram_last')
'''

#plot_EBs(xlim = (0,35), ylim = (0,5), filename = 'periodogram_all', minstrength = .3)
plot_EBs(xlim = (0,35), ylim = (0,5), filename = 'autocor_all', minstrength = .7) 
        
###