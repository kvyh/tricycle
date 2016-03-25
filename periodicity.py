import binaries
import numpy as np
import matplotlib.pyplot as plt
import data
from interpacf import interpolated_acf, dominant_period
from scipy.ndimage import gaussian_filter
from scipy import signal

#This code uses modules from https://github.com/StellarArmy/tricycle and
#https://github.com/bmorris3/interp-acf and requires the gatspy and MySQLdb
#packages to run properly

class Periodicity:
    '''
    Takes a list of kic numbers as inputs
    '''
    def __init__(self, kic_list):
        self.kics = kic_list
        
    def lightcurve(self, kic, curve_cut = True):
        binary = binaries.RealBinary(kic)
        
        time, flux, err, quarter = binary.curve_cut()
        plt.plot(time, flux)
        plt.show()
    
    def autocorrelate(self, kic):
        '''
        Takes a kic number and outputs the autocorrelation function and the best period
        '''
        binary = binaries.RealBinary(kic)
        time, flux, err, quarter = binary.curve_cut()
        flux -= np.mean(flux)
        lag, acf = interpolated_acf(time, flux)
        detected_period = binary.p_orb
        strength = 0
        try:
            detected_period, strength = dominant_period(lag, acf, min=0.2, max=45.0)
        except TypeError:
            print 'TypeError for kic ' + str(kic)
        return lag, acf, detected_period, strength, binary.p_orb
        
    
    def autocorrelate_period(self, kics, write = False, append = ''):
        '''
        Makes a dictionary mapping each kic number to the strongest period, the strength of the period
        and the orbital period of the binary. If desired, it will also write the dictionary to a file
        ending with append.
        '''
        period_dict = {}
        for n, kic in enumerate(self.kics):
            lag, acf, detected_period, strength, p_orb = self.autocorrelate(kic)
            period_dict['kic'+str(kic)] = [detected_period, strength, p_orb]

        if write:
            data = period_dict
            kics = [x for x in data.iterkeys()]
            best_period = [x[0] for x in data.itervalues()]
            strength = [x[1] for x in data.itervalues()]
            p_orb = [x[2] for x in data.itervalues()]
            
            array = np.array([[kics[0], best_period[0], strength[0], p_orb[0]]])
            for n in range(1,len(kics)):

                array = np.append(array,[[kics[n], best_period[n], strength[n], p_orb[n]]],axis=0)
            np.savetxt('periods'+append, array, fmt ='%.18s')

        print 'done with autocorrelation'        
        return period_dict
    
    def modified_lomb_scargle(self, kic):
        '''
        Takes a kic number and outputs the periodogram and the best period
        '''
        binary = binaries.RealBinary(kic)
        period, power, best_period = binary.periodogram(period_range = (.05,45))
        return period, power, best_period, binary.p_orb
    
    
    def lomb_scargle_period(self, kics, write = False, append = ''):
        '''
        Makes a dictionary mapping each kic number to the strongest period, the strength of the period
        and the orbital period of the binary. If desired, it will also write the dictionary to a file
        ending with append.
        '''
        period_power = {}
        best_periods = {}
        for n,kic in enumerate(self.kics):
            period, power, best_period, p_orb = self.modified_lomb_scargle(kic)
            period_power['kic'+str(kic)] = [period, power, p_orb]
            where = np.abs(period - best_period).argmin()
            strength = power[where]
            best_periods['kic'+str(kic)] = [best_period, strength, p_orb]

        #for printing
        if write:
            data = best_periods
            kics = [x for x in data.iterkeys()]
            best_period = [x[0] for x in data.itervalues()]
            strength = [x[1] for x in data.itervalues()]
            p_orb = [x[2] for x in data.itervalues()]
            
            array = np.array([[kics[0], best_period[0], strength[0], p_orb[0]]])
            for n in range(1,len(kics)):

                array = np.append(array,[[kics[n], best_period[n], strength[n], p_orb[n]]],axis=0)
            np.savetxt('periods'+append, array, fmt ='%.18s')

        print 'done with periodograms'
        return best_periods
    
    def plot_periodogram(self, kic, l_s = True, autocor = False):
        '''
        Takes the kic number and plots the periodogram and autocorrelation function
        '''
        detected_period =0
        best_period = 0
        if l_s:
            period, power, best_period, p_orb = self.modified_lomb_scargle(kic)
            plt.plot(period, power, color = 'g', label = 'Lomb_Scargle')
            plt.axvline(x = best_period, ls = '--',color = 'c', label = 'Primary L-S period')
            
        if autocor:
                     
            lag, acf, detected_period, strength, p_orb = self.autocorrelate(kic)
            plt.plot(lag, acf/np.max(acf), label='ACF')
            sigma = 18/2.355
            truncate = 56/sigma
            smooth_acf = gaussian_filter(acf, sigma, truncate=truncate)
            plt.plot(lag, smooth_acf/np.max(smooth_acf), label='Smoothed ACF')
            plt.axvline(x = detected_period, ls='--', color='r', label='Primary ACF period')

        plt.axvline(x = p_orb, label = 'Orbital Period')    
        plt.legend()
        
        plt.title('strength of periods for kic'+str(kic))
        plt.xlim(0, max(detected_period,best_period)*10)
        plt.ylabel('strength')
        plt.xlabel('period')
        plt.show()
        
        
        
        
    def pullfile(self, filename = 'best_periods'):
        '''
        This opens files that are written by the autocorrelate_period and lomb_scargle_period functions. 
        It returns a dictionary mapping each kic number to the strongest period, the strength of the period
        and the orbital period of the binary if used to open such a file.
        '''
        f = open(filename)
        lines = f.read().split('\n')
        lines.pop()

        for i,line in enumerate(lines):
            lines[i]=line.split(' ')
        f.close()  


        best_periods = {}
        for t in range(len(lines)):
            best_periods[lines[t][0]] = [float(lines[t][1]), float(lines[t][2]), float(lines[t][3])]
           
        return best_periods
    
    def refined_periods(self, dict1, dict2, correlation = (.95,1.05)):
        '''
        Takes the dictionaries output by lomb_scargle_period and autocorrelate_period and compares
        and combines the two outputs.
        '''
        periods = {}
        bad_kics = []
        for key in dict1.iterkeys():
            for n in [.2,.25,.33,.5,1,2,3,4,5]:
                if n*dict1[key][0] > dict2[key][0]*correlation[0] and n*dict1[key][0] < dict2[key][0]*correlation[1]:
                    periods[key] = [min(dict1[key][0],dict2[key][0]), dict1[key][1], dict1[key][2]]
                if n*dict2[key][0]/dict2[key][2] > .99 and n*dict2[key][0]/dict2[key][2] < 1.01:
                    bad_kics.append(key)
                    periods[key] = []
        for kic in bad_kics:
            del periods[kic]
        return periods

    
    def plot_all(self, xlim = (0,200), ylim = (0,10), dictionary = None, inputs = False , filename = None, minstrength = 0, maxstrength = 1, maskintegers = True):
        '''
        Takes either a dictionary (as returned above), or a filename (as written above) and
        plots orbital period/rotational period on the y-axis, and orbital period on the x-axis
        given the constraints set.
        '''

        if inputs:
            data = dictionary        
        else:
            data = self.pullfile(filename)

        kic = [x for x in data.iterkeys()]
        best_period = [x[0] for x in data.itervalues()]
        strength = [x[1] for x in data.itervalues()]
        p_orb = [x[2] for x in data.itervalues()]
    
        #p_orb_rot = [x[0]/x[1] for x in np.array([p_orb,best_period])]
        #still neeed to figure out if this can be done
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
        
        print len(p_orb_rot)
        
        plt.plot(p_orb, p_orb_rot, linestyle = 'none', marker='o', markersize=1)
        if inputs:
            source = 'period dictionary'
        else:
            source = filename
        plt.title('Plotted from ' + source + '\nshowing EBs with strength between ' + str(minstrength) + ' and ' + str(maxstrength)) 
        plt.ylabel('orbital period/rotation period')
        plt.xlabel('orbital period')
        plt.ylim(ylim)
        plt.xlim(xlim)
        plt.show()
        
    def file_merge(self, files, newfile):
        '''
        This takes multiple files with the format kic best_period strength orbital_period and merges them
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
    
    def make_period_files(self):
        '''
        Creates files containing all of the EBs for the autocorrelation and Lomb-Scargle functions 
        '''
        length = len(self.kics)/100
        for n in range(length):
            self.autocorrelate_period(self.kics[100*n:100+100*n], write = True, append = '_autocor_'+ str(n))
            self.lomb_scargle_period(self.kics[100*n:100+100*n], write = True, append = '_periodogram_' + str(n))
        self.autocorrelate_period(self.kics[length*100:], write = True, append = '_autocor_last')
        self.lomb_scargle_period(self.kics[length*100:], write = True, append = '_periodogram_last')

        files = []
        for x in range(0, length):
            files.append('periods_autocor_'+str(x))
        files.append('periods_autocor_last')
        self.file_merge(files, 'autocor_all')

        files = []
        for x in range(0,length):
            files.append('periods_periodogram_'+str(x))
        files.append('periods_periodogram_last')
        self.file_merge(files, 'periodogram_all')
        return None
