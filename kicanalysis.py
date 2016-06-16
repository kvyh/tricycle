import binaries
import numpy as np
import matplotlib.pyplot as plt
import data
import periodicity2
from scipy import signal

#This code uses modules from https://github.com/StellarArmy/tricycle and
#https://github.com/bmorris3/interp-acf and requires the gatspy and MySQLdb
#packages to run properly

class kic_analyze():
    def __init__(self, kic_list):
        '''
        
        '''
        self.kics = kic_list
        self.autocor_results = {}
        self.periodogram_results = {}
        self.interpolate = True
        self.eclipsewidth = .7
    
    def plot_periodogram(self, kic, l_s = True, autocor = False):
        '''
        
        '''
        time, flux, fluxerr, p_orb = self.get_info(kic)
        periods = periodicity2.Periodicity(time, flux, fluxerr)
        periods.Lomb_scargle()
        period = periods.period_power_ls[0]
        power = periods.period_power_ls[1]
        plt.plot(period,power,label='periodogram')
        plt.axvline(x=p_orb,label='orbital period',color = 'k')
        period_list = self.potential_targets[kic]
        plt.axhline(xmin = min(period_list),xmax = max(period_list),y=.5,label='target',color = 'g')
        plt.xlim(0,45)
        plt.legend()
        plt.show()
        pass
        
    
    def get_info(self, kic):
        '''
        
        '''
        curve = binaries.RealBinary(kic)
        time, flux, fluxerr, cadence, quarter, quality = data.loadlc_db(kic)
        if self.interpolate:
            time_cut, flux_cut, err_cut, quarter_cut = curve.interpolate(widthfactor = self.eclipsewidth)
        else:
            time_cut, flux_cut, err_cut, quarter_cut = curve.curve_cut(widthfactor = self.eclipsewidth)
        return time_cut, flux_cut, err_cut, curve.p_orb
    
    def rotation_periods(self, interpolate = True, eclipsewidth = .7):
        '''
        
        '''
        self.interpolate = interpolate
        self.eclipsewidth = eclipsewidth
        self.periodogram_results = {}
        self.autocor_results = {}
        for kic in self.kics:
            time, flux, fluxerr, p_orb = self.get_info(kic)
            periods = periodicity2.Periodicity(time, flux, fluxerr)
            periodogram_result, autocor_result = periods.results()
            periodogram_result.append(p_orb)
            autocor_result.append(p_orb)
            self.periodogram_results[kic] = periodogram_result
            self.autocor_results[kic] = autocor_result
        return self.periodogram_results, self.autocor_results
    
    def list_potentials(self):
        #work in progress
        self.interpolate = True
        self.eclipsewidth = .7
        self.potential_targets = {}
        for kic in self.kics:
            time, flux, fluxerr, p_orb = self.get_info(kic)
            periods = periodicity2.Periodicity(time, flux, fluxerr)
            interested = periods.targeted_LS(target_range = (p_orb/.85,p_orb/.7), kic = kic)
            if interested:
                self.potential_targets[kic] = interested
        return self.potential_targets.keys()
    
    def lightcurve(self, kic):
        '''
        
        '''
        curve = binaries.RealBinary(kic)
        self.interpolate = True
        time, flux, fluxerr, cadence, quarter, quality = data.loadlc_db(kic)
        time_cut, flux_cut, err_cut, p_orb = self.get_info(kic)
        #plt.figure(1, figsize = (16,6))
        #plt.subplot(122)
        plt.plot(time_cut,flux_cut,label = 'eclipses cut')
        plt.legend()
        
        #plt.subplot(121)
        plt.plot(time,flux, label = 'with eclipses')
        plt.legend()
        plt.show()
            
    def readfile(self, filename):
        kic, period, strength, p_orb = np.genfromtxt(filename, unpack = True)

        tempdict = {}
        for x in range(len(kic)):
            tempdict[int(kic[x])] = [period[x], strength[x], p_orb[x]]
        
        if filename[0] == 'a':
            self.autocor_results = tempdict
        else:
            self.periodogram_results = tempdict
        return tempdict
    
    def better_periods(self, files = None):
        '''
        uses the data contained in files or output by self.rotation_periods
        '''
        if not files == None:
            for filename in files:
                self.readfile(filename)
        if len(self.autocor_results)<1 or len(self.periodogram_results)<1:
            print 'period files are needed'
        self.better_periods = {}
        for kic in self.autocor_results.iterkeys():
            auto = self.autocor_results[kic]
            perio = self.periodogram_results[kic]
            ratio = self.autocor_results[kic][0]/self.periodogram_results[kic][0]
            for scalar in [.2,.25,.33,.5,1,2,3,4,5]:
                if .98<auto[0]*scalar/perio[0] and auto[0]*scalar/perio[0]<1.02:
                    self.better_periods[kic] = [min(auto[0],perio[0]),auto[1],auto[2]]
                    #using the minimum will remove higher harmonics, but may show lower harmonics
        return self.better_periods
                
    def test_data(self):
        '''
        uses better_periods
        '''
        one=[]
        point8=[]
        others=[]
        for kic,k in self.better_periods.iteritems():
            if .99 < k[2]/k[0] and 1.01>k[2]/k[0]:
                one.append(kic)
            elif .7 < k[2]/k[0] and .85 > k[2]/k[0]:
                point8.append(kic)
            else:
                others.append(kic)
        testkics = one[:15]+point8[:15]+others[:15]
        for kic in testkics:
            time, flux, fluxerr, p_orb = self.get_info(kic)
            periods = periodicity2.Periodicity(time, flux, fluxerr)
            phase_amp = periods.test_period(self.better_periods[kic][0],'kic')
            print np.std([x[0] for x in phase_amp.itervalues()])
        
    def histogram(self, dictionary_array):
        '''
        
        '''
        #if len(self.autocor_results) < 1:
        #    self.find_periods()
        it = np.nditer(dictionary_array, flags=['refs_ok','external_loop'], order='F')
        plt.figure(1,figsize = (16,10))
        n = 0
        colors = ['r','c','g','m','b','k']
        for dictionary, name in it:
            ratio = [x[2]/x[0] for x in dictionary.itervalues()]
            plt.subplot()
            plt.hist(ratio, range=(0,5), bins = 200,label = str(name),histtype='bar', ec = colors[n],fill=False, alpha = .5)
            n+=1
        plt.legend()
        plt.show()
        
    def histogram2D(self, dictionary_array):
        '''
        needs work on array size
        '''
        #if len(self.autocor_results) < 1:
        #    self.find_periods()
        it = np.nditer(dictionary_array, flags=['refs_ok','external_loop'], order='F')
        plt.figure(figsize = (16,10))
        n = 0
        colors = ['r','c','g','m','b','k']
        for dictionary, name in it:
            ratio = [x[2]/x[0] for x in dictionary.itervalues()]
            p_orb =[x[2] for x in dictionary.itervalues()]           
            plt.subplot(2,2,n+1)
            plt.title(str(name))
            plt.hist2d(p_orb,ratio, range=[[0,15],[0,2.5]], bins = [80,100], cmax = 10)
            n+=1
        plt.legend()
        plt.show()            
    
    def plot_results(self, dictionary = None, autocorrelation = True, xlim = (0,200), ylim = (0,10) ,minstrength = 0, maxstrength = 1):
        '''
        
        '''
        if len(self.autocor_results) == 1:
            self.find_periods()
        if dictionary == None:
            if autocorrelation:
                data = self.autocorrelation_results
            else:
                data = self.periodogram_results
        else:
            data = dictionary
        
        kic = np.array([x for x in data.iterkeys()])
        best_period = np.array([x[0] for x in data.itervalues()])
        strength = np.array([x[1] for x in data.itervalues()])
        p_orb = np.array([x[2] for x in data.itervalues()])
    
        p_orb_rot = []
        for i in range(len(kic)):
            p_orb_rot.append(p_orb[i]/best_period[i])

        mask = np.where((strength >= minstrength) & (strength <= maxstrength))
        #mask2 = np.where(strength <= maxstrength)
        #mask = mask & mask2
        
        p_orb_rot = np.array(p_orb_rot)[mask]
        p_orb = p_orb[mask]
        strength = strength[mask]
        
        print len(p_orb_rot)
        
        plt.plot(p_orb, p_orb_rot, linestyle = 'none', marker='o', markersize=1)
        plt.title('Eclipsing binaries with strength between ' + str(minstrength) + ' and ' + str(maxstrength)) 
        plt.ylabel('orbital period/rotation period')
        plt.xlabel('orbital period')
        plt.ylim(ylim)
        plt.xlim(xlim)
        plt.show()
    
    def write(self, data, base, append):
        '''
        
        '''
        kics = [x for x in data.iterkeys()]
        best_period = [x[0] for x in data.itervalues()]
        strength = [x[1] for x in data.itervalues()]
        p_orb = [x[2] for x in data.itervalues()]
        
        array = np.array([[kics[0], best_period[0], strength[0], p_orb[0]]])
               
        for n in range(1,len(kics)):
            array = np.append(array,[[kics[n], best_period[n], strength[n], p_orb[n]]],axis=0)
        np.savetxt(base+'_periods_all'+append, array, fmt ='%.18s')
    
    def write_results(self, append):
        '''
        
        '''
        if len(self.autocor_results) > 0:
            pass
        else:
            self.find_periods()
        self.write(self.autocor_results, 'autocor', append)
        self.write(self.periodogram_results, 'periodogram', append)
        

    def results(self):    
        '''
        
        '''
        print 'autocorrelation results:', self.autocor_results
        print 'periodogram results:', self.periodogram_results

        
def main():
    allkics = data.select_kics()
    ka = kic_analyze(allkics[:])
    for x in [.75,.8,.9,1]:
        ka.rotation_periods(interpolate = True, eclipsewidth = x)
        ka.write_results('_interp_'+str(x))
        ka.rotation_periods(interpolate = False, eclipsewidth = x)
        ka.write_results('_cut_'+str(x))
    
    #ka.readfile()
    #ka.find_periods()
    #ka.results()
    #ka.lightcurve(allkics[3])
    #ka.plot_results(periodograms = False)
    
if __name__ == '__main__':
    main()  
