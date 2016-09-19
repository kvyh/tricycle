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
        """
        Parameters
        ----------
        kic_list : list of numbers
            list of kic numbers for the objects that will be analyzed
        """
        self.kics = kic_list
        self.autocor_results = {}
        self.periodogram_results = {}
        self.interpolate = True
        self.eclipsewidth = .7
    
    def plot_periodogram(self, kic, l_s=True, autocor=False):
        '''
        plots the periodogram of the object with kic number ``kic``
        '''
        time, flux, fluxerr, p_orb = self.get_info(kic)
        periods = periodicity2.Periodicity(time, flux, fluxerr)
        periods.Lomb_scargle()
        period = periods.period_power_ls[0]
        power = periods.period_power_ls[1]
        plt.plot(period, power, label='periodogram')
        plt.axvline(x=p_orb, label='orbital period', color='k')
        #period_list = self.potential_targets[kic]
        #plt.axhline(xmin=min(period_list), xmax=max(period_list), y=.5, label='target', color='g')
        plt.title(str(kic))
        plt.xlim(0, 45)
        plt.legend()
        plt.show()

    def get_info(self, kic):
        """
        Parameters
        ----------
        kic : int
            kic number of the target

        Returns
        -------
        time array, flux array, error array, orbital period for the target
        with the eclipses cut out (and interpolated over if self.interpolate)
        """
        curve = binaries.RealBinary(kic)
        time, flux, fluxerr, cadence, quarter, quality = data.loadlc_db(kic)
        if self.interpolate:
            time_cut, flux_cut, err_cut, quarter_cut = curve.interpolate(widthfactor=self.eclipsewidth)
        else:
            time_cut, flux_cut, err_cut, quarter_cut = curve.curve_cut(widthfactor=self.eclipsewidth)
        return time_cut, flux_cut, err_cut, curve.p_orb
    
    def rotation_periods(self, interpolate=True, eclipse_width=.7):
        """
        runs both the autocorrelation function and periodogram on each
        target in ``self.kics`` and records the results in ``self.autocor_results``
        and ``self.periodogram_results``
        Parameters
        ----------
        interpolate : bool
            if the eclipses should be cut, or cut and interpolated over
        eclipse_width : float in (0,1)
            how much of the eclipse should be cut (high values can cause errors
            in the functions if there are targets whose lightcurve consists almost
            entirely of eclipses (e.g. very short period binaries))

        Returns
        -------
        periodogram dictionary, autocorrelation dictionary
        of the form {kic#: [rotation period, strength, orbital period]}
        """
        self.interpolate = interpolate
        self.eclipsewidth = eclipse_width
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

    def l_s_multi_period(self, interpolate=True, eclipse_width=.7):
        """
        as above, but finds multiple periods in the periodogram

        Returns
        -------
        peirodogram dictionary of the form {kic#: [rot period, strength,
        rot period, strength, rot period, strength, rot period, strength,
        rot period, strength, orbital period
        """
        self.interpolate = interpolate
        self.eclipsewidth = eclipse_width
        self.periodogram_multi_results = {}
        for kic in self.kics:
            time, flux, fluxerr, p_orb = self.get_info(kic)
            periods = periodicity2.Periodicity(time, flux, fluxerr)
            periods.multi_period_l_s()
            periodogram_results = periods.periodogram_results
            results = []
            for result in periodogram_results:
                results.append(result[0])
                results.append(result[1])
            results.append(p_orb)
            self.periodogram_multi_results[kic] = results
        return self.periodogram_multi_results

    def read_multi_period(self, filename):
        """
        reads a file containing multiple strong periods per object
        """
        array = np.genfromtxt(filename)
        tempdict = {}
        for line in array:
            tempdict[line[0]] = [[line[1], line[2]], [line[3], line[4]],
                                 [line[5], line[6]], [line[7], line[8]],
                                 [line[9], line[10]], line[11]]
        self.periodogram_multi_results = tempdict
        return tempdict


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
        plots the lightcurve with and without eclipses on the same graph
        '''
        curve = binaries.RealBinary(kic)
        self.interpolate = True
        time, flux, fluxerr, cadence, quarter, quality = data.loadlc_db(kic)
        time_cut, flux_cut, err_cut, p_orb = self.get_info(kic)
        #plt.figure(1, figsize = (16,6))
        #plt.subplot(122)
        plt.plot(time_cut, flux_cut, label='eclipses cut')
        plt.legend()
        
        #plt.subplot(121)
        plt.plot(time, flux, label='with eclipses')
        plt.legend()
        plt.title(str(kic))
        #plt.show()
            
    def readfile(self, filename):
        """
        reads in a file of the same type that the ``self.write`` method creates
        and stores it in the appropriate ``self.autocor_results`` or
        ``self.periodogram_results``
        """
        kic, period, strength, p_orb = np.genfromtxt(filename, unpack=True)

        tempdict = {}
        for x in range(len(kic)):
            tempdict[int(kic[x])] = [period[x], strength[x], p_orb[x]]
        
        if 'autocor' in filename:
            self.autocor_results = tempdict
            print 'autocor read'
        elif 'periodo' in filename:
            self.periodogram_results = tempdict
            print 'periodogram read'
        elif 'better' in filename:
            self.better_period = tempdict
            print 'better read'
        return tempdict
    
    def better_periods(self, files=None):
        """
        uses the data contained in files or the output of self.rotation_periods
        and compares the results from the autocorrelation and the periodogram,
        returning only targets where they match
        """
        if not files == None:
            for filename in files:
                self.readfile(filename)
        if len(self.autocor_results)<1 or len(self.periodogram_results)<1:
            print 'period files are needed'
        self.better_period = {}
        for kic in self.autocor_results.iterkeys():
            auto = self.autocor_results[kic]
            p_gram = self.periodogram_results[kic]
            ratio = self.autocor_results[kic][0]/self.periodogram_results[kic][0]
            # check if they are the same, or either one registered a harmonic
            for scalar in [.2, .25, .33, .5, 1, 2, 3, 4, 5]:
                if .98 < ratio*scalar < 1.02:
                    self.better_period[kic] = [min(auto[0], p_gram[0]), auto[1], auto[2]]
                    #using the minimum will remove higher harmonics, but may show lower harmonics
        return self.better_period
                
    def test_data(self):
        '''
        uses better_periods to restrict targets to those where we are confident
        in the found period, then prints the deviation in the phases found by
        test_period in periodicity2
        '''
        one=[]
        point8=[]
        others=[]
        for kic,k in self.better_period.iteritems():
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
            phase_amp = periods.test_period(self.better_period[kic][0],'kic')
            print np.std([x[0] for x in phase_amp.itervalues()])
        
    def histogram(self, dictionary_array):
        """
        plots a histogram with bins based on the
        orbit period/rotation period ratio of targets

        Parameters
        ----------
        dictionary_array : array of dictionaries
            array of dictionaries with values of
            [rotation period, strength, orbit period]
        """
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
        """
        plots a histogram with bins based on the
        orbit period/rotation period ratio of targets

        Parameters
        ----------
        dictionary_array : array of dicts
            array of dictionaries with values of
            [rotation period, strength, orbit period]
        """
        #if len(self.autocor_results) < 1:
        #    self.find_periods()
        it = np.nditer(dictionary_array, flags=['refs_ok','external_loop'], order='F')
        plt.figure(figsize = (16,10))
        n = 0
        for dictionary, name in it:
            ratio = [x[2]/x[0] for x in dictionary.itervalues()]
            p_orb =[x[2] for x in dictionary.itervalues()]           
            plt.subplot(2,2,n+1)
            plt.title(str(name))
            # uses cmax to prevent the high-density bins around ratio==1 from
            # skewing the color scale (if they are left in, bins with 3 objects
            # have the same color as those with 10 objects)
            plt.hist2d(p_orb,ratio, range=[[0,15],[0,2.5]], bins = [80,100], cmax = 10)
            n+=1
        plt.legend()
        plt.show()            
    
    def plot_results(self, dictionary=None, autocorrelation=True, xlim=(0, 200),
                     ylim=(0, 10), minstrength=0, maxstrength=1):
        """

        Parameters
        ----------
        dictionary : dict
            dictionary of form {kic#: [rotation period, strength, orbital period]}
        autocorrelation : bool
            if True, uses ``self.autocorrelation_results``, if false
            uses ``self.periodogram_results``
        xlim : tuple
            x limits for the plot
        ylim : tuple
            y limits for the plot
        minstrength : float in (0,1)
            lower bound on allowed strength (strength is the value returned by
            the period-finding function that describes the likelyhood of the period)
            note: this is not a confidence percent
        maxstrength : float in (0,1)
            upper bound on allowed strength (strength is the value returned by
            the period-finding function that describes the likelyhood of the period)
            note: this is not a confidence percent
        Returns
        -------
        scatterplot of the (orbit period, ratio) of each object
        """
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
        
        print '# of points: ' + str(len(p_orb_rot))
        
        plt.plot(p_orb, p_orb_rot, linestyle = 'none', marker='o', markersize=1)
        plt.title('Eclipsing binaries with strength between ' + str(minstrength) + ' and ' + str(maxstrength)) 
        plt.ylabel('orbital period/rotation period')
        plt.xlabel('orbital period')
        plt.ylim(ylim)
        plt.xlim(xlim)
        #plt.show()
    
    def write(self, data, base, append):
        """
        creates a text file with the information stored in ``data`` that
        has a filename of ``base``_periods_all``append``

        Parameters
        ----------
        data : dict
            dictionary of style {kic#: [rotation period, strength, orbital period]}
        base : string
            start of the filename
        append : string
            end of the filename
        """


        kics = [x for x in data.iterkeys()]
        if len(data[kics[0]]) == 3:
            best_period = [x[0] for x in data.itervalues()]
            strength = [x[1] for x in data.itervalues()]
            p_orb = [x[2] for x in data.itervalues()]

            array = np.array([[kics[0], best_period[0], strength[0], p_orb[0]]])

            for n in range(1,len(kics)):
                array = np.append(array,[[kics[n], best_period[n], strength[n], p_orb[n]]],axis=0)
            np.savetxt(base+'_periods_all'+append, array, fmt ='%.18s')

        else:
            x = np.zeros((len(data[kics[0]])+1,len(kics)))
            x[0] = kics
            for n in range(len(data[kics[0]])):
                x[n+1] = [a[n] for a in data.itervalues()]
            array = x.T
            np.savetxt(base+'_periods_all'+append, array, fmt = '%.18s')

    def write_results(self, append):
        """
        writes self.autocor_results and self.periodogram_results
        to files that are named autocor_periods_all``append`` and
        periodogram_periods_all``append``

        Parameters
        ----------
        append : string
            ending for the filename
        """
        if len(self.autocor_results) > 0:
            pass
        else:
            self.find_periods()
        self.write(self.autocor_results, 'autocor', append)
        self.write(self.periodogram_results, 'periodogram', append)

    def results(self):    
        '''
        prints the contents of self.autocor_results and
        self.periodogram_results
        '''
        print 'autocorrelation results:', self.autocor_results
        print 'periodogram results:', self.periodogram_results
        return self.autocor_results, self.periodogram_results

        
def main():
    # this was written to create files with the results of the
    # periodogram and autocorrelation with and without interpolation
    # over the removed eclipses. It failed while running 0.75, likely
    # due to having too much of the lightcurve removed when the eclipses
    # were removed
    allkics = data.select_kics()
    ka = kic_analyze(allkics[:])
    for x in [0.3, 0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8]:
        ka.rotation_periods(interpolate = True, eclipse_width= x)
        ka.write_results('_interp_2_'+str(x))
        ka.rotation_periods(interpolate = False, eclipse_width= x)
        ka.write_results('_cut_2_'+str(x))
    
    #ka.readfile()
    #ka.find_periods()
    #ka.results()
    #ka.lightcurve(allkics[3])
    #ka.plot_results(periodograms = False)
    
if __name__ == '__main__':
    main()  
