import numpy as np
import matplotlib.pyplot as plt
import data
import file_merge as fm
from gatspy.periodic import LombScargleFast
from interpacf import interpolated_acf, dominant_period
from scipy.ndimage import gaussian_filter
from scipy import signal
from scipy import optimize

#This code uses modules from https://github.com/StellarArmy/tricycle and
#https://github.com/bmorris3/interp-acf and requires the gatspy and MySQLdb
#packages to run properly


class Periodicity():
    '''
    Takes a lightcurve as input and contains methods
    for analysis of the lightcurve
    '''
    def __init__(self, time, flux, fluxerr):
        """

        Parameters
        ----------
        time : array of times
            The times at which the data were taken
        flux : array of fluxes
            The data associated with the times
        fluxerr : array of flux errors
            The errors in the data
        """
        self.time = time
        self.flux = flux
        self.err = fluxerr
        self.autocor_results = []
        self.periodogram_results = []

    def autocorrelate(self, min_days=0.2, max_days=45.0):
        '''
        This takes the ``lightcurve`` and runs it through interp-acf's
        autocorrelation function to find the most probable period
        within the bounds of ``min_days`` and ``max_days``
        '''
        flux = self.flux
        time = self.time
        err = self.err        
        # as described in the interpacf code:
        flux -= np.mean(self.flux)
        lag, acf = interpolated_acf(time, flux)
        detected_period = -0.1
        strength = -1.0
        try:
            detected_period, strength = dominant_period(lag, acf, min_days, max_days)
        except TypeError:
            print 'TypeError during autocorrelation'
        self.autocor_results = [detected_period, strength]
            
    def Lomb_scargle(self, oversampling=5, period_range=(.05, 45)):
        '''
        Uses gatspy's ``LombScargleFast`` to create a periodogram, stores
        the period, power values in ``self.period_power_ls`` and the strongest
        period, and its strength, in ``self.periodogram_results``.
        '''
        # Copied from lightcurve.py from Tricycle and the gatspy documentation
        
        time=self.time
        flux=self.flux
        err=self.err
        
        model = LombScargleFast().fit(time, flux, err)
        period, power = model.periodogram_auto(oversampling=oversampling)
        self.period_power_ls = [period, power]
        model.optimizer.period_range = period_range
        detected_period = model.best_period
        where = np.abs(period - detected_period).argmin()
        strength = power[where]
        self.periodogram_results = [detected_period, strength]

    def multi_period_l_s(self, oversampling=5, period_range=(.05, 45)):
        time = self.time
        flux = self.flux
        err = self.err
        model = LombScargleFast().fit(time, flux, err)
        period, power = model.periodogram_auto(oversampling=oversampling)
        self.period_power_ls = [period, power]
        model.optimizer.period_range = period_range
        model.optimizer.quiet = True
        periods, scores = model.optimizer.find_best_periods(model, return_scores=True)
        self.periodogram_results = [[periods[i], scores[i]] for i in range(len(periods))]
        
    def targeted_LS(self, oversampling=5, target_range=None, kic='no kic given'):
        '''
        target_range should be a tuple. For tricycle equal to (p_orb/.85,p_orb/.7)
        '''
        from scipy.signal import argrelextrema
        time = self.time
        flux = self.flux
        err = self.err
        
        model = LombScargleFast().fit(time, flux, err)
        period, power = model.periodogram_auto(oversampling=oversampling)
        powmax = max(power)
        maximask = argrelextrema(power, np.greater)
        minmask = np.where(power < 9.*powmax/10)
        validmaxes = (np.array(np.intersect1d(minmask, maximask)),)
        self.targeted_ls = np.array([period[validmaxes],power[validmaxes]])
        interested = [a[0] for a in np.transpose(self.targeted_ls) if a[0]<target_range[1] and a[0]>target_range[0] and a[1]>.05]
        if len(interested)>0 and min(interested)>.1 and min(interested)<45:
            #print kic + 'may have a period at: ' +str(interested)
            return interested
        
    def test_period(self, test_period, kic):
        """
        takes a period and a kic number and tries to fit the target's
        lightcurve to a cosine function of that period. It performs this
        fit on each ``test_period``-length interval and uses the results
        of the previous fit as defaults.

        returns
        -------
        an array with [phase offset, amplitude] of the fitted
        function for each interval.
        """
        int_time = self.time//test_period
        phase_amp = {}
        par = [.005, 0]
        x = 0
        for time in range(int(min(int_time))+1, int(max(int_time))):
            if len(np.where(int_time == float(time))[0]) > 1:
                mask = np.where(int_time == float(time))
                #if len(mask[0])>100:
                #    print len(mask[0])
                func = lambda par, t: par[0]*np.cos(2*np.pi/test_period*t+par[1]%(2*np.pi))
                funcerr = lambda par, t, x: func(par,t)-x
                P0 = par
                P0[1]=P0[1]%(2*np.pi)
                par, success = optimize.leastsq(funcerr,P0,args=(self.time[mask],self.flux[mask]-np.median(self.flux[mask])))
                if par[0] >= 1:
                    phase_amp[time] = [par[1] % (2*np.pi), par[0]]
                else:
                    phase_amp[time] = [(par[1]+np.pi) % (2*np.pi), -par[0]]
                if x == 1 or 20 or 50:
                    time = np.linspace(self.time[mask].min(), self.time[mask].max(), 100)
                    plt.plot(self.time[mask], self.flux[mask] - 1., linestyle='None', marker='o')
                    plt.plot(time, func(par, time))
                    plt.title(kic)
                    plt.show()
                x+=1
            #else:
            #    print 'optimize being stupid between',time*testperiod,',',(time+1)*testperiod
        
        return phase_amp

    def results(self):
        self.autocorrelate()
        self.Lomb_scargle()
        return self.periodogram_results, self.autocor_results

        
        
 