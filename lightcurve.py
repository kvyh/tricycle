from gatspy.periodic import LombScargleFast
import matplotlib.pyplot as plt

class LightCurve(object):
    """
    Eclipsing binary light curve with system parameters.

    Parent class with analysis methods.

    Parameters
    ----------
    time : ndarray
        The observation times, in days.
    flux : ndarray
        The system flux. May be in physical or relative units.
    err : ndarray
        The flux errors. Should have same units as `flux`.
    p_orb : float
        Orbital period in days.
    t_0 : float
        Reference time of mid-eclipse in days.
    p_depth : float
        Depth of primary eclipse.
    s_depth : float
        Depth of secondary eclipse.
    p_width : float
        Width of primary eclipse in phase.
    s_width : float
        Width of secondary eclipse in phase.
    sep : float
        Seperation between the primary eclipse in secondary eclipse in phase
    kic : int
        Kepler Input Catalog (KIC) ID number.

    """
    def __init__(self, time, flux, err, p_orb, t_0, p_depth, s_depth, p_width,
                 s_width, sep, kic):

        self.time = time
        self.flux = flux
        self.err = err
        self.p_orb = p_orb
        self.t_0 = t_0
        self.p_depth = p_depth
        self.s_depth = s_depth
        self.p_width = p_width
        self.s_width = s_width
        self.sep = sep
        self.kic = kic

    def get_period(self, units='days'):
        """
        Return the orbital period.

        Parameters
        ----------
        units : {'days', 'hours', 'minutes'}
            Desired units for period.

        Returns
        -------
        period : float
            The orbital period in the specified units.

        """
        if units == 'days':
            period = self.p_orb
        elif units == 'hours':
            period = self.p_orb * 24.
        elif units == 'minutes':
            period = self.p_orb * 24. * 60.
        else:
            raise ValueError('Invalid choice of units.')

        return period

    def curve_cut(self):
        """This a function that takes parameters of a star and gives you the time and flux with the eclipses cut out
    
        Parameters:
        Time
        Flux
        bjdo
        period
        pwidth
        swidth
        separation
        
        return:
            timecut and fluxcut"""
        # Convert BJD to Kepler date
        t_0 = self.t_0 - 54833

        phase=((self.time- t_0) % self.p_orb) / self.p_orb
        mask= ((phase > self.p_width / 2.) &
               (phase < 1 - self.p_width / 2.)) & \
              ((phase > self.sep + self.s_width / 2.) |
               (phase < self.sep - self.s_width / 2.))
               
        timecut= self.time[mask]
        fluxcut= self.flux[mask]
        errcut = self.err[mask]

        return timecut, fluxcut, errcut

    def periodogram(self, time, flux, err, p_fold=None, plt_color='k',
                    max_days=100.0, oversampling=5, plot=False, cut_eclipses=True, best_period = True, period_range = (.05,45)):
        """
        Plot a the light curve, periodogram, and phase-folded light curve.

        The orbital period and possible aliases are indicated in red on the
        periodogram.

        Parameters
        ----------
        time : array_like
            Observation times in days.
        flux : array_like
            Fluxes.
        err : array_like
            Flux errors.
        p_fold : float, optional
            Plot the light curve folded at this period (in days), and indicate
            the period and likely aliases on the periodogram plot.
        plt_color : str, optional
            The line and scatter plot color.
        max_days : float, optional
            The maximum number of days to plot in the light curve and
            periodogram.
        oversampling: int, optional
            The oversampling factor for the periodogram.

        """
        if cut_eclipses:
            time, flux, err = self.curve_cut()
            
        else:
            time=self.time
            flux=self.flux
            error=self.error
            
        model = LombScargleFast().fit(time, flux, err)
        period, power = model.periodogram_auto(oversampling=oversampling)
        
         #new stuff to output the best period
        if best_period:
            model.optimizer.period_range = period_range
            Best_period = model.best_period
            
        if plot:
            fig1, (ax1, ax2) = plt.subplots(nrows=2, figsize=(7, 14))
    
            ax1.plot(time, flux, color=plt_color)
            ax1.set_xlim(time.min(), time.min() + max_days)
            ax1.set_xlabel('Time (days)')
            ax1.set_ylabel('Relative Flux')
    
            ax2.plot(period, power, color=plt_color)
            ax2.set_xlim(0.1, max_days)
            ax2.set_xlabel('Period (days)')
            ax2.set_ylabel('Power')
    
            # Plot some the most likely aliases of eclipse period.
            factors = [0.5, 1, 2, 3, 4]
            for factor in factors:
                ax2.axvline(self.p_orb / factor, color='r')
    
            if p_fold is not None:
                for factor in factors:
                    ax2.axvline(p_fold / factor, color='b')
    
            fig1.suptitle('KIC {0:d} --- p_orb = {1:3.5f} days'.
                          format(self.kic, self.p_orb))
    
            if p_fold is not None:
                fig2, ax3 = plt.subplots()
                phase = (time % p_fold) / p_fold
                ax3.scatter(phase, flux, color=plt_color, s=0.1)
                ax3.set_xlim(0, 1)
                ax3.set_xlabel('Phase')
                ax3.set_ylabel('Relative Flux')
    
            plt.show()
        #I don't know how it would react to returning nonexistent stuff
        if best_period:
            return period, power, Best_period
        else:
            return period, power    
