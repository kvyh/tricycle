import numpy as np

from lightcurve import LightCurve
import data
import binaries
import matplotlib.pyplot as plt
import pandas as pd

class ModelBinary(LightCurve):
    """
    Model eclipsing binary light curve.

    See parent LightCurve class documentation of shared parameters.

    Parameters
    ----------
    p_rot : float
        Rotation period of star occulted during primary eclipse.
    s_rot : float
        Rotation period of star occulted during secondary eclipse.
    p_amp : float
        Relative spot amplitude of star occulted during primary eclipse.
    s_amp : float
        Relative spot amplitude of star occulted during secondary eclipse.
    length : float
        Length of desired light curve in days.
    sig: float
        Standard deviation of errors, in relative flux units.
    sc : bool
        Default is for short cadence. Set to False for long cadence.

    """
    def __init__(self, p_orb=7.07, t_0=0.0, p_depth=0.23, s_depth=0.17,
                 p_width=0.02, s_width=0.02, p_rot=1.51, s_rot=0.80,
                 p_amp=0.03, s_amp=0.02, length=100.0, sig=0.02, sc=True):

        time, flux, err = self._make_lc(p_orb, t_0, p_width * p_orb, p_depth,
                                        p_rot, s_rot, p_amp, s_amp, length,
                                        sig, sc)

        super(ModelBinary, self).__init__(time, flux, err, p_orb, t_0, p_depth,
                                          s_depth, p_width, s_width)

    # TODO: Light curve model should account for both primary and secondary
    # eclipse properties.
    @staticmethod
    def _make_lc(p_orb, t_0, e_dur, depth, p_rot1, p_rot2, amp_1, amp_2,
                 length, sig, sc):
        """
        Make a light curve.

        """
        # One minute exposure time for short cadence, 30 minute for long.
        if sc:
            exptime = 1.0 / 60.0 / 24.0
        else:
            exptime = 30.0 / 60.0 / 24.0
        time = np.arange(0, length, exptime)

        # Make the starspot signals.
        f_1 = -amp_1 * np.sin(time * (2 * np.pi) / p_rot1)
        f_2 = -amp_2 * np.sin(time * (2 * np.pi) / p_rot2)

        # Make the eclipse signal as a step function.
        in_e = np.where((((time - t_0) % p_orb) / p_orb < (e_dur / p_orb)))
        f_e = np.zeros_like(time)
        f_e[in_e] = -depth

        # Make random noise.
        noise = np.random.random(len(time)) * sig

        # Add signals together.
        flux = 1 + f_e + f_1 + f_2 + noise

        return time, flux, noise

class RealBinary(LightCurve):
    """ 
    It's an eclipsing binary from Kepler's data
    """
    def __init__(self,kic):
        #This calls the data file
        catfile='villanova-db.csv'
        df = pd.read_csv(catfile)
    
        #Calls the row and stores it into a new data file that only has the row's information
        dfnew = df[df['KIC']==kic]
        #We need dfnew to be greater than zero
        #Calls each individual componet from the row
        if len(dfnew) ==0:
            print "no"
        else:
            print 'yasss'
        
        kicval=dfnew.KIC.values[0]
        KOI=dfnew.KOI.values[0]
        mult=dfnew.Mult.values[0]
        period=dfnew.period.values[0]
        swidth=dfnew.swidth.values[0]
        bjd0= dfnew.bjd0.values[0]
        pwidth=dfnew.pwidth.values[0]
        pdepth=dfnew.pdepth.values[0]
        sdepth=dfnew.sdepth.values[0]
        sep=dfnew.sep.values[0]
        time, flux, fluxerr, cadence, quarter, quality = data.loadlc_db(kic)
        super(RealBinary, self).__init__(time, flux, fluxerr, period, bjd0, pdepth,
                                    sdepth, pwidth, swidth, sep)
        