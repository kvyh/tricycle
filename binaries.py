import numpy as np

from lightcurve import LightCurve


class ModelBinary(LightCurve):
    """
    Model eclipsing binary light curve.

    See parent LightCurve class documentation of shared parameters.

    Parameters
    ----------
    p_rot1 : float
        Rotation period of star 1 in days.
    p_rot2 : float
        Rotation period of star 2 in days.
    amp_1 : float
        Relative spot amplitude of star 1.
    amp 2 : float
        Relative spot amplitude of star 2.
    length : float
        Length of desired light curve in days.
    sig: float
        Standard deviation of errors, in relative flux units.
    sc : bool
        Default is for short cadence. Set to False for long cadence.

    """
    def __init__(self, p_orb=7.07, t_0=0.0, e_dur=0.27, depth=0.5, p_rot1=1.51,
                 p_rot2=0.80, amp_1=0.03, amp_2=0.02, length=100.0, sig=0.02,
                 sc=True):

        time, flux, err = self.__make_lc(p_orb, t_0, e_dur, depth, p_rot1,
                                         p_rot2, amp_1, amp_2, length, sig, sc)

        super(ModelBinary, self).__init__(time, flux, err, p_orb, t_0, e_dur,
                                          depth)

    @staticmethod
    def __make_lc(p_orb, t_0, e_dur, depth, p_rot1, p_rot2, amp_1, amp_2,
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
