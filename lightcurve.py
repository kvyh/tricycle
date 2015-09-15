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
    e_dur : float
        Eclipse duration in days.
    depth : float
        Relative eclipse depth.

    """
    def __init__(self, time, flux, err, p_orb, t_0, e_dur, depth):

        self.time = time
        self.flux = flux
        self.err = err
        self.p_orb = p_orb
        self.t_0 = t_0
        self.e_dur = e_dur
        self.depth = depth

    def get_period(self, units='days'):
        """
        Return the orbital period.

        Parameters
        ----------
        units : str
            Desired units for period: days, hours, or minutes.

        Returns
        -------
        period : float
            The orbital period in the specified units.

        """
        if units == 'days':
            per = self.p_orb
        elif units == 'hours':
            per = self.p_orb*24.
        elif units == 'minutes':
            per = self.p_orb*24.*60.
        else:
            raise ValueError('Invalid choice of units.')

        return per
