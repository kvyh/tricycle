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

    """
    def __init__(self, time, flux, err, p_orb, t_0, p_depth, s_depth, p_width,
                 s_width):

        self.time = time
        self.flux = flux
        self.err = err
        self.p_orb = p_orb
        self.t_0 = t_0
        self.p_depth = p_depth
        self.s_depth = s_depth
        self.p_width = p_width
        self.s_width = s_width

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
