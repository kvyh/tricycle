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

    Attributes
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

    def cut_eclipses(self, window=1.0):
        """
        Return light curve with eclipses removed.

        Parameters
        ----------
        window : float, optional
            Width of the window around the eclipse to cut, in fractions of
            the eclipse duration. (Default: 1.0)

        Returns
        -------
        time_cut : ndarray
            Observation times.
        flux : ndarray
            Fluxes.
        err : ndarray
            Flux errors.

        """
        raise Exception('This module has not yet been written!')

    def periodogram(self, eclipses=False):
        """
        Return periodogram for light curve.

        Parameters
        ----------
        eclipses : bool, optional
            Include eclipses if True. (Default: False)

        Returns
        -------
        period : ndarray
            An array of periods.
        power : ndarray
            The power at the corresponding periods.

        """
        raise Exception('This module has not yet been written!')
