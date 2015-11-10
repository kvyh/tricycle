#!/usr/bin/env python
# -*- coding: utf-8 -*-
from collections import OrderedDict
import numpy as np
import jktebop_f2py

class BopIt(object):
    """
    Simulate an eclipsing binary using JKTEBOP.

    Parameters
    ----------
    sine_star : tuple, optional
        The star that each sine belongs to. -1 for primary, -2 for secondary.
    sine_params : tuple, optional
        For each sine, provide the time of zero phase (HJD), period (days),
        and amplitude (relative flux) as (hjd_0, period, amp).
    ldtype : tuple, optional
        Type of limb darkening.

    Other Parameters
    ----------------
    csbr : float, optional
        Central surface brightness ratio. (Default: 1.0)
    sfsr : float, optional
        Sum of fractional stellar radii, i.e. scaled by semimajor axis.
        (Default: 0.1)
    rsr : float, optional
        Ratio of stellar radii. (Default: 1.0)
    ld_a :
    ld_b :
    inc :
    ecc :
    omg_bar :
    gd_a :
    gd_b :
    rl_a :
    rl_b :
    mr :
    tll :
    l3 :
    pcf :
    lsf :
    irs :
    p_orb :
    t_0 :
    ld2_a :
    ld3_a :
    ld4_a :
    ld2_b :
    ld3_b :
    ld4_b :
    rv_a :
    rv_b :
    sys_a :
    sys_b :

    Examples
    --------
    >>> bop_it = BopIt()
    >>> time, flux = bop_it.light_it()

    """
    def __init__(self, sine_star=None, sine_params=None, ldtype=(1, 1), **kwargs):

        self.ldtype = ldtype

        # An ordered dictionary of optional keywords for JKTEBOP.
        self.params = OrderedDict([('csbr', 1.0),
                                  ('sfsr', 0.1),
                                  ('rsr', 1.0),
                                  ('ld_a', 0.0),
                                  ('ld_b', 0.0),
                                  ('inc', 90.0),
                                  ('ecc', 0.0),
                                  ('omg_bar', 0.0),
                                  ('gd_a', 1.0),
                                  ('gd_b', 1.0),
                                  ('rl_a', 0.0),
                                  ('rl_b', 0.0),
                                  ('mr', 1.0),
                                  ('tll', 0.0),
                                  ('l3', 0.0),
                                  ('pcf', 0.0),
                                  ('lsf', 0.0),
                                  ('irs', 1.0),
                                  ('p_orb', 3.35),
                                  ('t_0', 54833.0),
                                  ('ld2_a', 0.0),
                                  ('ld3_a', 0.0),
                                  ('ld4_a', 0.0),
                                  ('ld2_b', 0.0),
                                  ('ld3_b', 0.0),
                                  ('ld4_b', 0.0),
                                  ('rv_a', 0.0),
                                  ('rv_b', 0.0),
                                  ('sys_a', 0.0),
                                  ('sys_b', 0.0)])

        # Parse any input options.
        for key, value in self.params.iteritems():
            self.params[key] = kwargs.pop(key, value)

        # Make sure that we didn't get any unknown options.
        if len(kwargs):
            raise TypeError("__init__() got an unexpected keyword argument "
                            "'{0}'".format(kwargs.keys()[0]))

        # Make parameter array for JKTEBOP.
        self.vv = np.zeros(138, dtype=float)
        for ii, param in enumerate(self.params.values()):
            self.vv[ii] = param

        # Deal with sine parameters.
        self.psine = np.zeros(9)
        self.nsine = 0

        if sine_star is not None and sine_params is not None:

            if len(sine_star) is not len(sine_params):
                raise ValueError('Number of sines does not equal number of '
                                 'parameter sets.')

            # Add sine parameters to parameter array.
            for idx, sine_params in enumerate(sine_params):
                self.vv[30 + idx * 3:33 + idx * 3] = sine_params

            self.psine[0:len(sine_star)] = sine_star
            self.nsine = len(sine_star)

        elif sine_star is None and sine_params is not None:
            raise ValueError("Use 'sine_star' to specify which star(s) to "
                             "apply sine(s) to")
        elif sine_star is not None and sine_params is None:
            raise ValueError("'sine_params' not given")

    def light_it(self, length=90.0, sigma=None, sc=False):
        """
        Produce a Kepler-like light curve.

        Parameters
        ----------
        length : float, optional
            Length of the desired light curve in days. (Default: 90.0)
        sigma : float, optional
            Standard deviation of the noise, in relative flux units.
        sc : bool, optional
            Default is long cadence (30 min. data). Set to True for short
            cadence (1 min data).

        Returns
        -------
        time : array_like
            Observation times in HJD.
        flux : array_like
            Relative flux of EB.

        """
        # Additional filler for JKTEBOP.
        vary = np.zeros(138, dtype=float)
        npoly = 0
        ppoly = np.zeros(9)
        dtype = 1
        l_a = 0.0
        l_b = 0.0
        num_int = 1
        n_interval = 0

        # Set cadence.
        if sc:
            dt = 1.0 / 24.0 / 60.0
        else:
            dt = 1.0 / 24.0 / 2.

        # Create time and flux arrays.
        time = np.arange(self.params['t_0'], self.params['t_0'] + length + dt,
                         dt)
        flux = np.empty_like(time)

        # Loop through time and run JKTEBOP.
        for ii, tt in enumerate(time):
            mag = jktebop_f2py.getmodel(self.vv, vary, self.ldtype, self.nsine,
                                        self.psine, npoly, ppoly, tt, dtype,
                                        l_a, l_b, num_int, n_interval)
            flux[ii] = 10.0 ** (mag / - 2.5)

        # Add noise.
        if sigma is not None:
            flux += np.random.normal(0, sigma, len(flux))

        return time, flux