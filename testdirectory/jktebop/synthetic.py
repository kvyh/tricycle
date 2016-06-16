#!/usr/bin/env python
# -*- coding: utf-8 -*-
from collections import OrderedDict
import random
import signal

import numpy as np

import jktebop_f2py


class BopIt(object):
    """
    Simulate an eclipsing binary using JKTEBOP.

    JKTEBOP written by John Southworth
    http://www.astro.keele.ac.uk/jkt/codes/jktebop.html

    A description of the parameters are below. Consult the orignal FORTRAN code
    for more extensive documentation. Code copied from Ian Crossfield for
    setting up FORTRAN interface and to Dan Forman-Mackey for handling
    keyword arguments.

    Parameters
    ----------
    sine_star : tuple, optional
        The star that each sine belongs to. -1 for star A, -2 for star B.
    sine_params : tuple, optional
        For each sine, provide the time of zero phase (HJD), period (days),
        and amplitude (relative flux) as (hjd_0, period, amplitude).
    ldtype : 2-sequence of ints, optional
        Type of limb darkening for stars (A, B). (Default: (4, 4))
        1: linear
        2: linear + logarithmic
        3: linear + square root
        4: linear + quadratic
        5: linear + cubic
        6: Claret's four parameter form

    Other Parameters
    ----------------
    csbr : float, optional
        Central surface brightness ratio. (Default: 1.0)
    sfsr : float, optional
        Sum of fractional stellar radii, i.e. scaled by semimajor axis.
        (Default: 0.1)
    rsr : float, optional
        Ratio of stellar radii. (Default: 1.0)
    ld_a : float, optional
        Linear limb darkening coefficient for star A. Default of 0.3985 is from
        Sing et al. (2010) for quadratic law, and a 5750K dwarf.
    ld_b : float, optional
        Linear limb darkening coefficient for star B. (Default: 0.3985)
    inc : float, optional
        Orbital inclination in degrees. (Default: 90.0)
    ecc : float, optional
        Orbital eccentricity. (Default: 0.0)
    omg_bar : float, optional
        Longitude of periapsis in degrees. (Default: 0.0)
    gd_a : float, optional
        Gravity darkening of star A. (Default: 0.0)
    gd_b : float, optional
        Gravity darkening of star B. (Default: 0.0)
    rl_a : float, optional
        Reflected light for star A. (Default: 0.0)
    rl_b : float, optional
        Reflected light for star B. (Default: 0.0)
    mr : float, optional
        Mass ratio (M_b / M_a). Only used to calculate deformations. Set to
        -1 for fixed spherical stars. (Default: 1.0)
    tll : float, optional
        Tidal lead/lag angle (degrees). (Default: 0.0)
    l3 : float, optional
        Third light in units where (L_A + L_B + L_C = 1) (Default: 0.0)
    pcf : float, optional
        Phase correction factor (i.e. phase of primary eclipse). (Default: 0.0)
    lsf : float, optional
        Light scaling factor in magnitudes. (Default: 0.0)
    irs : float, optional
        Integration ring size in degrees. (Default: 1.0)
    p_orb : float, optional
        Orbital period in days. (Default: 3.35)
    t_0 : float, optional
        Reference time of primary eclipse in HJD. (Default: 54833.0)
        BJD = 2454833.0 is the zero point of the Kepler data.
    ld2_a : float, optional
        Limb darkening coefficient 2 for star A. (Default: 0.2586)
        Default is from Sing et al. as described above.
    ld3_a : float, optional
        Limb darkening coefficient 3 for star A. (Default: 0.0)
    ld4_a : float, optional
        Limb darkening coefficient 4 for star A. (Default: 0.0)
    ld2_b : float, optional
        Limb darkening coefficient 2 for star B. (Default: 0.2586)
    ld3_b : float, optional
        Limb darkening coefficient 3 for star B. (Default: 0.0)
    ld4_b : float, optional
        Limb darkening coefficient 4 for star B. (Default: 0.0)
    rv_a : float, optional
        Radial velocity amplitude of star A in km/s. (Default: 0.0)
    rv_b : float, optional
        Radial velocity amplitude of star B in km/s. (Default: 0.0)
    sys_a : float, optional
        Systemic velocity of star A in km/s. (Default: 0.0)
    sys_b : float, optional
        Systemic velocity of star B in km/s. (Default: 0.0)

    Examples
    --------
    Simulate a 15.5 day period binary.
    Star A: 1% modulation with a period of 3.71 days.
    Star B: 1.5% modulation with a period of 1.21 days.

    >>> sine_star = (-1, -2)
    >>> sine_params = ((54830, 3.71, 0.01), (54831, 1.21, 0.015))
    >>> bop_it = BopIt(sfsr=0.2, p_orb=15.5, csbr=0.8,
    >>>                sine_star=sine_star, sine_params=sine_params)

    Make a light curve with some noise.
    >>> phase, time, flux = bop_it.light_it(sigma=0.005)

    """
    def __init__(self, sine_star=None, sine_params=None, ldtype=(4, 4),
                 **kwargs):

        self.ldtype = ldtype

        # An ordered dictionary of optional keywords for JKTEBOP.
        self.params = OrderedDict([('csbr', 1.0),
                                  ('sfsr', 0.1),
                                  ('rsr', 1.0),
                                  ('ld_a', 0.3985),
                                  ('ld_b', 0.3985),
                                  ('inc', 90.0),
                                  ('ecc', 0.0),
                                  ('omg_bar', 0.0),
                                  ('gd_a', 0.0),
                                  ('gd_b', 0.0),
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
                                  ('ld2_a', 0.2586),
                                  ('ld3_a', 0.0),
                                  ('ld4_a', 0.0),
                                  ('ld2_b', 0.2586),
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
                self.vv[(30 + idx * 3):(33 + idx * 3)] = sine_params

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
        phase : array_like
            Orbital phase.
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
        phase = ((time - self.params['t_0']) % self.params['p_orb']) / \
            self.params['p_orb']
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

        return phase, time, flux

    @staticmethod
    def _bop_it(start_time=6, repeats=5):
        """
        A silly Easter egg method that (poorly) emulates the game Bop It!

        Parameters
        ----------
        start_time : int, optional
            The starting value of the timeout in seconds. (Default: 6)
        repeats : int, optional
            The number of times to repeat a timeout value before decreasing
            by one second.

        Examples
        --------
        Type the same text as the on-screen prompt. Press Enter to start.

        Bop it!
        Bop it!

        Twist it!
        Twist it!

        Pull it!
        Twist it!

        Aaaahhh!!! You lose.

        """
        commands = ['Bop it!', 'Twist it!', 'Pull it!']

        def interrupted():
            raise StandardError

        signal.signal(signal.SIGALRM, interrupted)

        def text_input(correct):
            try:
                    print '\n{}'.format(correct)
                    foo = raw_input()
                    if foo == correct:
                        return True
                    else:
                        return False
            except StandardError:
                    return False

        alarm_times = np.arange(repeats, repeats * (start_time + 1))[::-1] / \
            repeats

        raw_input('\n Type the same text as the on-screen prompt.'
                  'Press Enter to start.')

        lost = False
        for alarm_time in alarm_times:
            draw = random.choice(commands)
            signal.alarm(alarm_time)
            okay = text_input(draw)

            if not okay:
                print "\n\n Aaaahhh!!! You lose.\n"
                lost = True
                break

        if not lost:
            print 'You win!'
