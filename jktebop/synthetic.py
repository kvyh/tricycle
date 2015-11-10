#!/usr/bin/env python
# -*- coding: utf-8 -*-
from collections import OrderedDict
import numpy as np

class BOPit(object):
    """

    """
    def __init__(self, sines=None, sine_params=None, ldtyp_a=1, ldtype_b=1,
                 **kwargs):
        self.params = OrderedDict([('csbr', 1.0),
                                  ('sfr', 0.1),
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

        if sines is not None and sine_params is not None:

            if len(sines) is not len(sine_params):
                raise TypeError('Number of sines does not equal number of '
                                'parameter sets.')

            # Add sine parameters to parameter array.
            for idx, sine_params in enumerate(sine_params):
                self.vv[30 + idx * 3:33 + idx * 3] = sine_params

            self.psine[0:len(sines)] = sines
            self.nsine = len(sines)

        #

        # Additional filler for JKTEBOP
        vary = np.zeros(138, dtype=float)
        npoly = 0
        ppoly = np.zeros(9)
        dtype = 1

    def getmodel_args(self):
        pass
    
    def light_curve(length=90.0, sc=False):
        pass