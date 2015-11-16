#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

import pandas as pd
import numpy as np
from scipy.interpolate import interp1d


def phase_fit(phase, flux, order=2, iters=1000, step=0.01, knots=None,
              find_knots=True, find_step=True):
    """
    Run polyfit on a phased lightcurve to find the eclipses.

    This is a somewhat ad hoc Python wrapper around the C code.
    http://phoebe-project.org/1.0/?q=node/103

    Polyfit attempts to divide the lightcurve into four segments and fit
    polynomials to each. The segments are the primary and secondary eclipse,
    and the phases in between.

    Parameters
    ----------
    phase : array_like
        The input orbital phase.
    flux : array_light
        The input fluxes in relative units.
    order : int, optional
        The order of the polynomial fits. (Default: 2)
    iters : int, optional
        The number of iterations. (Default: 1000)
    step : float, optional
        Step size for random knot displacement.
    knots : str, optional
        An explicit list of knots.
    find_knots : bool, optional
        Default is to find knots automatically.
    find_step : bool, optional
        Default is to find step size automatically.

    Returns
    -------
    params : 6-sequence
        The parameters of the eclipses.
        p_width: width of primary eclipse in phase.
        s_width: width of secondary eclipse.
        p_depth: depth of primary eclipse in relative flux.
        s_depth: depth of secondary eclipse.
        sep: separation between primary and secondary eclipse in phase.
        p_phase: Phase of the primary eclipse.

    References
    ----------
    Prša et al. (2008), ApJ 687, 542

    """
    # Polyfit expects phase from -0.5 to 0.5
    phase = phase - 0.5

    # Parse arguments to construct polyfit command.
    command = "polyfit -o {} -i {}".format(order, iters)

    if find_step:
        command += " --find-step"
    else:
        command += " -s {}".format(step)

    if find_knots:
        command += " --find-knots"
    elif knots is not None:
        command += " -k {}".format(knots)
    elif knots is None:
        raise ValueError("Use 'knots' to provide an explict list of knots.")

    command += " lc.dat > lc.out"

    # Write light curve to text file line by line.
    text_file = open("lc.dat", "w")

    for ph, fl in zip(phase, flux):
        text_file.write("{:3.6f} {:3.6f} \n".format(ph, fl))

    text_file.close()

    os.system(command)

    df = pd.read_table('lc.out', delim_whitespace=True, comment='#',
                       names=('phase', 'flux'))

    # We want phase from 0 to 1
    df.phase += 0.5

    # Interpolate the fit light curve to find flux at mid-eclipse.
    interpolated_lc = interp1d(df.phase, df.flux)

    ff = open("lc.out")
    lines = ff.readlines()

    widths = np.empty(4)
    midpoints = np.empty(4)
    midfluxes = np.empty(4)

    # FIXME: This pulls parameters from specific lines and columns in file.
    # This is stable as long as the output of polyfit.c doesn't change,
    # but a more robust solution is probably desirable.
    for ii, row in enumerate([19, 20, 21, 22]):
        ph1 = float(lines[row][18:24]) + 0.5
        ph2 = float(lines[row][26:32]) + 0.5

        # Egress should be at greater phase than ingress.
        if ph2 < ph1:
            ph2 += 1.0

        mid = (ph2 + ph1) / 2.0
        # Midpoint should be between 0 and 1.
        if mid >= 1.0:
            mid -= 1.0

        widths[ii] = ph2 - ph1
        midpoints[ii] = mid
        midfluxes[ii] = interpolated_lc(mid)

    # Sort by flux at midpoint to find primary and secondary eclipses.
    sort_indices = np.argsort(midfluxes)
    p_width, s_width = widths[sort_indices][:2]
    p_depth, s_depth = midfluxes[sort_indices][:2]
    sep = np.diff(midpoints[sort_indices][:2])[0]
    p_phase = midpoints[sort_indices][0]

    # Separation should be > 0.
    if sep < 0:
        sep += 1.0

    params = (p_width, s_width, p_depth, s_depth, sep, p_phase)

    os.remove('lc.dat')
    os.remove('lc.out')

    return params, df.phase, df.flux
