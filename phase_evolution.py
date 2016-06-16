import matplotlib.pyplot as plt
import numpy as np

from scipy import stats, interpolate


def detrend(time, flux, quarter, order=3):
    """
    Detrend each quarter separately using polynomial fit.
    Parameters
    ----------
    time: array_like
        Observation times.
    flux: array_like
        Fluxes.
    quarter: array_like
        Kepler quarter.
    order: int, optional
        The order of the polynomial fit. Default: 3
    Returns
    -------
    flux_detrended : numpy.ndarray
        The detrended, median normalized flux.
    """
    # "Empty" array to hold detrended fluxes
    flux_detrended = np.zeros_like(flux)

    for q in np.unique(quarter):
        # Only select data in given quarter
        mask = quarter == q

        # Compute polynomial fit
        poly = np.polyfit(time[mask], flux[mask], order)
        z = np.poly1d(poly)

        # Subtract fit and median normalize
        flux_detrended[mask] = (flux[mask] - z(time[mask])) / \
            np.median(flux[mask])

    return flux_detrended


def median_lc_corr(time, flux, p_fold, t0=0., delta_phase=0.01, cad_min=3,
                   plot=False):
    """
    The cross-correlation with the median phase-folded light curve.
    For a given period, compute the binned, median phase-folded light curve.
    Then compute the cross-correlation with each successive cycle of light
    curve at that period.
    Parameters
    ----------
    time: array_like
        Observation times.
    flux: array_like
        Fluxes.
    p_fold: float
        The period at which to fold the light curve.
    t0: float, optional
        The reference time, e.g., time of primary eclipse. Default: 0.
    delta_phase: float, optional
        The phase bin width. Default: 0.01
    cad_min: int, optional
        Exclude light curve sections with fewer cadences than `cad_min`.
    plot : bool, optional
        Set to True to plot phase-folded light curve and cross-correlation
    Returns
    -------
    cycle_num : numpy.ndarray
        Cycle number.
    corr : numpy.ndarray
        The cross-correlation
    """
    # Calculate the phase.
    phase = ((time - t0) % p_fold) / p_fold
    # Calculate the cycle number.
    cycle = ((time - t0) // p_fold).astype(int)
    # Start at zero
    cycle -= cycle.min()

    # Computed binned median
    bins = np.arange(0., 1. + delta_phase, delta_phase)
    binned_med, bin_edges = stats.binned_statistic(phase, flux,
                                                   statistic="median",
                                                   bins=bins)[:2]

    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    # Wrap around at the beginning and end of the binned light curve
    # using the fact that f(phase) = f(phase + 1)
    wrap = ([bin_centers[-1] - 1.], bin_centers, [1. + bin_centers[0]])
    bin_centers = np.concatenate(wrap)
    wrap = ([binned_med[-1]], binned_med, [binned_med[0]])
    binned_med = np.concatenate(wrap)

    # Linear interpolation of binned light curve.
    interp = interpolate.interp1d(bin_centers, binned_med)

    # Only use cycles with more cadences than `cad_min`.
    cycle_num = np.arange(cycle.max() + 1)[np.bincount(cycle) > cad_min]

    # Empty array to hold cross-correlation
    corr = np.zeros_like(cycle_num, dtype=float)

    for i, n in enumerate(cycle_num):

        mask = cycle == n

        p = phase[mask]
        f = flux[mask]
        f_i = interp(p)

        # Normalize light curves
        a = (f - np.mean(f)) / (np.std(f) * len(f))
        v = (f_i - np.mean(f_i)) / np.std(f_i)

        corr[i] = np.correlate(a, v)

    if plot:

        fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10, 5))

        ax1.scatter(phase, flux, color="k", s=0.1)
        p = np.linspace(0, 1, 1000)
        ax1.plot(p, interp(p), color="r")

        ax1.set_xlim(0, 1)
        ymax = np.percentile(np.abs(flux), 98)
        ax1.set_ylim(-ymax, ymax)
        ax1.set_xlabel("Phase")
        ax1.set_ylabel("Normalized Flux")

        ax2.plot(cycle_num * p_fold, corr, "-k", markersize=3)
        ax2.set_ylim(-1.1, 1.1)
        ax2.set_xlabel("Time")
        ax2.set_ylabel("Normalized Cross-correlation")

        plt.show()

    return cycle_num, corr