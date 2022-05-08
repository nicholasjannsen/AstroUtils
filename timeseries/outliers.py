#!/usr/bin/env python3

import numpy as np
from au.operations import convolution
from au.plots import 


def outliers(s, n=2, threshold=5e-4, plot=False):
    """
    This function can be used to locate bad data points using a data convolution.
    For the median filter, instead of deleting bad data these are replaced by a median value.
    @param [ary1D] s      : Input signal
    @param [int  ] n      : Car-box size for colvolution utility
    @param [float] cutoff : Threshold for a bad data point (optional)
    @param [bold ] plot   : Plot result with option "true" (optional)
    return [ary1D] s_new  : Data corrected for bad data
    """
    # Finding dif:
    s_med = convolution('median', s, n)
    dif0  = s/s_med - 1

    # Replace median signal if outside cutoff region:
    above    = np.where(dif0 > threshold)[:][:]
    s[above] = s_med[above]
    below    = np.where(dif0 < -threshold)[:][:]
    s[below] = s_med[below]

    # Consistency check for median replacement:
    s_med = convolution('median', s, n)
    dif1  = s/s_med - 1

    # Plot outliers:
    if plot is True:
        pt.plot_locate(t, dif0, dif1, cutoff, n, 'Time (days)', 'dif')

    # Output:
    return np.vstack([t, S]).T

