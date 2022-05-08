#!/usr/bin/env python3

import sys
import math
import time
import numpy as np
import matplotlib.pyplot as plt
from au.arteritmic import convolution



def outliers(signal, n=1, cutoff=5e-4, plot=False):
    """
    This function can be used to locate bad data points using a "moving median" filter. For the median filter, instead of deleting bad data these are replaced by a median value.
    @param [ary1D] s      : Input signal
    @param [int  ] n      : Car-box size for colvolution utility
    @param [float] cutoff : Threshold for a bad data point (optional)
    @param [bold ] plot   : Plot result with option "true" (optional)
    return [ary1D] s_new  : Data corrected for bad data
    """
    # Finding dif:
    s_med = tt.moving('median', s, n)
    dif0  = s/s_med - 1
    
    # Replace median signal if outside cutoff region:
    above    = np.where(dif0>cutoff)[:][:]
    S[above] = S_med[above]
    below    = np.where(dif0<-cutoff)[:][:]
    S[below] = S_med[below]

    # Consistency check for median replacement:
    S_med = tt.moving('median', S, n)
    dif1  = S/S_med - 1
    
    # Plot outliers:
    if plot==1:
        pt.plot_locate(t, dif0, dif1, cutoff, n, 'Time (days)', 'dif')

    # Output:
    return np.vstack([t, S]).T



def jumps(data, gapsize, plot=None):
    """
    This function corrects for jumps in the data larger than "gapsize". Normally ajacent datapoint to a jump is effected, hence, a specified number of points of each side of a jump can be removed.
    ------------INPUT:
    data             : Data containing [time, signal].
    gapsize          : Signal difference at which there will be considered as a jump.
    plot             : Plots the original and corrected data. If plot==1 a plot is made (optional)
    -----------OUTPUT:
    [t, S]           : Data corrected for jumps."""
    print('-------------------------- jumps')
    
    # Unfold the data:
    t      = data[:,0]
    S      = data[:,1].copy()

    # Find distances/difference between data points:
    S_diff = np.diff(S)
    
    # Find gaps:
    index = np.where(abs(S_diff)>gapsize)[0][:]
    
    # Move the data when a jump:
    for i in index:
        S[i + 1:] -= S_diff[i]
    
    # Plot outliers:
    if plot==1:
        pt.plot_gapsize(S_diff, 'diff', 'log(N  of  diff)')
        pt.plot_jumps(t, data[:,1], S, 'Time (days)', 'Signal')

    # Output:
    return np.vstack([t, S]).T
