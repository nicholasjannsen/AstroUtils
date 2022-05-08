#!/usr/bin/env python3

import sys
import math
import time
import numpy as np
import matplotlib.pyplot as plt
from au.arteritmic import convolution


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
