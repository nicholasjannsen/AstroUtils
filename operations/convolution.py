#!/usr/bin/env python3

import numpy as np


def convolution(signalIn, convType, numBox):
    """
    This utility makes the proper convolution to a signal dataset.
    Notice: the carbox size here is twice what is default by numpy.

    @param [str] convType   : Filter is either: median or mean
    @param [ary] signalIn   : Signal needed for processing
    @param [int] numBox     : Integer used as car-box size for convolution

    @return [ary] signalOut : Convolved signal
    """

    # Constants
    n     = numBox
    S     = signalIn.copy()   # Avoid overwritting the input signal
    S_new = np.zeros(len(S))
    nzero = np.zeros(2*n+1)   # Optimization constant

    # Select convolution
    if convType == 'log':     conv = np.log10
    if convType == 'mean':    conv = np.mean
    if convType == 'median':  conv = np.median
    if convType == 'sin':     conv = np.sin
    if convType == 'arcsin':  conv = np.arcsin
    if convType == 'arcsinh': conv = np.arcsinh
    if convType == 'cos':     conv = np.cos
    if convType == 'arcsin':  conv = np.arcsin
    if convType == 'arcsinh': conv = np.arcsinh

    # Interval: d[n, 1+n, ... , N-1, N-n]
    for i in range(len(S)-2*n):
        S_new[n+i] = conv(S[range((n+i)-n, (n+i)+n+1)])

    for i in range(n):
        # Interval: d[-n, -(n-1), ... , n-1, n] - Low end of data
        low = nzero
        low[range(n-i)] = S[0]*np.ones(n-i)
        low[-(n+1+i):]  = S[range(0, n+1+i)]
        S_new[i]        = conv(low)

        # Interval: d[N-n, N-(n-1), ... , N+(n-1), N+n] - High end of data
        high = nzero
        high[range(n+1+i)] = S[range(len(S)-(n+i+1), len(S))]
        high[-(n-i):]      = S[-1]*np.ones(n-i)
        S_new[len(S)-1-i]  = conv(high)

    # Output
    signalOut = S_new
    return signalOut
