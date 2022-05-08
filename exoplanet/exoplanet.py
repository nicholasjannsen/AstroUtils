#!/usr/bin/env python3

import sys
import math
import time
import numpy as np
import matplotlib.pyplot as plt


def model(t, P, phi, dT, A):
    """
    Simple exoplanet transit model.
    @param [ary] t     : Time series with time (t) in [days] and signal (S)
    @param [flt] P     : Period
    @param [flt] phi   : Phase
    @param [flt] dT    : Transit duration
    @param [flt] A     : Scale amplitude
    return [ary] Model : Exoplanet model
    """
    # Model created:
    Model = np.zeros(len(t))
    N_t   = np.ceil(t[-1]/P)   # Max number of transits
    for m in np.arange(N_t+1):
        model        = np.abs(t - m*P - phi) < dT/2.
        Model[model] = A
    return Model


def cc_coefficient(x, y):
    """
    Function to find the cross-correlation coefficient between two datasets.
    @param [ary] x : Signal from 1st dataset
    @param [ary] y : Signal from 2nd dataset
    return [ary] r : Correlation coefficient
    """
    cor  = np.sum( (x-np.mean(x)) * (y-np.mean(y)) )
    norm = sqrt( np.sum((x-np.mean(x))**2) * np.sum((x-np.mean(x))**2) )
    r    = cor/norm
    return r


def crosscor(data, P_int, phi_int=None, model_const=None, save=None): 
    """
    This function is used to create a transit model for the Cross-Correlation (CC). To create the model the subrutine called 'model' is used. To perform the CC the CC coefficient is also needed and this is calculated in the subroutine 'cc_coefficients'.
    @param [ary2D] data           : Time series with time (t) in [days] and signal (S).
    @param [flout] P_int          : Periode interval to perform CC [hours].
    @param [float] phi_int        : Phase interval to perform CC [hours] (optional).
    @param [float] model_const    : Settings for model [dT, dP, dPhi] all in [min] (optional).
    @param [bold ] save           : If 1 the correlation is saved to data file with user defined name (optional)
    """

    start_time = time.time()      # Take time

    # Load functions:
    from Timeseries_Tools import model, cc_coefficient
    from Plot_Tools import SURF, plot_period_cc

    # Unpack data:
    t = data[:,0]*24.*60.       # Time [min]
    S = data[:,1].copy()        # Save copy
 
    # Normalized signal MUST be used:
    x = (1-S)-mean(1-S)
 
    # Model parameters:
    if model_const==None:
        dT   = 120.                 # Transit duration  [min]
        dP   = 5.                   # Period resolution [min]
        dphi = 30.                  # Phase  resolution [min]
        A    = 1.                   # Scale amplitude
    else:
        dT   = model_const[0]       # [min]
        dP   = model_const[1]       # [min]
        dphi = model_const[2]       # [min]
        A    = model_const[3] 

    # Parameter space:
    P_int   = array(P_int)*60.             # Period interval [min]
    if phi_int==None:
        phi_int = [0, int(P_int[1]+dphi)]  # Phase interval [min] (Default: all possible)
    else:
        phi_int = array(phi_int)*60.       # Phase interval  [min]

    # Tested period-grid:
    P     = arange(P_int[0], P_int[1], dP)           # Period range [min]
    phi   = arange(phi_int[0], phi_int[1], dphi)     # Phase  range [min]

    # Dimentions:
    N = len(P) 
    M = len(phi) 
    print 'CC matrix = N(P) x M(phi) = {} x {} = {}'.format(N, M, N*M)

    # Preparations for cross-correlation:
    cc    = zeros((N, M))                           # N*M matrix with zeros for amplitudes
    Model = zeros(len(t))                           # 0-array for model data     

    # Range of highest number of transits given by P_int:
    N_max = range(1, int(t[-1]/P_int[0])+2)

    # Perform Cross-correlation:
    for i in range(N):
        for j in range(M):
            y        = model(t, P[i], phi[j], dT, A)
            r_cc     = cc_coefficient(x, y)
            cc[i, j] = r_cc
        # Keep track on time:
        if i==round(N*0.25):
            print ('Done 25  procent --- %s seconds ---' % (time.time() - start_time))       
        if i==round(N*0.50):
            print ('Done 50  procent --- %s seconds ---' % (time.time() - start_time))       
        if i==round(N*0.75):
            print ('Done 75  procent --- %s seconds ---' % (time.time() - start_time))       
        if i==round(N-1):
            print ('Done 100 procent --- %s seconds ---' % (time.time() - start_time))       

    # Best P and Phi by cross-correlation:
    cc_max   = max(cc)
    cc_max_i = where(cc==cc_max)
    h = 60.;     P_hour = P[cc_max_i[0][0]]/h; phi_hour = phi[cc_max_i[1][0]]/h
    d = 60.*24.; P_days = P[cc_max_i[0][0]]/d; phi_days = phi[cc_max_i[1][0]]/d
    print 'Best Period: {:.6f} hours and {:.6f} days'.format(P_hour, P_days)
    print 'Best Phase : {:.6f} hours and {:.6f} days'.format(phi_hour, phi_days)
    print 'P resolution  : {:.4f}'.format(dP/(60.*24.))
    print 'phi resolution: {:.4f}'.format(dphi/(60.*24.))

    # Save result?
    if save!=None: 
        savetxt(save[0], P)
        savetxt(save[1], phi)
        savetxt(save[2], cc)

    return


def autocor(data, npeaks, cutoff_x, cutoff_y, plot=None):
    """----------------------------------- FUNCTION ---------------------------------------:    
    # This function performs a auto-correlation on the data signal.
    #-----------INPUT:
    # data           : Time series time [mins] and signal. 
    # npeaks         : Number of peaks to be found in ACF by the 'argrelmax' function from Scipy.
    # cutoff_x       : Time range where peaks should be used.
    # cutoff_y       : ACF range where peaks should be used.
    # plot           : If plot==1 the data is plotted."""
    print '------------------------------------------------------------------------- autocor'

    # Load functions:
    from scipy.interpolate import griddata     # Linear interpolation
    from scipy.signal import argrelmax         # Peak location
    from scipy import stats                    # Uncertainty to linear fit
    from Plot_Tools import PLOT

    # Avoid overwritting data:
    t = data[:,0].copy()*24*60         # [min]
    S = 1 - data[:,1].copy()    # Normalized signal MUST be used

    # Plot corrected timeseries:
    if plot==1:
        data0 = vstack([t/(60.*24.), S]).T
        PLOT([data0], 'b-', '$t$ [days]', '$1-S$', 'Correected data')

    # Interpolate to uniform grid:
    dt = median(diff(t))
    tt = arange(t[0], t[-1]+dt, dt)
    SS = griddata(t, S, tt, method='nearest')   # Function from scipy libiary

    # Remove points below 0:
    SS[SS<0] = 0.0

    # Plot interpolated data:
    if plot==1:
        data1 = vstack([tt/(60.*24.), SS]).T
        PLOT([data0, data1], ['b-','k.'], '$t$ [days]', '$1-S_{grid}$', 'Interpolation')

    # Prepare auto-correlation:
    N   = len(tt)
    acf = zeros(N)

    # Perform auto-correlation:
    for i in range(1, N):
        S_stat  = SS[i:]  - mean(SS[i:])                   # Stationary grid to be correlated
        S_move  = SS[:-i] - mean(SS[:-i])                  # Moving the stationary grid by i
        acf[i]  = sum(S_stat*S_move)/(float(N)-float(i))   # Correlates the two grides    
    # Auto-correlation for i=0:
    S_stat  = SS - mean(SS)
    acf[0]  = (1./float(N))*sum(S_stat**2)                 # Correlates the two grides

    # Find peaks in ACF and find location and height:
    peaks_i = argrelmax(acf, order=npeaks)
    peaks_x = tt[peaks_i]
    peaks_y = acf[peaks_i]

    # Peaks inside y-range:
    peaks_keep = peaks_y > cutoff_y[0]
    peaks_x    = peaks_x[peaks_keep]
    peaks_y    = peaks_y[peaks_keep]
    peaks_keep = peaks_y < cutoff_y[1]
    peaks_x    = peaks_x[peaks_keep]
    peaks_y    = peaks_y[peaks_keep]
    peak0      = peaks_x[0]

    # Peaks inside x-range:
    peaks_keep = peaks_x > cutoff_x[0]*24.*60.  # cutoff_y is given in days
    peaks_x    = peaks_x[peaks_keep]
    peaks_y    = peaks_y[peaks_keep]
    peaks_keep = peaks_x < cutoff_x[1]*24.*60. 
    peaks_x    = peaks_x[peaks_keep]
    peaks_y    = peaks_y[peaks_keep]

    # Plot autocorrelation:
    if plot==1:
        data0 = vstack([tt/(60.*24.), acf]).T
        data1 = vstack([peaks_x/(60.*24.), peaks_y]).T
        PLOT([data0, data1], ['b-', 'k.'], '$t_{grid}$ [days]', '$ACF$', 'Auto-correlation')

    # Measured periods:
    P_mea = peaks_x - t[0]          # Remove time-offset
    P     = mean(diff(peaks_x))     # Mean period value

    # This is an comtempolary solution: When the cutoff_x[0]>0 the linear plot do no work since
    # it needs the 1. detected transit, hence, one have to count N*P to the first peak that's
    # inside the range cutoff_x. Here 6 period is evident before the first inside cutoff_x:
    #P_mea = P_mea - 6*P

    # Uncertainty:
    P_std = std(diff(peaks_x))

    # Calculate periods:
    P_cal = zeros(len(P_mea))
    for i in range(1, len(P_mea)+1):
        P_cal[i-1] = i*P-P

    # Find period by linear regression:
    coef, stats = np.polynomial.polynomial.polyfit(P_cal, P_mea, 1, full=True)
    b     = coef[0]
    a     = coef[1]
    p     = arange(0., P_cal[-1], 1.)
    P_fit = a*p + b

    # Write out the estimated period:
    print 'Best: P = {} +/- {} hours'.format(b/60., P_std/60.)
    print 'Best: P = {} +/- {} days'.format(b/(60.*24.), P_std/(60.*24.))
    print 'Mean: P = {} +/- {} days'.format(P/(60.*24.), P_std/(60.*24.))

    # Plot linear fit:
    if plot==1:
        data0 = vstack([P_cal/(60.*24.), P_mea/(60.*24.)]).T
        data1 = vstack([p/(60.*24.), P_fit/(60.*24.)]).T
        PLOT([data0, data1], ['go','k-'], \
             r'Calculated: $N \times P-P$ [days]', r'Measured: $N \times P$ [days]', \
             'ACF Period: ${:.4f} \pm {:.4f}$ days'.format(b/(60.*24.), P_std/(60.*24.)))

    return
