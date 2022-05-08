#!/usr/bin/env python3

import sys
import math
import time
import numpy as np
import matplotlib.pyplot as plt
import Timeseries_Tools as tt
import Plot_Tools as pt


def power(data, f_interval=None, f_resolution=None, sampling=None, w_column=None):
    """
    This function calculate the power spectrum of a time varying data set.

    data           : Is a matrix with a column of time, signal, and possible weights.
    f_interval     : Frequency interval containing [fmin, fmax] (optional).
    f_resolution   : Frequency resolution (optional).
    sampling       : Sampling (oversampling >1) of the data (optional).
    w_column       : Column number of weight (optional). w_column = 0 (calculates without weights)

    Pf_data        : Frequencies (f) and power (P).
    P_comp         : Components for defining P (alpha and beta).
    f_interval     : Frequency interval may be usefull later.
    f_resolution   : Frequency resolution is handy when is it not specified
    """

    # Checking data structure:
    if ( np.size(data[0,:])<2   and                  # 2 columns minimum
         np.size(data[:,0])<2   and                  # 2 rows minimum 
         isinstance(data,float)):                    # Text
         sys.exit('Error: Wrong data structure')     # Abort

    # Replacing nan by 0:
    if np.sum(np.sum(data))==np.nan:
        data[data==nan]=0
    # Replacing inf by 0:
    if np.sum(np.sum(data))==np.inf or np.sum(np.sum(data))==-np.inf:
        data[data== np.inf]=0
        data[data==-np.inf]=0
       
    # Splitting data:
    t = data[:,0]                              # [days]
    S = data[:,1]-np.mean(np.array(data[:,1]))       # [arbt. unit] Correcting signal

    # Handy definitions:
    dt_min = np.min(np.diff(t))     # Minimum time interval
    dt_max = np.max(np.diff(t))     # Maximum time interval
    dt_med = np.median(np.diff(t))  # Median time interval
    f_Nq = 1./(2*dt_min)            # Cyclic Nyquist frequency
    f_md = 1./(2*dt_med)            # Cyclic median frequency
    
    # Checking if data is ordered in time:
    if dt_min<0:
        sys.exit('Error: Data not ordered in time')

    # Checking the existance of f_interval and f_resolution:
    # If f_interval and f_resolution is not defined:
    if f_interval==None and f_resolution==None:
        if dt_max>2*dt_med or math.isnan(f_Nq):
            # print 'Using f_md         = %f c/d' %f_md
            f_interval = [0, f_md]
            f_resolution = 1./(t[-1]-t[0])
        else:
            # print 'Using f_Nq         = %f c/d' %f_Nq
            f_interval = [0, f_Nq]
            f_resolution = 1./(t[-1]-t[0])
    # If only f_resolution is defined:
    elif f_interval==None and f_resolution!=None and isinstance(f_resolution, float): 
        f_resolution = f_resolution
        if dt_max>2*dt_med or math.isnan(f_Nq):
            # print 'Using f_md         = %f c/d' %f_md
            f_interval = [0, f_md]
        else:
            # print 'Using f_Nq         = %f c/d' %f_Nq
            f_interval = [0, f_Nq]
    # If only f_interval is defined:
    elif f_resolution==None and len(f_interval)==2:
        f_interval   = f_interval
        f_resolution = 1./(t[-1]-t[0])
    # If both f_interval and f_resolution is defined:
    elif f_interval!=None and f_resolution!=None:
        f_interval   = f_interval
        f_resolution = f_resolution
    else:
        sys.exit('Error: Wrong input argument for "f_interval"')
    # print 'Using f_interval   = {} c/d'.format(f_interval)    # Use .format when having tuples
    # print 'Using f_resolution = {} c/d'.format(f_resolution)
   
    # The frequency interval is found:
    if sampling==None: 
        f = np.arange(f_interval[0], f_interval[1], f_resolution)
    elif sampling!=None:
        # print 'Using sampling     = %f' %sampling
        f = np.arange(f_interval[0], f_interval[1], f_resolution*sampling)

    # The power spectrum is calculated:
    s  = np.zeros(len(f))
    c  = np.zeros(len(f))
    ss = np.zeros(len(f)) 
    cc = np.zeros(len(f))
    sc = np.zeros(len(f))
    picon = 2*np.pi*t
    
    # Without weights
    if w_column==0 or np.size(data[0,:])==2:
        # print 'Calculating WITHOUT weights'
        for i in range(len(f)):
            f_i = f[i]
            sinsum = np.sin(f_i*picon)
            cossum = np.cos(f_i*picon)
            s[i]  = np.sum(S*sinsum)
            c[i]  = np.sum(S*cossum)
            ss[i] = np.sum(sinsum**2)
            cc[i] = np.sum(cossum**2)
            sc[i] = np.sum(sinsum*cossum)
            # Print compilation time to bash:
            pt.compilation(i, len(f), 'power: without weights')
        print # This function is needed for compilation 

    # Calculate with weights:
    if np.size(data[0,:])>2:
        # print 'Calculating WITH weights'
        if w_column==None:
            w = data[:,2]
        elif w_column!=None:
            w = data[:,w_column]
        for i in range(len(f)):
            f_i = f[i]
            sinsum = np.sin(f_i*picon)
            cossum = np.cos(f_i*picon)
            s[i]  = np.sum(w*S*sinsum)
            c[i]  = np.sum(w*S*cossum)
            ss[i] = np.sum(w*sinsum**2)
            cc[i] = np.sum(w*cossum**2)
            sc[i] = np.sum(w*sinsum*cossum)
            # Print compilation time to bash:
            pt.compilation(i, len(f), 'power: with weights')
        print 

    # Calculate alpha, beta, P, A
    alpha = (s*cc-c*sc)/(ss*cc-sc**2)
    beta  = (c*ss-s*sc)/(ss*cc-sc**2)
    P     = alpha**2+beta**2
    A     = np.sqrt(P)

    # Output:
    Pf_data = np.vstack([f, P]).T 
    P_comp  = np.vstack([alpha, beta]).T
    return Pf_data, P_comp, f_interval, f_resolution  



def window(data, f_interval=None, f_resolution=None, sampling=None, w_column=None):
    """
    This function calculate the so-called window function. The function calls the routine 'power'.
    It returns the power spectrum of the window function.

    data         : Is a matrix with a column of time, signal, and possible weights
    f_interval    : Frequency interval [fmin, fmax] (optional).
    f_resolution  : Frequency resolution (optional).
    oversampling  : Oversampling of the data (optional).
    w_column      : Column the weights are placed (optional).

    Pf_window     : Frequency and power of the window function.
    """

    # Avoid overwritting data:
    data0 = data.copy()

    f_range = round(f_interval[0]+(f_interval[1]-f_interval[0])/2)
    picon   = 2*np.pi*f_range*data[:,0]
    fsin    = np.sin(picon)
    fcos    = np.cos(picon)

    # Sinusoidal
    data0[:,1] = fsin
    Pf_power, _, _, _, = tt.power(data0, f_interval, f_resolution, sampling, w_column)
    f    = Pf_power[:,0]
    Psin = Pf_power[:,1]

    # Co-sinusoidal
    data0[:,1] = fcos
    Pf_power, _, _, _, = tt.power(data0, f_interval, f_resolution, sampling, w_column)
    f    = Pf_power[:,0]
    Pcos = Pf_power[:,1]

    # Output:
    P = 1./2*(Pcos+Psin)
    Pf_window = np.vstack([f, P]).T
    return Pf_window



def clean(data, N_peaks, f_interval=None, f_resolution=None, sampling=None, w_column=None):
    """
    This function indentify and select a sepcified number of highest valued peaks. The peaks are determined with a high accuracy and is then subtracted from the times series. As a output the routine returns this more "clean" time series.  
    -----------INPUT:
    data            : Is a matrix with a column of time, signal, and possible weights. 
    f_interval      : Frequency interval containing [fmin, fmax] (optional).
    f_resolution    : Frequency resolution (optional).
    sampling        : Sampling (oversampling >1) of the data (optional).
    w_column        : Column where the weights are placed (optional).
    ----------OUTPUT:
    St_clean        : Times (t) and a cleaned Signal (S).
    P_comp          : Components of P is alpha and beta.
    f_peaks         : Highest frequency peak value for subtracted peaks
    """

    # Avoid overwritting data:
    data0 = data.copy()

    # Standard frequency resolution:
    T = data0[-1,0]-data[0,0]
    if f_resolution==None:
        f_resolution = 1/T

    # Avoid 0 as input as not peaks are found:
    if f_interval[0]==0:
        f_interval = [f_resolution, f_interval[1]]

    # Constants:
    SAMPLING = 1
    f_RES    = 0.1*f_resolution     # Standard frequency resolution
    picon    = 2*np.pi*data0[:,0]      # Optimization constant
    f_peaks  = np.zeros(N_peaks)
    A_peaks  = np.zeros(N_peaks)

    for i in range(N_peaks):
        k = i+1
        print '%s. Peak' %k

        # 1. Iteration - start finding largest peak:
        Pf_power, _, _, _, = tt.power(data0, f_interval, f_resolution, sampling, w_column)
        f     = Pf_power[:,0];   P = Pf_power[:,1];  j = np.nanargmax(P)
        f_int = (f[j-1], f[j+1]) # Smaller f_int (Tuple instead of array for optimization)

        # Testing that the frequency resolution > sigma_f to continue:
        A_peak    = P[j]
        A_av      = np.mean(np.sqrt(P))
        sigma_a   = 0.8*A_av
        sigma_phi = sigma_a/A_peak
        sigma_f   = np.sqrt(3)*sigma_phi/(np.pi*T)

        if f_RES>sigma_f:

            # 2. Iteration: uses now f_res and so on..
            Pf_power, _, _, _, = tt.power(data0, f_int, f_RES, SAMPLING, w_column)
            f     = Pf_power[:,0];   P = Pf_power[:,1];  j = np.nanargmax(P)
            f_int = (f[j-1], f[j+1])

            # 3. Iteration: last
            Pf_power, P_comp, _, _, = tt.power(data0, f_int, f_RES, SAMPLING, w_column)
            f = Pf_power[:,0];  P = Pf_power[:,1];  j = np.nanargmax(P)
            fpicon = picon*f[j]  # Optimization constant
            alpha  = P_comp[:,0];  beta = P_comp[:,1]
            alpha0 = alpha[j]*np.sin(fpicon)
            beta0  = beta[j]* np.cos(fpicon)
            data0[:,1] = data0[:,1] - alpha0 - beta0
            f_peaks[i] = f[j]
            A_peaks[i] = np.sqrt(P[j])

    # Output:
    St_clean = data0
    print f_peaks, A_peaks
    return St_clean, f_peaks, A_peaks



def filters(data, f_interval, f_resolution=None, sampling=None, w_column=None):
    """
    This function takes a data set and a frequency interval and calcuates the power spectrum within this interval. The software can be used as a low-pass, band-pass, or a high-pass filter, which i solely determined by the frequency interval.

    data           : Is a matrix with a column of time, signal, and possible weights
    f_interval     : Frequency interval containing [fmin, fmax].
    f_resolution   : Frequency resolution (optional).
    sampling       : Sampling (oversampling >1) of the data (optional).
    w_column       : Column the weights are placed (optional).

    St_low_band    : Time series for low/band-pass filter.
    St_high        : Time series for high-pass filter
    """

    # Avoid overwritting data:
    data0 = data.copy()

    # Avoid 0 as input as not peaks are found:
    if f_interval[0]==0:
    f_interval = [f_resolution, f_interval[1]]

    # Calculates power spectrum:
    Pf_power, P_comp, _, _, = tt.power(data0, f_interval, f_resolution, sampling, w_column)
    t     = data0[:,0]
    f     = Pf_power[:,0]
    alpha = P_comp[:,0] 
    beta  = P_comp[:,1]

    # Calculates P_filter:
    P_filter = np.zeros(len(t))
    fpicon = 2*np.pi*f                    # Optimization constant
    for i in range(len(t)):
        tfpicon     = fpicon*t[i]         # Optimization constant
        alpha_sin   = alpha*np.sin(tfpicon)
        beta_cos    = beta* np.cos(tfpicon)
        P_filter[i] = np.sum(alpha_sin + beta_cos)

    # Calculates window function:
    Pf_window = tt.window(data0, f_interval, f_resolution, sampling)
    P_window  = Pf_window[:,1]

    # Bandpass/Lowpass and Highpass filter:
    S_low_band  = P_filter/np.sum(P_window)
    S_high      = data0[:,1]-S_low_band
    St_low_band = np.vstack([t, S_low_band]).T
    St_high     = np.vstack([t, S_high]).T
    return St_low_band, St_high





def stellar_noise(data, f_int, sampling, N, plot=None):
    """
    This function clean the power spectrum for long period variations from the stellar host or other long term variations. 
    ------------INPUT:
    data             : Containing [time, signal].
    f_int            : Frequency interval that should be cleaned/filtered.
    sampling         : Sampling of the data used in clean.
    N                :  Numbers of peaks to be removed by clean.
    plot             : Plot power spectrum and timeseries if plot is 1 (optional).
    -----------OUTPUT:
    St               : Corrected time series [time, Signal]. """
    print('-------------------------- stellar noise')
    
    # Run functions:
    f_int0   = [0, 1]
    Pf_power, _, _, f_res = tt.power(data, f_int0, None, sampling)
    St, _, _,             = tt.clean(data, N,  f_int, f_res)
    Pf, _, _, _,          = tt.power(St, f_int0, None, sampling)

    # Plot outliers:
    if plot==1:
        fx_int = [0, 0.5]; fy_int = [0, 2e6]
        pt.plot_stellarnoise_power(Pf_power, Pf, fx_int, fy_int,'Frequency $(c/d)$','Power')
        pt.plot_stellarnoise_time(data, St, 'Time (days)', 'Signal')

    # Output:
    return St



def slowtrend(data, gapsize, n=25, m=25, jump=None, plot=None):
    """
    This function corrects for slow trends in the time series.
    -----------INPUT:
    data            : Data containing [time, signal].
    gapsize         : Signal difference at which there will be considered as a jump.
    n and m         : Integer used in moving median and mean filter. n=m=25 if None (optional).
    jump            : If jump==1 the data is filtered with a n=2 median filter 'm' points of each                       side around a jump in time that is bigger than 'gapsize' (optional).
    ----------OUTPUT:
    data_new        : Data corrected for slow trends."""
    print('-------------------------- slowtrend')

    # Splitting data:
    t = data[:,0]
    S = data[:,1].copy()
    
    # Median and mean filters:
    S_medi  = tt.moving('median', S,      n)
    S_medi2 = tt.moving('median', S,      2)
    S_mean  = tt.moving('mean'  , S_medi, m)

    # Use median filter m points of each side around jumps:
    if jump==1:
        t_diff = np.diff(t)                              # Find difference between data points
        index  = np.where(np.abs(t_diff)>gapsize)[0][:]  # Find gap indices
        mcon   = range(-m, 1+m)                          # Optimization constant range
        for i in index:
            k = mcon+(1+i)*np.ones(len(mcon))  
            for j in k:
                S_mean[int(j)] = S_medi2[int(j)]
    
    # Correction:
    S_new  = S/S_mean
 
    # Data:
    data1 = np.vstack([t, S_medi]).T 
    data2 = np.vstack([t, S_mean]).T
    data3 = np.vstack([t, S_new]).T
    
    # Plot outliers:
    if plot==1:
        pt.plot_slowtrend(data, data1, data2, data3, n, m, 'Time (days)', 'Signal')
        pt.PLOT(data3[:,0], [data3[:,1]], 'b-', '$t$ [days]', '$S$', 'Corrected Lightcurve')

    #Output:
    return data3

