>"""
SOFTWARE DESCRIPTION:
---------------------

Written April 2018 -- Nicholas Jannsen
Typeset in Python 2.7

In short this software performs aperture photometry for every star that have a image coordinate asigned. The first coordinate needs to belong to the target star for which a final light curve is desired. However, the flux and magnitude is calculated for every star, and all stars can be used as reference to correct for seeing variations during the total observation. Such a corrections is often refered to 'field (aperture) photometry' as the combined field of stars is used to correct for weather conditions during the night. The main function that needs to be called is 'aperture_photometry'.

Input paramters:
---------------------
path        (string): Directory path to your data 
LF_name     (string): Name of the Light Frames (LF) without end-number extentions
coor_name   (string): Name of stellar coordinate textfile
coor_target (list)  : Coordinates of target star '[x, y]' (int) 
R           (list)  : Apertures [r (int), r_min (int), r_max (int)]
dm_int      (list)  : Seeing box to choose stars within mag range and lowest spread in magn variations:
                      [mag_min , mag_max , dm_min, dm_max] (int, float)
plot        (int)   : Plot if you like: 
                      'plot=1': shows only final corrected lightcurve
                      'plot=2': shows (1) fits image with stellar coordinates, (2) trumpet, (3) corrected
                      'plot=3': shows (1) SNR as function of aperture size (2) apertures ontop of fits image
save        (int)   : If 'save=1' a textfile named 'final_lightcurve.txt' with [time, mag].
---------------------

How to run the code:
--------------------
First 
Name of stellar coordinates. With this string the software look for a textfile where all the coordinates is within (x, y). The software 'starfinder' within this software can be used to find and save all stellar coordinates. This function takes the ............


Aperture size of stellar aperture 'r', and minimum and maximum sky aperture radius 'r_min' and 'r_max', respectively.

Name of the image data/fits files. The software is written such that if your images is named e.g. 'image2018-XX-XXTXX-XX-XX.fits' you can simply use 'image' (if no other data in the same folder is named exactly this).

 The software only uses circular apertures and to select a good aperture size for the target star just set 'plot=3' to plot the Signal-to-Noise Ratio (SNR) as a function of aperture radius. The best aperture is the usually the local maximum of this plot. This aperture is then also used for every reference star.




Packages needed to run the code:
--------------------------------
"""
from numpy import inf, nan, sin, cos, pi, sqrt, diff, std, diag, argmin, log10, meshgrid
from numpy import mean, median, nanargmax, zeros, ones, ceil, delete, shape, roll, nonzero
from numpy import arange, array, size, vstack, hstack, copy, loadtxt, where, savetxt, linspace, shape
from numpy import min, max, sum, float, round, int
import math, sys, time, glob
import numpy as np
# Functions:
import scipy 
import scipy.ndimage as snd
from matplotlib.colors import LogNorm
from scipy.misc import imsave
from astropy.io import fits
from astropy.time import Time
# Plots:
from astropy.io import fits
import matplotlib.pyplot as plt
from Plot_Tools import FITS, MAG, plot_trumpet


def aperture_photometry(path, LF_name, coor_name, coor_target, R, dm_int=None, plot=None, save=None):
    
    #start_time = time.time()
    #from Photometry import aperture, SNR, center_of_flux
 
    #--- DATA ---#
    LF_files0 = np.sort(glob.glob('{}{}*'.format(path, LF_name[0])))
    LF_files1 = np.sort(glob.glob('{}{}*'.format(path, LF_name[1])))
    LF0       = fits.getdata(str(LF_files0[0]))
    LF1       = fits.getdata(str(LF_files1[0]))
    
    coor_stars = loadtxt('{}{}.txt'.format(path, coor_name))
    x, y   = coor_stars[:,1], coor_stars[:,0] # Stellar coordinates
    xs, ys = coor_target[0],  coor_target[1]  # Target  coordinates
    n = len(LF_files1)                        # Number of images

    #--- CONSTANTS ---#
    gain  = 0.73          # Gain of camera: electrons pr ADU (ADU = counts from object)- ajust to camera!
    ron   = 3.3           # Read out noise - ajust to camera!
    con   = 25            # Magnitude constant
    r     = R[0]          # Aperture radius = R[1]          # Minimum aperture radius to measure sky flux
    r_max = R[2]          # Maximum aperture radius to measure sky flux
    N     = len(x)        # Number of stars (start)
    h, w  = shape(LF1)     # Image dimentions (height, width)
    
    #--- STELLAR COORDINATES ---#
    # Make target star the first star:
    d = zeros(N)
    for i in range(N): d[i] = sqrt((xs-x[i])**2+(ys-y[i])**2) # Simple pythagorean geometry
    index_min = argmin(d)    # Find target star by geometry
    x = roll(x, -index_min)  # Shift data to target star index
    y = roll(y, -index_min)
    # Check if stellar coordinates are closer than r_max:
    i_x1 = where(x<r_max)[0]; i_x2 = where(x>w-r_max)[0]
    i_y1 = where(y<r_max)[0]; i_y2 = where(y>h-r_max)[0]
    i_xy = hstack([i_x1, i_x2, i_y1, i_y2])
    # Discard these coordinates:
    X = delete(x, i_xy)
    Y = delete(y, i_xy)
    N = len(X)               # New number of stars
    # Plot if you like:
    if plot==2:
        FITS(LF1,'linear'); plt.plot(x,y,'ro',mfc='none'); plt.plot(x[0],y[0],'b+');\
        plt.plot(x[i_xy], y[i_xy],'y*'); plt.show()

    #--- FIND OPTIMAL APERTURE RADIUS ---#
    if plot==3:
        SNR_r = optimal_aperture(LF1, X[0], Y[0], R, gain, ron, plot)
    
    #--- PHOTOMETRY ---#
    # Find fluxes:
    flux_star = zeros((n,N))
    SNR_i     = zeros((n,N))
    time      = zeros(n) 
    for i in range(n): # Loop over all images:
        hdulist = fits.open(str(LF_files0[i]))
        time[i] = Time(hdulist[0].header['DATE'], scale='utc').jd  # Scaling to utc time and Julian Date
        LF      = fits.getdata(str(LF_files1[i]))
        for j in range(N): # Loop over all stars and find flux (using same radius):
            flux_sky, n_pix_star, flux_star[i][j] = aperture(LF, X[j], Y[j], R, plot)
            SNR_i[i][j] = SNR(flux_sky, n_pix_star, flux_star[i][j], gain, ron)

    
    #--- SEEING CORRECTION ---#:
    mag_star = -2.5*log10(flux_star[:,0])+con          # From flux to magnitude for target star
    mag_ref  = -2.5*log10(flux_star[:,range(1,N)])+con # From flux to magnitude for reference stars
    # Find all delta-mag (dm):
    #print mag_star; sys.exit()
    dm  = zeros((n,N-1))                             
    for i in range(n):
        dm[i] = mag_ref[0,:]-mag_ref[i,:]              # dm from 0th frame to every other
    # Make trumpet plot:
    dm_med = median(dm,axis=0)
    # Using all stars:
    mag_all = mag_star + median(dm, axis=1)            # Using all stars

    # Which stars should be used to correct with:
    if dm_int==None:
        dm_int = [-0.01, 0.01, 3, 10]
    
    # Using interval:
    im_x = where((mag_ref[0,:]>dm_int[0])*(mag_ref[0,:]<dm_int[1]))[0]
    im_y = where((dm_med>dm_int[2])*(dm_med<dm_int[3]))[0]
    im_xy = list(set(im_x).intersection(im_y))
    mag_int = mag_star + median(dm[:,im_xy], axis=1)

    # Plot results:
    if plot==1 or plot==2:
        plot_trumpet(mag_ref, dm_med, dm_int, im_xy, '$m_i$', '$m_0\, -\, m_i$', 'Trumpet - Seeing Corretion')
        plt.show()
        MAG(time, [mag_star, mag_all, mag_int], ['.', 'r.', 'g*'], '$t$ (hours) - UTC JD', '$m$')
        plt.show()
        
    return


###############################################################################################

def star_finder(path, LF_name, sigma=2, plot=None, save=None):
    """
    This function finds all stellar object that is more lumious than TOL. 
    ----------INPUT:
    path           : Directory path to data
    N              : Number of frames
    LF_name        : Filename of Light Frame (LF)
    method         : Method to be used: Centroids='cen', PSF='psf'
    sigma          : Number of standard deviations used for flux criterion
    plot           : Plot==1 display found stars
    save           : save==1 saves a corrdinates [x, y]
    ---------OUTPUT:
    x_cen, y_cen   : x and y coordinates of found stars.
    """

    # Packages:
    import time, sys, glob
    import numpy as np
    from astropy.io import fits
    import scipy.ndimage as snd
    import matplotlib.pyplot as plt
    from Plot_Tools import FITS

    #-----------
    # Load data:
    #-----------
    
    LF_files = np.sort(glob.glob('{}{}*'.format(path, LF_name)))
    LF0      = fits.getdata(str(LF_files[0]))
    h, w     = np.shape(LF0) # Heigh and width
    n        = len(LF_files) # Number of images
    
    # Use same number of stars in all frames from statistic:
    threshold = mean(LF0) + sigma*std(LF0)
    labels, m = snd.label(LF0>threshold, ones((3,3)))
    print 'Using %i stars to shift images' %m
    
    #------------------
    # Align image data:
    #------------------

    LF_align = np.zeros((n, h, w))
    for i in range(n):
        LF_i    = fits.getdata(str(LF_files[i])) 
        cen     = snd.center_of_mass(LF_i, labels, range(m))
        x_cen_i = np.array(cen)[:,0]
        y_cen_i = np.array(cen)[:,1]

        if i==0: 
            # Initial mean values:
            x_ini = np.median(x_cen_i, axis=0)
            y_ini = np.median(y_cen_i, axis=0)
            
        # Frames are shifted using and interpolation of the frames:
        x_shift = np.median(x_cen_i, axis=0)
        y_shift = np.median(y_cen_i, axis=0)
        x_diff = np.diff([x_ini, x_shift])
        y_diff = np.diff([y_ini, y_shift])
        LF_align[i] = snd.interpolation.shift(LF_i, [x_diff, y_diff])

        # Save images if you like:
        if save==1:
            fits.writeto(('{}{}_align_%03d.fits'.format(path, LF_name) %i), LF_align[i], overwrite=True)

    #-------------
    # Coordinates:
    #-------------

    LF    = np.median(LF_align, axis=0)
    cen   = snd.center_of_mass(LF, labels, range(m))
    x_cen = np.array(cen)[:,0]
    y_cen = np.array(cen)[:,1]

    if plot==1: FITS(LF, 'linear', 2); plt.plot(y_cen, x_cen, 'r+'); plt.show()
    if save==1: np.savetxt('{}{}'.format(path, 'star_coor.txt'), np.vstack([x_cen, y_cen]).T)





def aperture(LF, x, y, R, plot):
    """
    This function find stellar flux with a circular aperture.
    ----------INPUT:
    LF             : Light Frame (LF)
    x, y           : Stellar coordinate
    R              : [r, r_min, r_max]: aperture, inner and outer radius for background
    plot           : Plot if you like. 
    ----------------"""
    # Packages:
    from numpy import meshgrid, ogrid, sqrt, mean, median, argmax, sum, cos, sin, nonzero
    # Split varibles:
    r, r_min, r_max = R[0], R[1], R[2]
    x, y = int(x), int(y)
    # Make meshgrid maks:
    x_grid, y_grid = meshgrid(range(x-r_max,x+r_max), range(y-r_max,y+r_max))
    square = LF[y_grid, x_grid]                    # Square frame of each image
    c      = sqrt((x_grid-x)**2 + (y_grid-y)**2)   # Equation of a circle
    # Background:
    c_sky = ((c>=r_min)*(c<=r_max))*square         # r_max^2 dimensional array
    sky_mean = mean(c_sky[nonzero(c_sky)])         # Mean sky flux
    sky_medi = median(c_sky[nonzero(c_sky)])       # Median sky flux
    flux_sky = 3*sky_medi - 2*sky_mean             # Sky flux
    # Star:
    c_star     = (c<=r)*square
    star       = c_star[nonzero(c_star)]-flux_sky  # Stellar corrected pixels
    n_pix_star = sum(c<=r)                         # Number of used star pixels
    flux_star  = sum(star)                         # Flux from star

    # Plot zoom-in on star with aperture: 
    if plot==3:
        FITS(square, 'linear')
        t = linspace(0,2*pi)
        x_off = (2*r_max+1)/2
        y_off = (2*r_max+1)/2
        plt.plot(x_off,y_off,'+r')
        plt.plot(r*cos(t)+x_off, r*sin(t)+y_off, 'g-')
        plt.plot(r_min*cos(t)+x_off, r_min*sin(t)+y_off, 'b-')
        plt.plot(r_max*cos(t)+x_off, r_max*sin(t)+y_off, 'b-')
        plt.show() # Shows only figure ones in loop
        sys.exit()
        
    return flux_sky, n_pix_star, flux_star
  


def SNR(flux_sky, n_pix_star, flux_star, gain, ron):
    """------------------------------------------- FUNCTION ------------------------------------------------:
    This function calculates the Signal-to-Noise Ratio (SNR). 
    ----------INPUT:
    flux_sky       : Flux from each sky pixel
    n_pix_star     : Number of pixels used to find star flux
    flux_star      : Flux from stellar object
    gain           : Gain CCD (e-/ADU)
    ron            : Read out noise (ron) (e-)
    ----------------------------------------------------------------------------- SNR """
    SNR = (gain*flux_star/sqrt(gain*flux_star + n_pix_star*gain*flux_sky + n_pix_star*ron**2))
    return SNR



def optimal_aperture(LF, x, y, R, gain, ron, plot):
    """------------------------------------------- FUNCTION ------------------------------------------------:
    This function 
    ----------INPUT:
    images         : 
    ---------OUTPUT:
    master_flat    : """
    #print '--------------------------------------------------------------------- optimal_aperture'
    # Constants:
    r_min, r_max = R[1], R[2]
    r_i = range(r_max)
    # Find SNR as a funcion of aperture radius:
    SNR_r = zeros(r_max)
    for r in r_i:
        x, y = int(x), int(y)
        # Make meshgrid maks:
        x_grid, y_grid = meshgrid(range(x-r_max-2,x+r_max+2), range(y-r_max-2,y+r_max+2))
        square = LF[y_grid, x_grid]                    # Square frame of each image
        c      = sqrt((x_grid-x)**2 + (y_grid-y)**2)   # Equation of a circle
        # Background:
        c_sky = ((c>=r_min)*(c<=r_max))*square         # r_max^2 dimensional array
        sky_mean = mean(c_sky[nonzero(c_sky)])         # Mean sky flux
        sky_medi = median(c_sky[nonzero(c_sky)])       # Median sky flux
        flux_sky = 3*sky_medi - 2*sky_mean             # Sky flux
        # Star:
        c_star     = (c<=r)*square
        star       = c_star[nonzero(c_star)]-flux_sky  # Stellar corrected pixels
        n_pix_star = sum(c<=r)                         # Number of used star pixels
        flux_star  = sum(star)                         # Flux from star
        # SNR:
        SNR_r[r] = SNR(flux_sky, n_pix_star, flux_star, gain, ron)
    # Plotting:
    plt.plot(range(r_max), SNR_r, '*')
    plt.xlabel('Aperture Radius')
    plt.ylabel('SNR')
    plt.show()   
    return SNR_r



def center_of_flux(LF, x, y, R):
    """This function finds the center of flux for all desired stellar object"""
    # Packages:
    from numpy import zeros, max, sum, round, meshgrid, argmax
    # Extract data:
    r, r_max = R[0], R[2]
    x, y     = x.astype(int), y.astype(int)
    # Outer: loops over all stars:
    center_x = zeros(len(x))
    center_y = zeros(len(y))
    for i in range(len(x)):
        x_grid, y_grid = meshgrid(range(x[i]-r_max,x[i]+r_max), range(y[i]-r_max,y[i]+r_max))
        square = LF[y_grid, x_grid]      # Square frame of each image
        # Test:
        FITS(square,'linear'); plt.plot(x_grid-x[i]+r_max, y_grid-y[i]+r_max, 'b*'); plt.show()
        # Inner: loops over all pixels: 
        for j in range(2*r):
            x_max_i = max(square)                       # Row and value of max pixel
            pixel_i = square.index(x_max_i)             # Index of max value pixel
            _, y_max_i       = max(square[x_max_i, :])   # Column for max pixel                   
            x_max[j]         = x_max_i                   # Rows gathers in a list
            y_max[j]         = y_max_i                   # Columns gathers in a list  
            pixel[j]         = pixel_i                   # Max pixels gathers in a list
            square[x_max_i, y_max_i] = 0                 # Max pixel is set to zero to find the nextr
        # Flux center is found:
        flux      = sum(pixel)                           # Total flux from pixels
        center_x  = 1/flux*dot(pixel,x_max)              # Flux-center in x
        center_y  = 1/flux*dot(pixel,y_max)              # Flux-center in y
    # Center of fluc for all stars:
    cen = [center_x-r, center_y-r]
    return cen





def plot_trumpet(mag_ref, dm_med, dm_int, im_xy, xlabel, ylabel, title):
    fig, ax = plt.subplots(1,1)
    plt.plot(mag_ref[0,:], dm_med, '.')
    ax.add_patch(patches.Rectangle((dm_int[0], dm_int[2]), dm_int[1]-dm_int[0], dm_int[3]-dm_int[2],\
                                   fill=False)) # Make empy box
    plt.plot(mag_ref[0,im_xy], dm_med[im_xy], 'g.')
    plot_settings(fig, ax, xlabel, ylabel, None, title)
    
