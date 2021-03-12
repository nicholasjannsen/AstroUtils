"""
CLASS DESCRIBTION:
------------------
Written by Nicholas Jannsen,  2018

This 
"""

# Loading packages:
from numpy import min, max
import math, sys, time
import numpy as np
# Functions:
import glob, pyfits, pylab, scipy
import matplotlib.pyplot as plt
import scipy.ndimage as snd
from scipy.misc import imsave
from astropy.io   import fits
from astropy.time import Time
# Own function:
from Plot_Tools import plot_settings, FITS, MAG, plot_trumpet

# For plotting:
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


###########################################################################################################
#                                            DEFINE CLASS                                                 #
###########################################################################################################

class Photometry(object):
    # INITILIZE THE CLASSE: 
    def __init__(self, path, LF_name, plot, save):

        # DEFINE GLOBAL VARIABLES (DGV)

        # Customized information:
        self.path    = path       # Directory path to data
        self.LF_name = LF_name    # Name of Light Frame (LF)
        self.plot    = plot       # Plot if 1
        self.save    = save       # Save if 1

        # Header information:
        self.LF_files = np.sort(glob.glob('{}{}*'.format(self.path, self.LF_name)))
        hdulist  = fits.open(str(self.LF_files[0]))
        header0 = hdulist[0].header
        try: # For now these has to be done manually
            #--- ORO ---#
            # self.target  = header0['OBJECT']
            # self.w       = header0['NAXIS1']   # Height of image 
            # self.h       = header0['NAXIS2']   # Width  of image
            # self.filt    = header0['FILTER']
            # self.exptime = header0['EXPTIME']
            # self.time    = Time(header0['TIME'], scale='utc').jd # Scaling to utc time and Julian Date 2000
            #--- IAC80 ---#
            self.target  = header0['OBJECT']
            self.w       = header0['NAXIS1']   # Height of image 
            self.h       = header0['NAXIS2']   # Width  of image
            self.filt    = header0['INSFILTE']
            self.exptime = header0['EXPTIME']
            self.time    = Time(header0['DATE'], scale='utc').jd # Scaling to utc time and Julian Date 2000

        except NameError: print('ERROR: FITS HEADER INFORMATION DO NOT MATCH PROGRAM VALUES'); return

###########################################################################################################
#                                            MAIN FUNCTIONS                                               #
###########################################################################################################

    def aperture_photometry(self, coor_name, coor_target, R, dm_int=None, plot=None, save=None):
        """------------------------------------------- MAIN FUNCTION --------------------------------------
        This function 
        ----------INPUT:
        images         : 
        ---------OUTPUT:
        master_flat    : """
        print '--------------------------------------------------------------------- aperture_photometry'
        #start_time = time.time()
        from Photometry import aperture, SNR, center_of_flux

        #--- CONSTANTS ---#
        gain  = 0.73          # Gain of camera: electrons pr ADU (ADU = counts from object)- ajust to camera!
        ron   = 3.3           # Read out noise - ajust to camera!
        con   = 25            # Magnitude constant
        r     = R[0]          # Aperture radius = R[1]          # Minimum aperture radius to measure sky flux
        rmin  = R[1]
        rmax  = R[2]          # Maximum aperture radius to measure sky flux

        #--- STELLAR COORDINATES ---#

        # Load stellar coordinates:
        if type(coor_name)==str:  coor_stars = np.loadtxt('{}{}.txt'.format(self.path, coor_name))
        if type(coor_name)==list: coor_stars = coor_stars 
        else: print('ERROR: STELLAR COORDINATES HAS A WRONG FORMAT') 
        x, y   = coor_stars[:,1], coor_stars[:,0] # Stellar coordinates
        xt, yt = coor_target[0],  coor_target[1]  # Target  coordinates
        N = len(x)                                # Number of stars used
        
        # Make target star the first star:
        d = np.zeros(N)
        for i in range(N): d[i] = np.sqrt((xt-x[i])**2+(yt-y[i])**2)  # Simple pythagorean geometry
        index_min = np.argmin(d)                                      # Find target star by geometry
        x = np.roll(x, -index_min)                                    # Shift data to target star index
        y = np.roll(y, -index_min)
        # Check if stellar coordinates are closer than r_max:
        i_x1 = np.where(x<rmax)[0]; i_x2 = np.where(x>self.w-rmax)[0]
        i_y1 = np.where(y<rmax)[0]; i_y2 = np.where(y>self.h-rmax)[0]
        i_xy = np.hstack([i_x1, i_x2, i_y1, i_y2])
        # Discard these coordinates:
        X = np.delete(x, i_xy)
        Y = np.delete(y, i_xy)
        N = len(X) # New number of stars
        # Check coordinates if you like:
        if plot==2:
            LF = pyfits.getdata(str(self.LF_files[0]))
            fig, ax = plt.subplots(1,1)
            FITS(LF,'linear')
            plt.plot(x,y,'ro',mfc='none'); plt.plot(x[0],y[0],'b+'); plt.plot(x[i_xy], y[i_xy],'y*')
            plot_settings(fig, ax)
            plt.show()

        #--- FIND OPTIMAL APERTURE RADIUS ---#
        if plot==3:
            LF = pyfits.getdata(str(self.LF_files[0]))
            SNR_r = optimal_aperture(LF, X[0], Y[0], R, gain, ron, plot)

        #--- PHOTOMETRY ---#
        # Find fluxes:
        flux_star = np.zeros((n,N))
        SNR_i     = np.zeros((n,N))
        time      = np.zeros(n) 
        for i in range(n): # Loop over all images:
            hdulist = fits.open('{}{}_%03d.fits'.format(path, LF_name) %i)
            LF      = fits.getdata('{}{}_%03d.fits'.format(path, LF_name) %i)
            time[i] = Time(hdulist[0].header['date'], scale='utc').jd  # Scaling to utc time and Julian Date
            for j in range(N): # Loop over all stars and find flux (using same radius):
                flux_sky, n_pix_star, flux_star[i][j] = aperture(LF, X[j], Y[j], R, plot)
                SNR_i[i][j] = SNR(flux_sky, n_pix_star, flux_star[i][j], gain, ron)

        #--- SEEING CORRECTION ---#:
        mag_star = -2.5*log10(flux_star[:,0])+con          # From flux to magnitude for target star
        mag_ref  = -2.5*log10(flux_star[:,range(1,N)])+con # From flux to magnitude for reference stars
        # Find all delta-mag (dm):
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
        # Plot trumpet:
        if plot==1 or plot==2:
            plot_trumpet(mag_ref, dm_med, dm_int, im_xy, '$m_i$', '$m_0\, -\, m_i$', \
                         'Trumpet - Seeing Corretion')
            plt.show()
            MAG(time, [mag_star, mag_all, mag_int], ['.', 'r.', 'g*'], '$t$ (days) - UTC JD', '$m$'); plt.show()

        return

###########################################################################################################
#                                            SUB FUNCTIONS                                                #
###########################################################################################################

    def aperture(LF, x, y, R, plot):
        """------------------------------------------- FUNCTION ------------------------------------------------:
        This function find stellar flux with a circular aperture.
        ----------INPUT:
        LF             : Light Frame (LF)
        x, y           : Stellar coordinate
        R              : [r, r_min, r_max]: aperture, inner and outer radius for background
        plot           : Plot if you like. 
        ----------------------------------------------------------------------------- aperture """
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


    def aperture(self, LF, x, y, aperture, background):
        """
        This function calculate the stellar and sky background flux either using a highly elliptical aperture
        or a aperture that traces the the Center Of Flux (COF) to fit the startrails. The routine returns 
        stellar flux 'flux_star', the sky background flux 'flux_sky', and the number of stellar pixels
        'n_pix_star'.
        """
        # Aperture only handle integers:
        x, y = int(x), int(y)

        # Grid for aperture:
        x_grid, y_grid = meshgrid(range(self.w), range(self.h))
        grid = LF[y_grid, x_grid]  # Grid
        xx = x_grid - x            # Displacement for stellar x coor
        yy = y_grid - y            # Displacement for stellar y coor

        # The aperture:
        a     = aperture[1]                # Star radius/semi-minor axis
        b     = aperture[2]                # Star semi-major axis (circle if a = b)
        r_min = aperture[3]                # min radius
        r_max = aperture[4]                # 
        phi   = math.radians(aperture[5])  # Tilt angle: [0:180] deg

        # Global background:
        if background=='global':
            flux_sky = self.global_sky_background(LF)    

        #--- CIRCULAR OR ELLIPTIC APERTURE ---#
        if aperture[0]=='ellipse':

            # Parametrasation and rotation of the ellipse (1 = x^2/a^2 + y^2/b^2)
            phi = phi - pi/2 # The ellipse have an offset from unity circle 
            EE_star    = (xx*cos(phi)+yy*sin(phi))**2/a**2         + (xx*sin(phi)-yy*cos(phi))**2/b**2
            EE_sky_min = (xx*cos(phi)+yy*sin(phi))**2/(a+r_min)**2 + (xx*sin(phi)-yy*cos(phi))**2/(b+r_min)**2
            EE_sky_max = (xx*cos(phi)+yy*sin(phi))**2/(a+r_max)**2 + (xx*sin(phi)-yy*cos(phi))**2/(b+r_max)**2

            # Local background:
            if background=='local':
                sky_grid = ((EE_sky_min>1)*(EE_sky_max<1))*grid   # Sky background determined by width q 
                sky_pixs = sky_grid[nonzero(sky_grid)]            # Sky pixel values
                flux_sky = 3*median(sky_pixs) - 2*mean(sky_pixs)  # Robust sky background flux

            # Star: 
            star_grid  = (EE_star<=1)*grid
            star_pixs  = star_grid[nonzero(star_grid)]-flux_sky  # Stellar corrected pixels
            n_pix_star = sum(EE_star<=1)                         # Number of used star pixels
            flux_star  = sum(star_pixs)                          # Flux from star

        #--- TRACE APERTURE ---#
        if aperture[0]=='trace':

            # Star initial:
            CC_star    = sqrt(xx**2 + yy**2) - a        # Start frame used to trace from:
            star_grid  = (CC_star<=1)*grid              # Star images
            star_pixs  = star_grid[nonzero(star_grid)]  # Stellar pixels
            n_pix_star = len(star_pixs)                 # Number of pixels in circle

            # Loop-step in x or y depends on phi:
            # Step in x and finds y centroid for each step:
            if 0<=phi<=pi/4 or pi*3/4<=phi<=pi:
                step = 'x step'
                x_step = range(b)
                y_step = zeros(b)
            # Step in y and finds x centroid for each step:
            if pi/4<phi<pi*3/4:
                step = 'y step'
                x_step = zeros(b)
                y_step = range(b)

            # Loop over trace step:
            x_cen = zeros(b)
            y_cen = zeros(b)
            star_true = zeros((b, self.h, self.w))
            sky_true  = zeros((b, self.h, self.w))
            for i in range(b):
                # Find Center Of Flux (COF):
                y_cen[i], x_cen[i] = self.center_of_flux(star_img, n_pix_star)
                if step=='x step': XX = x + x_step[i]; YY = y_cen[i] + y_step[i]
                if step=='y step': YY = y + y_step[i]; XX = x_cen[i] + x_step[i]
                # A circular aperture is used to trace with:
                CC_star = sqrt((x_grid-XX)**2 + (y_grid-YY)**2) - a
                star_true[i] = (CC_star<=1)                      # Array of true and false statement: star
                star_img     = star_true[i]*grid                 # Image: star
                n_pix_star   = len(star_img[nonzero(star_img)])  # Number of pixels: star
                # If local background:
                if background=='local':
                    CC_sky  = sqrt((x_grid-XX)**2 + (y_grid-YY)**2) - a_sky
                    sky_true[i] = ((CC_star>1)*(CC_sky<1))*grid  # Array of true and false statement: sky

            # Local sky flux:
            if background=='local':
                star_false = np.logical_not(np.sum(star_true, axis=0))  # Invert to get rid of star
                sky_true_x = np.sum(sky_true, axis=0) > 0
                sky_img  = sky_true_x.astype(np.int)*star_false*LF      # Here bol: True*False = False
                sky      = sky_img[nonzero(sky_img)]
                flux_sky = 3*median(sky) - 2*mean(sky)                  # Robust sky flux

            # Stellar flux: 
            star_true_x = np.sum(star_true, axis=0) > 0         # Bool array of stellar pixels
            star_img    = star_true_x.astype(np.int)*LF         # Stellar image
            star        = star_img[nonzero(star_img)]-flux_sky  # Stellar corrected pixels
            n_star_pix  = len(star)
            flux_star   = sum(star_img)

        # PLOT IF YOU LIKE:
        if self.plot==1:
            # Elliptic aperture with local background:
            if aperture[0]=='ellipse':
                from plot_tools import plot_ellipse
                FITS(LF, 'linear', 2)
                plot_ellipse(a,     b,     math.degrees(phi), x, y, 'g')  # Stellar aperture
                plot_ellipse(a_sky, b_sky, math.degrees(phi), x, y, 'm')  # Background aperture
                plt.show()
            # Box aperture with local background:
            if aperture[0]=='trace':
                 t = linspace(0, 2*pi)
                 FITS(LF, 'linear', 2)
                 [plt.plot(a_sky*cos(t)+(x+x_step[i]), a_sky*sin(t)+y_cen[i], 'b-') for i in range(b)]
                 [plt.plot(a*cos(t)+(x+x_step[i]), a*sin(t)+y_cen[i], 'g-') for i in range(b)]
                 plt.show()

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
        print '--------------------------------------------------------------------- optimal_aperture'
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

   
# def center_of_flux(self, LF, n_pix):
#         """
#         This function finds the center of flux for all desired stellar object. Here LF is the masked image
#         thus every pixel is set to zero except the star) and n_pix is the number of pixels one wish to use
#         in order to find the COF.
#         """
#         # Loops over all pixels:
#         LF_copy  = copy(LF)     # Copy to avoid overwriting
#         flux_max = zeros(n_pix)
#         x_max = zeros(n_pix)
#         y_max = zeros(n_pix)
#         pixel = zeros(n_pix)
#         for j in range(n_pix):
#             flux_max[j] = np.max(LF_copy)               # Maximum value for array
#             max_dex = np.where(LF_copy == flux_max[j])  # Find row, column for min value
#             x_max[j] = max_dex[0][0]                    # max for x coordinate
#             y_max[j] = max_dex[1][0]                    # max for y coordinate
#             pixel[j] = j
#             # Min pixel is et to max in order to find the next min:
#             LF_copy[int(y_max[j]), int(x_max[j])] = 0 

#         # Flux center is found:
#         flux = sum(pixel)
#         cen_x = 1/flux*dot(pixel,x_max)              # Flux-center in x
#         cen_y = 1/flux*dot(pixel,y_max)              # Flux-center in y
#         return cen_x, cen_y


 # def global_sky_background(self, LF):
 #        """
 #        This function calculates the global sky background flux. 200 minimum pixels are used as background
 #        as a standard. As the background flux may not be uniform across CCD in space, the light frame is
 #        divided into 's' number of subframes, thus, in each subframe a total number of 'n' number of pixels
 #        are found across the sky and used in the sky background flux.
 #        ----------INPUT:
 #        LF             : A single Light Frame (LF)
 #        ---------OUTPUT:
 #        flux_sky       : Global sky background flux
 #        """
 #        # Variables:
 #        s     = 9                             # Number of subframes (CHANGE IF NEEDED!) E.g. 4, 9, 16 etc. 
 #        n     = self.h*self.w/(self.h+self.w) # Number of pixels used in subframes scales with image dim 
 #        nrows = self.h/(s/2)                  # Numbers of rows in each subframe
 #        ncols = self.w/(s/2)                  # Numbers of columns in each subframe

 #        # Reshape light frame into subframe:
 #        LF_sub = (LF.reshape(self.h//nrows, nrows, -1, ncols).swapaxes(1,2).reshape(-1, nrows, ncols))

 #        # Loop over all subframes:
 #        min_val = np.zeros((s,n))
 #        for i in range(s):
 #            # Loops over all pixels:
 #            for j in range(n):
 #                min_val[i,j] = np.min(LF_sub[i])                # Minimum value for array
 #                min_dex = np.where(LF_sub[i] == min_val[i,j])   # Find row, column for min value
 #                # Min pixel is set to max in order to find the next min:
 #                LF_sub[i, min_dex[0][0], min_dex[1][0]] = np.max(LF_sub[i]) 

 #        # Flux:
 #        flux_sky = 3*median(min_val) - 2*mean(min_val)  # Mean flux from pixels
 #        return flux_sky
