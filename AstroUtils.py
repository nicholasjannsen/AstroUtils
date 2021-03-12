# Packages:
import math, sys, time, scipy, glob
from PIL import Image
import numpy as np
import pylab
import heapq

# Astropy:
from astropy.io import fits
from astropy.time import Time

# SciPy:
import scipy.misc as smi
import scipy.ndimage as snd

# Matplotlib:
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# Own functions:
import Plot_Tools as pt
np.set_printoptions(suppress=True, formatter={'float_kind':'{:7.3f}'.format}, linewidth=100)

#########################################################################################
#                                 IMAGE REDUCTION                                       #
#########################################################################################

def image_reduction(path, date, telescope='not', plot=0, save=0):
    """
    This function takes the directory path as input and load all
    image files. Based on the header reduction is automatically
    handled by the software if only clear filter images is present
    or other filters are available as well. The software returns
    saves the fully calibrated light frames with the originally
    filename extension added with a 'C' for calibrated in front.
    A master bias, dark current, flat(s) are saved as well with
    the names 'master_BF.fits', 'master_DF.fits', master_FF(-'filter').fits.
    If all calibration images already exist the routine will simply
    pass to the plots.
    ----------------------------
                INPUT          :
    ----------------------------
    path               (string): Directory path to data
    telescope          (string): Name of telescope
    plot                  (int): If plot=1: plots all relevant frames
    save                  (int): If save=1: saves all corrected light frames in seperate files
    ----------------------------
               OUTPUT          :
    ----------------------------
    Master_BF            (fits): Master bias frame
    Master_DF            (fits): Master dark frame
    Master_FF(-'filter') (fits): Master flat(s) frame
    C*.fits              (fits): Calibrated light frame
    """
    #------------------------------------------
    # TEST IF CALIBRATION IMAGES ALREADY EXIST:
    #------------------------------------------
    try:
        #print('{}../calib/master'.format(path))
        #sys.exit()
        BF = fits.open('{}../calib/master_'.format(path))

    #------------------------------------------
    # LOAD DATA AND EXTRACT HEADER INFORMATION:
    #------------------------------------------
    except IOError:
        # Load all files from object and calib folder:
        print('Loading image files... \n')
        files_object = np.sort(glob.glob('{}*'.format(path)))
        files_calibs = np.sort(glob.glob('{}../calibs/{}/*'.format(path, date)))
        hdu_object = np.array([fits.open(str(files)) for files in files_object])
        hdu_calibs = np.array([fits.open(str(files)) for files in files_calibs])

        # Separate files from each other based on telescope and header:
        type_object = [hdu_object[i][0].header['IMAGETYP'] for i in range(len(files_object))]
        type_calibs = [hdu_calibs[i][0].header['IMAGETYP'] for i in range(len(files_calibs))]
        if telescope is 'oro28': img_type = ['bias', 'dark', 'flat', 'object']
        if telescope is 'oro50': img_type = ['bias', 'dark', 'flat', 'object']
        if telescope is 'not':   img_type = ['BIAS',  None , 'FLAT,SKY', 'OBJECT']
        BF_dex = np.where(np.array(type_calibs)==img_type[0])[0]
        DF_dex = np.where(np.array(type_calibs)==img_type[1])[0]
        FF_dex = np.where(np.array(type_calibs)==img_type[2])[0]
        LF_dex = np.where(np.array(type_object)==img_type[3])[0]

        # Load images into arrays:
        BF_i = np.array([fits.getdata(str(files_calibs[i])) for i in BF_dex])
        DF_i = np.array([fits.getdata(str(files_calibs[i])) for i in DF_dex])
        FF_i = np.array([fits.getdata(str(files_calibs[i])) for i in FF_dex])

        print(BF_i)
        sys.exit()

        # Image constants:
        n    = len(files_object)
        h, w = np.shape(hdu_object[0][0].data)

        # Find exposure times:
        DF_exp = [hdu_calibs[DF_dex[i]][0].header['EXPTIME'] for i in range(len(DF_dex))]
        FF_exp = [hdu_calibs[FF_dex[i]][0].header['EXPTIME'] for i in range(len(FF_dex))]
        LF_exp = [hdu_object[LF_dex[i]][0].header['EXPTIME'] for i in range(len(LF_dex))]

        # Test if different filters exist for the light frames:
        if telescope is 'not':
            # ALFOSC filters:
            FF_filt_type_alfosc = np.array([hdu_calibs[i][0].header['ALFLTNM'] for i in FF_dex], dtype='<U6')
            LF_filt_type_alfosc = np.array([hdu_object[i][0].header['ALFLTNM'] for i in LF_dex], dtype='<U6')
            # FASU A filters:
            FF_filt_type_fasua  = np.array([hdu_calibs[i][0].header['FAFLTNM'] for i in FF_dex], dtype='<U6')
            LF_filt_type_fasua  = np.array([hdu_object[i][0].header['FAFLTNM'] for i in LF_dex], dtype='<U6')
            # FASU B filters:
            FF_filt_type_fasub  = np.array([hdu_calibs[i][0].header['FBFLTNM'] for i in FF_dex], dtype='<U6')
            LF_filt_type_fasub  = np.array([hdu_object[i][0].header['FBFLTNM'] for i in LF_dex], dtype='<U6')
            # Correct to right filter (OBS the order here is very specific!):
            FF_filt = np.char.add(FF_filt_type_alfosc, FF_filt_type_fasua)
            FF_filt = np.char.add(FF_filt, FF_filt_type_fasub)
            FF_filt = np.array([str(FF_filt[i]).replace('OpenOpenOpen','Clear') for i in range(len(FF_filt))])
            FF_filt = np.array([str(FF_filt[i]).replace('Open', '')       for i in range(len(FF_filt))])
            FF_filt = np.array([str(FF_filt[i]).replace('Halp 6', 'Halp') for i in range(len(FF_filt))])
            FF_filt = np.array([str(FF_filt[i]).replace('Hbet 4', 'Hbet') for i in range(len(FF_filt))])
            FF_filt = np.array([str(FF_filt[i]).replace("'_SDS", '_SDSS') for i in range(len(FF_filt))])
            FF_filt = np.array([str(FF_filt[i]).replace('int', 'SDSS')    for i in range(len(FF_filt))])
            FF_filt = np.array([str(FF_filt[i]).replace('[', '')          for i in range(len(FF_filt))])
            FF_filt = np.array([str(FF_filt[i]).replace(']', '')          for i in range(len(FF_filt))])
            FF_filt = np.array([str(FF_filt[i]).replace(' ', '')          for i in range(len(FF_filt))])
            #------
            LF_filt = np.char.add(LF_filt_type_alfosc, LF_filt_type_fasua)
            LF_filt = np.char.add(LF_filt, LF_filt_type_fasub)
            LF_filt = np.array([str(LF_filt[i]).replace('OpenOpenOpen','Clear') for i in range(len(LF_filt))])
            LF_filt = np.array([str(LF_filt[i]).replace('Open', '')       for i in range(len(LF_filt))])
            LF_filt = np.array([str(LF_filt[i]).replace('Halp 6', 'Halp') for i in range(len(LF_filt))])
            LF_filt = np.array([str(LF_filt[i]).replace('Hbet 4', 'Hbet') for i in range(len(LF_filt))])
            LF_filt = np.array([str(LF_filt[i]).replace("'_SDS", '_SDSS') for i in range(len(LF_filt))])
            LF_filt = np.array([str(LF_filt[i]).replace('int', 'Bes')    for i in range(len(LF_filt))])
            LF_filt = np.array([str(LF_filt[i]).replace('[', '')          for i in range(len(LF_filt))])
            LF_filt = np.array([str(LF_filt[i]).replace(']', '')          for i in range(len(LF_filt))])
            LF_filt = np.array([str(LF_filt[i]).replace(' ', '')          for i in range(len(LF_filt))])

        # Same for other telescopes:
        else:
            FF_filt = [hdu_calibs[i][0].header['INSFILTE'] for i in FF_dex]
            LF_filt = [hdu_object[i][0].header['INSFILTE'] for i in LF_dex]
        # Find index of the same flats:
        FF_dif      = np.unique(FF_filt)
        FF_filt_dif = np.array([np.where(FF_filt==FF_dif[i])[0]  for i in range(len(FF_dif))])
        # Find index of flats filter corresponding to the light frames:
        FF_filt_dex = np.array([np.where(FF_dif==LF_filt[i])[0] for i in range(len(LF_dex))])

        # If the same filter appears several times: TODO!

        #---------------------------------
        # MAKE MASTER CALIBRATION FRAMES :
        #---------------------------------

        # Bias:
        if BF_dex.any():
            BF = np.median(BF_i, axis=0)

        # Dark:
        if DF_dex.any() and BF_dex.any() and FF_name.any():
            DF = (DF_exptimes/LF_exptimes) * np.median(DF_i, axis=0)

        # Bias and Dark:
        if BF_dex.any() and DF_dex.any() and not FF_name.any():
            DF = (DF_exptimes/LF_exptimes) * np.median(DF_i - BF, axis=0)

        # Flats (with/without Bias and Darks:
        if FF_dex.any():
            FF_k = np.zeros((len(FF_filt_dif), h, w))
            # Loop over each filter
            for i in range(len(FF_filt_dif)):
                FF_j = np.array([FF_i[j] for j in FF_filt_dif[i]])
                # Different cases:
                if not BF_dex.any() and not DF_dex.any(): FF = np.median(FF_j          , axis=0)
                if     BF_dex.any() and not DF_dex.any(): FF = np.median(FF_j - BF     , axis=0)
                if not BF_dex.any() and     DF_dex.any(): FF = np.median(FF_j - DF     , axis=0)
                if     BF_dex.any() and     DF_dex.any(): FF = np.median(FF_j - BF - DF, axis=0)
                # If negative pixel values:
                FF[FF<=0] = np.median(FF)
                # Normalization and store:
                FF = FF/np.median(FF)
                FF_k[i] = FF

    #---------------------
    # PERFORM CALIBRATION:
    #---------------------

    # Find index of flats filter corresponding to the light frames:
    CF_i = np.zeros((len(files_object), h, w))
    for i in range(n):
        # Open and close and singel image at a time:
        with fits.open(str(files_object[i])) as hdu_i:
            LF_i = hdu_i[1].data
            # Cut out image if ALFOSC images:
            if telescope is 'not': LF_i = LF_i[ borders[0]:borders[1], borders[2]:borders[3]]
            # Check for calib images:
            if not BF_dex.any(): BF   = np.zeros(np.shape(LF_i))
            if not DF_dex.any(): DF   = np.zeros(np.shape(LF_i))
            if not FF_dex.any(): FF_k = np.one(len(files_object), h, w)
            if FF_dex.any(): j = FF_filt_dex[i]
            else: j = i
            # Make calibration:
            #print(i)
            #print(j)
            CF = (LF_i - BF)/FF_k[j]
            CF_i[i] = CF

        # Print compilation time to bash:
        pt.compilation(i, n, 'Image Reduction')
    print

    #----------------------------------------------------
    # Cutout borders if needed (cutout=[x1, x2, y1, y2]):
    #----------------------------------------------------
    if cutout==None: cutout = [0, len(BF[:,0]), 0, len(BF[:,1])]
    BF = BF[:, cutout[2]:cutout[3], cutout[0]:cutout[1]]
    DF = DF[:, cutout[2]:cutout[3], cutout[0]:cutout[1]]
    FF_k = [FF_k[i, cutout[2]:cutout[3], cutout[0]:cutout[1]] for i in range(len(FF_i))]
    SF_i = [SF_i[i, cutout[2]:cutout[3], cutout[0]:cutout[1]] for i in range(len(SF_i))]

    #-------------
    # SAVE IMAGES:
    #-------------
    if save==1:
        # Save master bias:
        if BF_dex.any():
            fits.writeto('{}../calibs/master_BF.fits'.format(path), \
                         BF, header=hdu_calibs[BF_dex[0]][0].header, overwrite=True)
        # Save master dark:
        if DF_dex.any():
            fits.writeto('{}../calibs/master_DF.fits'.format(path), \
                         DF, header=hdu_calibs[DF_dex[0]][0].header, overwrite=True)
        # Save master flats:
        if FF_dex.any():
            for i in range(len(FF_k)):
                fits.writeto('{}../calibs/master_FF_{}.fits'.format(path, FF_dif[i]), \
                             FF_k[i], header=hdu_calibs[FF_filt_dif[i][0]][0].header, overwrite=True)
        # Save science images:
        for i in range(len(CF_i)):
            fits.writeto('{}calib_{}.fits'.format(path, LF_filt[i]), \
                         CF_i[i], header=hdu_object[i][0].header, overwrite=True)

    #-------------
    # PLOT IMAGES:
    #-------------
    if plot==1:
        # if BF_dex.any():
        #     pt.FITS(BF, 'linear', 5); plt.show()
        #     print('Bias -- mean: {:2f}, std: {:2f}'.format(np.mean(BF), np.std(BF)))
        # if DF_dex.any():
        #     pt.FITS(DF, 'linear'); plt.show()
        #     print('Dark -- mean: {:2f}, std: {:2f}'.format(np.mean(DF), np.std(DF)))
        # if FF_dex.any(): pt.FITS(FF, 'linear', 5); plt.show()
        # pt.FITS(LF_i,    'linear', 2); plt.show()
        for i in range(n):
            plt.figure()
            pt.FITS(CF_i[i], 'linear', 2)
            plt.show()

    return

#########################################################################################
#                          REMOVE DEAD/HOT/COSMIC PIXELS                                #
#########################################################################################

def pixel_outliers(path, LF_name, sigma=2, plot=1, save=1):
    """
    This function takes all loaded flat-frames, dark-frames, bias-frames and combine them to a one
    master-flat-field image.
    ----------INPUT:
    path           : Directory path to data.
    LF_name        : Name of Light Frames (LF) except end number.
    FF_name        : Name of Flat  Frames (FF) except end number.
    DF_name        : Name of Dark  Frames (DF) except end number.
    BF_name        : Name of Bias  Frames (BF) except end number.
    N              : Number of [LF, FF, DF, BF] frames.
    plot           : If plot=1: plots all relevant frames..
    save           : If save=1: saves all corrected light frames in seperate files.
    ---------OUTPUT:
    CF_i           : Cube of Corrected Frames (CF).
    """
    # Load data:
    LF_files = np.sort(glob.glob('{}{}*'.format(path, LF_name)))
    hdu_i    = np.array([fits.open(str(files)) for files in LF_files])
    # Constants:
    n = len(LF_files)
    if LF_name=='calib': h, w = np.shape(hdu_i[0][0].data)
    else: h, w = np.shape(hdu_i[0].data)

    # ALFOSC filters:
    filt_type_alfosc = np.array([hdu_i[i][0].header['ALFLTNM'] for i in range(n)], dtype='<U5')
    filt_type_fasua  = np.array([hdu_i[i][0].header['FAFLTNM'] for i in range(n)], dtype='<U5')
    filt_type_fasub  = np.array([hdu_i[i][0].header['FBFLTNM'] for i in range(n)], dtype='<U5')
    # Correct to right filter:
    filt = np.char.add(filt_type_alfosc, filt_type_fasua)
    filt = np.char.add(filt, filt_type_fasub)
    filt = np.array([str(filt[i]).replace('Open', '') for i in range(len(filt))])
    filt = np.array([str(filt[i]).replace(' ', 'a')   for i in range(len(filt))])

    #------------------------------
    # REMOVE HOT/DEAD/COMIC PIXELS:
    #------------------------------

    COF_hotpix = np.array([])
    for i in range(n):
        with fits.open(str(LF_files[i])) as hdu_i:
            LF_i = hdu_i[0].data

            # HOT PIXELS:
            # Locate pixels structures above threshold:
            threshold    = np.mean(LF_i) + sigma*np.std(LF_i)
            thresmap_hot = np.where(LF_i > threshold, 1, 0)
            labels, N_struc = snd.label(thresmap_hot, structure=np.ones((3,3)))
            sums = snd.sum(thresmap_hot, labels, range(1, N_struc+1)).astype(int)
            # Find hot pixels as small structures:
            min_pix = 1
            hotpix = np.where(sums <= min_pix)[0] + 1
            # Find hot pixel coordinates:
            hotmap = np.zeros(np.shape(LF_i))
            for pix in hotpix: hotmap = hotmap + np.where(labels==pix, 1, 0)
            labels, N_struc = snd.label(hotmap, structure=np.ones((3,3)))
            coor_hot = np.asarray(snd.center_of_mass(LF_i, labels, range(1, N_struc+1))).astype(int)
            # Replace hot pixels with mean-box-value:
            for i in range(len(coor_hot)):
                b1, b2 = 1, 2
                x, y = coor_hot[i,0], coor_hot[i,1]
                xg1, yg1 = np.meshgrid(range(x-int(b1),x+int(b1)), range(y-int(b1),y+int(b1)))
                xg2, yg2 = np.meshgrid(range(x-int(b2),x+int(b2)), range(y-int(b2),y+int(b2)))
                flux_box_diff = np.mean(LF_i[xg1, yg1])/np.mean(LF_i[xg2, yg2])
                dex = np.where(flux_box_diff<1)
                print(dex)

                #print(flux_box/((box+1)**2))
            print(np.median(LF_i))
            sys.exit()

                #LF_i[coor_hot[i,0], coor_hot[i,1]] = np.median(LF_i[x_box, y_box])
            sys.exit()

            # #r_hot   = np.sqrt(sums[hotpix-1]/np.pi)
            plt.figure()
            pt.FITS(LF_i, 'linear', 2)
            plt.scatter(coor_hot[:,1], coor_hot[:,0], marker='s', facecolors='none', edgecolors='r')
            plt.title('{} stars are idenfied'.format(N_struc))
            plt.xlabel('X (pixels)'); plt.ylabel('Y (pixels)')
            plt.show()
            sys.exit()

            # Save each image:
            if save==1:
                fits.writeto('{}hotpix_{}.fits'.format(path, filt[i]), \
                             CF_i[i], header=hdu_i[i][0].header, overwrite=True)

        # Print compilation time to bash:
        pt.compilation(i, n, 'Hot/Dead pixel removal')
    print #----------------------------------------------------------------------
    print('\n {} Hot/dead pixel found'.format(len((hot_pixels))))

        # # Find good threshold to locate pixels:
        # blur_i = snd.median_filter(CF_i[i], size=2)
        # diff_i = CF_i[i] - blur_i
        # threshold =  np.mean(diff_i) + sigma*np.std(diff_i)
        # # Find the hot pixels (ignoring the edges):
        # hot_pixels = np.nonzero((abs(diff_i[1:-1, 1:-1])>threshold))
        # hot_pixels = np.array(hot_pixels) + 1  # +1 because 1. row and column was ignored
        # for y,x in zip(hot_pixels[0], hot_pixels[1]):
        #     CF_i[i,y,x] = blur_i[y,x]

        # # PIXEL AT EDGES (BUT NOT CORNERS):
        # # Left and right sides:
        # for index in range(1, h-1):
        #     # Left side:
        #     med  = np.median(CF_i[i, index-1:index+2, 0:2])
        #     diff = abs(CF_i[i, index, 0] - med)
        #     if diff>threshold: 
        #         hot_pixels = np.hstack((hot_pixels, [[index],[0]]))
        #         CF_i[i, index, 0] = med
        #         # Right side:
        #     med  = np.median(CF_i[i, index-1:index+2, -2:])
        #     diff = abs(CF_i[i, index,-1] - med)
        #     if diff>threshold: 
        #         hot_pixels = np.hstack((hot_pixels, [[index], [w-1]] ))
        #         CF_i[i, index, -1] = med

        # # TOP AND BOTTOM:
        # for index in range(1, w-1):
        #     # Bottom:
        #     med  = np.median(CF_i[i, 0:2, index-1:index+2])
        #     diff = abs(CF_i[i, 0, index] - med)
        #     if diff>threshold: 
        #         hot_pixels = np.hstack((hot_pixels, [[0], [index]] ))
        #         CF_i[i, 0, index] = med
        #         # Top:
        #     med  = np.median(CF_i[i, -2:, index-1:index+2])
        #     diff = abs(CF_i[i, -1, index] - med)
        #     if diff>threshold: 
        #         hot_pixels = np.hstack((hot_pixels, [[h-1], [index]]))
        #         CF_i[i, -1, index] = med

        # # LASTLY CORNERS:
        # # Bottom left:
        # med  = np.median(CF_i[i, 0:2, 0:2])
        # diff = abs(CF_i[i, 0, 0] - med)
        # if diff>threshold: 
        #     hot_pixels = np.hstack((hot_pixels, [[0], [0]] ))
        #     CF_i[i, 0, 0] = med
        #     # Bottom right:
        # med  = np.median(CF_i[i, 0:2, -2:])
        # diff = abs(CF_i[i, 0, -1] - med)
        # if diff>threshold: 
        #     hot_pixels = np.hstack((hot_pixels, [[0], [w-1]] ))
        #     CF_i[i, 0, -1] = med
        #     # Top left:
        # med  = np.median(CF_i[i, -2:, 0:2])
        # diff = abs(CF_i[i, -1, 0] - med)
        # if diff>threshold: 
        #     hot_pixels = np.hstack((hot_pixels, [[h-1], [0]] ))
        #     CF_i[i, -1, 0] = med
        #     # Top right:
        # med  = np.median(CF_i[i, -2:, -2:])
        # diff = abs(CF_i[i, -1, -1] - med)
        # if diff>threshold: 
        #     hot_pixels = np.hstack((hot_pixels, [[h-1], [w-1]] ))
        #     CF_i[i, -1, -1] = med

    #     # Save each image:
    #     if save==1:
    #         fits.writeto('{}hotpix_{}.fits'.format(path, filt[i]), \
    #                      CF_i[i], header=hdu_i[i][0].header, overwrite=True)

    #     # Print compilation time to bash:
    #     pt.compilation(i, n, 'Hot/Dead pixel removal')
    # print #----------------------------------------------------------------------
    # print('\n {} Hot/dead pixel found'.format(len((hot_pixels))))

    print(np.shape(hot_pixels))
    
z    # Plot if you like:
    if plot==1:
        pt.FITS(LF_i[0], 'linear', 2)
        plt.scatter(hot_pixels[0], hot_pixels[1], facecolors='none', edgecolors='r')
        plt.show()
        
    return CF_i 

#########################################################################################
#                                 STAR FINDER                                           #
#########################################################################################

def star_finder(pixel_array, min_pix=7, sigma=5, plot=1):       
    """
    Function to find the coordinates of the stars in the image
    ----------------------------
                INPUT          :
    ----------------------------
    path               (string): Path to data
    pixelArray         : Input array with pixel values
    sigma              : Number of standard deviations used as criteria to determine where the stars are
    min_pix            : Minimum number of pixels (above threshold) in a structure for it to count as a star
    ----------------------------
               OUTPUT          :
    ----------------------------
    CoM                : (y,x) coordinates of the Center of Mass
    radius             : An estimate of the radius in pixels
    N_stars            : Number of stars detected
    """
    # Height and width:
    h, w = np.shape(pixel_array)

    # FIND CENTER OF FLUX FOR STARS ABOVE THRESHOLD:
    # Define threshold as a number of standard deviations above the mean:
    threshold = np.mean(pixel_array) + sigma*np.std(pixel_array)
    # Find all pixels above the threshold:
    above_threshold = np.where(pixel_array > threshold, 1, 0)
    # Label the structures (where starmap = 1 that are adjacent to others):
    labels, N_stars = snd.label(above_threshold, structure = np.ones((3,3)))
    # Sum the number of elements in each structure:
    sums = snd.sum(above_threshold, labels, range(1,N_stars+1))
    # Choose only structures with more than min_pix elements (+1 index mismatch):
    stars = np.where(sums > min_pix)[0] + 1
    # Define starmap as 0 where there are no stars and 1 where there are stars:
    starmap = np.zeros(np.shape(pixel_array))
    for star in stars:
        starmap = starmap + np.where(labels == star, 1, 0)
    # Label all the structures again:
    labels, N_stars = snd.label(starmap, structure = np.ones((3,3)))
    # Find the center of mass of all the stars found above
    COF = snd.center_of_mass(pixel_array, labels, range(1,N_stars+1))
    # Estimate the radius of the star in pixels
    radius = np.sqrt(sums[stars-1]/np.pi)
    
    # FIND AND REMOVE ANY OVERLAPPING STARS (REMOVE THE SMALLER OF THE TWO):
    # Only do so if there are more than one star
    if len(COF) > 1:
        indices = check_distance(COF, COF, 20, mask=True)
        indices_to_remove = []
        for index in indices:
            if radius[index[0]] < radius[index[1]]:
                indices_to_remove.append(index[0])
            else:
                indices_to_remove.append(index[1])
        for index in sorted(indices_to_remove, reverse=True):
            del COF[index]
            radius = np.delete(radius, index)
        N_stars = N_stars - len(indices_to_remove)

    # Make as numpy array:
    COF = np.asarray(COF)

    # REMOVE STARS CLOSER THAN STELLAR STRUCTURE TO THE EDGE:
    # border = min_pix
    # i_x1 = np.where(COF[:,0]<border)[0]; i_x2 = np.where(COF[:,0]>w-border)[0]
    # i_y1 = np.where(COF[:,1]<border)[0]; i_y2 = np.where(COF[:,1]>h-border)[0]
    # i_xy = np.hstack([i_x1, i_x2, i_y1, i_y2])
    # X = np.delete(COF[:,0], i_xy)
    # Y = np.delete(COF[:,1], i_xy)
    
    #-------------------- PLOT ILLUSTRATION --------------------:
    if plot==1:
        plt.figure()
        pt.FITS(pixel_array, 'linear')
        plt.scatter(COF[:,1], COF[:,0], s=radius*6, facecolors='none', edgecolors='r')
        plt.title('{} stars are idenfied'.format(N_stars))
        plt.xlabel('X (pixels)'); plt.ylabel('Y (pixels)')
        plt.show()
    #-----------------------------------------------------------:        
    return COF, radius, N_stars


# This is a nother program to be tested!
def check_distance(array1, array2, threshold, mask=False, plot=1):
    """
    Function to check the distance between all the points in two arrays of coordinates
    """
    # Import needed package:
    from scipy.spatial.distance import cdist
    # cdist returns 1 if match or 0 if new star:
    dist = cdist(array1, array2)
    # 
    if mask is True:
        mdist = np.ma.masked_where(np.tril(dist)==0, dist)
        indices = np.nonzero(mdist < threshold)
    #
    elif mask is False: indices = np.nonzero(dist < threshold) 
    return np.transpose(indices)


#########################################################################################
#                                 MATCH STELLAR CATALOGS                                #
#########################################################################################

def match_coordinates(array1, array2, threshold=10, plot=0):
    """
    This function match two set of coordinates. This is done by a purely
    geometrical technique and looking at the histogram. It finds the
    minimum distance from i'th array1 star to every other array2 star.
    ----------INPUT:
    array1         : first catalog : Usually be catalog from the observation.
    array2         : second catalog: Usually be the CDS catalog.
    threshold      : Threshold with 10 pixels as default.
    plot           : plot=1 display the result.
    ---------OUTPUT:
    indices        : Rows with all the indices of matching
    """
    value_min = np.zeros(len(array1))
    index_min = {}

    # FIND MINIMUM DISTANCE WITH PYTHAGOREAN GEOMETRY:
    for i in range(len(array1)):
        d = np.sqrt( (array2[:,0]-array1[i,0])**2 + (array2[:,1]-array1[i,1])**2 )
        index_min[i] = np.argmin(d)
        value_min[i] = d[index_min[i]]
    # array1 indices of all stars:
    index_min = list(index_min.values())
    # find array1 stars within threshold:
    indices1 = np.where(value_min<threshold)[0]
    # Final list of matching array2 stars:
    indices2 = [index_min[i] for i in indices1]
    #-------------------- PLOT ILLUSTRATION --------------------
    if plot==1:
        # Plotting a histogram showing the best threshold:
        #pt.HIST(dis, 2000, 'Bins', '$\Delta N$', [0.0, 0.0025], [0, 1e3])
        # Plot coordinates and match:
        plt.scatter(array1[:,1], array1[:,0], marker='o', facecolors='none', edgecolors='r')
        plt.scatter(array2[:,1], array2[:,0], marker='o', facecolors='none', edgecolors='b')
        plt.scatter(array2[:,1][indices2], array2[:,0][indices2], marker='x', facecolors='k')
        plt.title('{} stars in common out of ({}, {}) available'.format(len(indices2), \
                                                                        len(array1), len(array2)))
        plt.xlabel('X (pixels)'); plt.ylabel('Y (pixels)')
        plt.show()
    #-----------------------------------------------------------
    return indices1, indices2



'''
Function to save a txt file with the x,y positions of the stars ordered after magnitude.
--------------INPUT:
Dataframe          : Dataframe with x and y positions with column names 'x' and 'y' sorted after magnitude or numpy array refCOF.
file_name          : Name of the output file
-------------OUTPUT:

---------------NOTE:
It saves a txt file with the positions that can be used on the Astrometry.net website to identify the stars.
It should also work with the COF of startrails.
'''

# def AstrometryNet (Input, file_name):
#     precision = 5
#     if isinstance(Input, pd.DataFrame):
#         with open(file_name, 'w') as text_file:
#             for row in Input.iterrows():
#                 print('%.*g\t%.*g' % (precision, row[1]['x'], precision, row[1]['y']), file = text_file)
#     elif isinstance(Input, np.ndarray):
#         with open(file_name, 'w') as text_file:
#             for row in Input:
#                 print('%.*g\t%.*g' % (precision, row[1], precision, row[0]), file = text_file)


#########################################################################################
#                                 IMAGE ALIGNMENT                                       #
#########################################################################################

def image_alignment(path, LF_name, method='starcoor', plot=0, save=1):
    """
    This function takes N number of frames and align them, thus make the post processing easier. This 
    function can both shift frames using stellar centroids.
    ----------INPUT:
    path           : Directory path to data
    LF_name        : Filename of Light Frame (LF)
    method         : Method to be used: 'centroid', 'fourier'
    mode           : Mode is used to align CRGB images just type 'CRGB'.
    plot           : Plot==1 shows 2D map of image positions relative to the first frame
    save           : save==1 saves a 3D cube of aligned frames
    ---------OUTPUT:
    LF_cube        : 3D cube of aligned frames
    """
    #-----------
    # LOAD DATA:
    #-----------
    LF_files = np.sort(glob.glob('{}{}*'.format(path, LF_name)))
    if len(LF_files)==0: print('ERROR: Wrong path or input name!'); sys.exit()
    n = len(LF_files)

    #----------------------
    # LOOP OVER ALL IMAGES:
    #----------------------
    for i in range(len(LF_files)):
        # Make string to save files:
        string_i = str(LF_files[i]).replace(LF_name, "final")
        # Open and close one image at a time:
        with fits.open(str(LF_files[i])) as hdu_i:
            LF_i = hdu_i[0].data

            #-----------------------
            #  STELLAR COORDINATES :
            #-----------------------
            if method=='starcoor':
                """
                This utility use the 'star_finder' to find all the stars within a serie of images. It then
                uses the identical stars identified in all images and compute the offset in (x, y) for all
                stars. The median offset value is the computed one image to the next. This list is the used
                to move the images with the 1st image as a default reference, but this is optional. 
                """
                # Use star_finder to find stellar coordinates:
                COF_i, R, N, = star_finder(LF_i, plot=plot)
                # First image is the reference:
                if i is 0:
                    COF_0 = COF_i
                    fits.writeto(string_i, LF_i, overwrite=True) 
                # MATCH COORDINATES AND SHIFT IMAGES:
                if i is not 0:
                    indices0, indices1 = match_coordinates(COF_0, COF_i, threshold=10, plot=plot)
                    shift = np.mean(COF_0[indices0] - COF_i[indices1], axis=0) 
                    LF_shift = snd.interpolation.shift(LF_i, shift)
                    fits.writeto(string_i, LF_shift, overwrite=True) 

            #-----------------------
            #         ALIPY        : TODO! do not work at the moment!
            #-----------------------
            if method=='alipy':
                """
                This utility uses....
                """
                # Import module:
                import alipy
                # Align image:
                identifications = alipy.ident.run(ref_image, images_to_align, visu=False)
                # Save aimge:
                fits.writeto(string_i, LF_shift, overwrite=True)
                
            #-----------------------
            #     TRIGONOMETRY     : TODO! Do not work properly at the moment. Hard to find the stars?
            #-----------------------
            if method=='astroalign':
                """
                This utility uses triangles made of stars to solve the alignment of the images.
                """
                # Import module:
                import astroalign
                # Align image:
                LF_shift = aa.register(LF_i[i], LF_i[0])
                # Save image:
                fits.writeto(string_i, LF_shift, overwrite=True) 

            #-----------------------
            #  IMAGE REGISTRATION  : (ONLY PYHTON 2.7)
            #-----------------------
            if method is 'imagereg':
                """
                This utility uses fourier analysis to shift the images. 
                """
                # Import module:
                from image_registration import chi2_shift
                from image_registration.fft_tools import shift
                # Align image:
                xoff, yoff, exoff, eyoff = chi2_shift(LF_i[0], LF_i[i])
                img_shift = shift.shiftnd(LF_i[i], -yoff, -xoff)
                LF_align[i] = img_shift
                # Save image:
                fits.writeto(string_i, LF_shift, overwrite=True)
             
        # Print to bash ------------------------------------------------------------------------
        pt.compilation(i, n, 'Image Alignment')
    print('')    
    return

#########################################################################################
#                                     RGB COMBINE                                       #
#########################################################################################

def rgb_combine(path, LF_name, scale, target, plot=0, save=1): 
    """
    This function takes RGB images and combie them to one sigle color image. All images are aligned before
    combined. One needs to specify which scale (linear, log, sqrt asinh) the image will be viewed in. The
    best result is obtained by performing image reduction (with flats, and darks) before using this 
    function. Each filter can be saved seperately in order to work with them in e.g. StarTools. 
    ----------INPUT:
    path           : Directory path to data
    R_name         : Filename of R filter images
    G_name         : Filename of G filter images
    B_name         : Filename of B filter images
    scale          : Scale used to plot images (linear, sqrt, log, asinh)
    plot           : Plot images: plot=1
    save           : Save images: save=1
    ---------OUTPUT:
    img            : Final color image 
    """
    # Load filenames:
    LF_files = np.sort(glob.glob('{}{}*'.format(path, LF_name)))
    # Contruct 3D cube of light frames:
    C = fits.getdata(str(LF_files[0]))
    R = fits.getdata(str(LF_files[1]))
    G = fits.getdata(str(LF_files[2]))
    B = fits.getdata(str(LF_files[3]))
    h, w = np.shape(C)
    # Defining how to scale the images (choosing min and max):
    ds = 2
    C_min, C_max = C.mean()-ds*C.std(), C.mean()+ds*C.std()
    R_min, R_max = R.mean()-ds*R.std(), R.mean()+ds*R.std()
    G_min, G_max = G.mean()-ds*G.std(), G.mean()+ds*G.std()
    B_min, B_max = B.mean()-ds*B.std(), B.mean()+ds*B.std()
    # Scale to view image:
    from Image_Scale import linear, sqrt, log, asinh
    if scale=='linear': scale=linear
    if scale=='sqrt':   scale=sqrt
    if scale=='log':    scale=log
    if scale=='asinh':  scale=asinh
    # Combining colors:
    img = np.zeros((h, w, 4))
    img[:,:,0] = scale(C, C_min, C_max)  # Clear
    img[:,:,1] = scale(R, R_min, R_max)  # Red
    img[:,:,2] = scale(G, G_min, G_max)  # Green
    img[:,:,3] = scale(B, B_min, B_max)  # Blue
    # Plot and save: 
    if plot==1: pylab.clf(); pylab.imshow(img, aspect='equal'); pylab.show()
    # Save png image:
    if save==1:imsave('{}CRGB_{}.png'.format(path, target), img)

    return

#########################################################################################
#                              APERTURE PHOTOMETRY                                      #
#########################################################################################


def aperture_photometry(path, LF_name, coor_stars, coor_target, aperture=None, dm_int=None,\
                        plot=None, save=None):
    """
    In short this software performs aperture photometry for every star that have a image
    coordinate asigned. The first coordinate needs to belong to the target star for which
    a final light curve is desired. However the flux and magnitude is calculated for every
    star, and all stars can be used as reference to correct for seeing variations during
    the total observation. Such a corrections is often refered to 'field (aperture) photometry'
    as the combined field of stars is used to correct for weather conditions during the night.
    The main function that needs to be called is 'aperture_photometry'.
    ----------------------------
    :           INPUT          :
    ----------------------------
    path             (string)  : Directory path to data
    LF_name          (string)  : Filename of Light Frame (LF)
    coor_stars       (2d array): List or name of file with stellar coordinates of reference stars: [x,y]
    apertures        (list)    : Apertures for star, inner and outer radius of background [r, r_min, r_max]
    dm_interval      (list)    : Limits for reference stars with lowest scatter [m_min, m_max, dm_min, dm_max]
    plot             (int)     : plot=1: timeseries; plot=2: coordinates; plot=3: optimal aperture
    save             (int)     : save==1: save timeseries [time, m]
    ----------------------------
    :          OUTPUT          :
    ----------------------------
    """
    # ----------
    # LOAD DATA:
    #-----------

    LF_files = np.sort(glob.glob('{}{}*'.format(path, LF_name)))
    LF0 = fits.getdata(str(LF_files[0]))
    # Header information:
    hdulist  = fits.open('{}'.format(str(LF_files[0])))
    n    = len(LF_files)
    h    = hdulist[0].header['NAXIS1']
    w    = hdulist[0].header['NAXIS2']
    filt = hdulist[0].header['INSFILTE']
    ron  = hdulist[0].header['RDNOISE']
    gain = hdulist[0].header['GAIN']
    # Constants:
    con  = 25
    R    = [aperture[1], aperture[2], aperture[3]]

    # Coordinates of target star:
    xs, ys = coor_target[0], coor_target[1]  
    # Coordinates of reference stars:
    if isinstance(coor_stars, str):
        coor_stars = np.loadtxt('{}{}.txt'.format(path, coor_stars))
        y, x = coor_stars[:,1], coor_stars[:,0] # Stellar coordinates
        N = len(x) # Number of stars (start)   
    else:
        x, y = coor_stars[:,1], coor_stars[:,0]
        N = len(x)       

    #---------------------
    # STELLAR COORDINATES:
    #---------------------
    
    # Make target star the first star:
    if coor_stars is not None:
        d = np.zeros(N)
        for i in range(N): d[i] = np.sqrt((xs-x[i])**2+(ys-y[i])**2) # Simple pythagorean geometry
        index_min = np.argmin(d)    # Find target star by geometry
        x = np.roll(x, -index_min)  # Shift data to target star index
        y = np.roll(y, -index_min)
        # Check if stellar coordinates are closer than r_max:
        i_x1 = np.where(x<R[2])[0]; i_x2 = np.where(x>w-R[2])[0]
        i_y1 = np.where(y<R[2])[0]; i_y2 = np.where(y>h-R[2])[0]
        i_xy = np.hstack([i_x1, i_x2, i_y1, i_y2])
        # Discard these coordinates:
        X = np.delete(x, i_xy)
        Y = np.delete(y, i_xy)
        N = len(x)                  # New number of stars
        # Plot if you like:
        if plot==2:
            pt.FITS(LF0, 'linear')
            plt.plot(y, x, 'ro', mfc='none', label='Ref stars')
            plt.plot(y[0], x[0], 'b+', label='Target star')
            plt.plot(y[i_xy], x[i_xy],'y*', label='Discarded')
            plt.legend(loc='lower right')
            plt.show()
        
    # If only the target star:
    if coor_stars is None: 
        X, Y = [xs], [ys]
        N = 1
            
    #------------------------------
    # FIND OPTIMAL APERTURE RADIUS:
    #------------------------------ 
    
    if plot==3 or aperture[0]==None:
        SNR_r = optimal_aperture(LF0, y[0], x[0], R, gain, ron)
    
    #------------
    # PHOTOMETRY:
    #------------

    flux_sky, n_star_pix, flux_star = aperture(LF0, X[0], Y[0], R, plot)
    sys.exit()
    # Find fluxes:
    flux_star = np.zeros((n,N))
    SNR_i     = np.zeros((n,N))
    time      = np.zeros(n) 
    # Loop over all images:
    for i in range(n):
        hdulist = fits.open(str(LF_files[i]))
        time[i] = Time(hdulist[0].header['DATE'], scale='utc').jd 
        LF = fits.getdata(str(LF_files[i]))
        print(np.shape(LF))
        # Loop over all stars and find flux (using same radius):
        for j in range(N):
            flux_sky, n_star_pix, flux_star[i][j] = aperture(LF, X[j], Y[j], R, plot)
                
            print(x[j])
            sys.exit()
            # flux_sky, n_pix_star, flux_star[i][j] = aperture(LF, x[j], y[j], R, plot)
            # SNR_i[i][j] = SNR(flux_sky, n_pix_star, flux_star[i][j], gain, ron)

        # Print to bash:
        plots.compilation(i, n, '')
    print
            
    # print time
    # print np.shape(flux_star)

    # plt.plot(time, flux_star[:,1], '.'); plt.show()
    # sys.exit()

    #-------------------
    # SEEING CORRECTION:
    #-------------------

    mag_star = -2.5*np.log10(flux_star[:,0])+con          # From flux to magnitude for target star

    print(mag_star)
    print(time)
    # sys.exit()
    
    # if coor_stars is not None:
    #     mag_ref  = -2.5*np.log10(flux_star[:,range(1,N)])+con # From flux to magnitude for reference stars
    #     # Find all delta-mag (dm):
    #     #print mag_star; sys.exit()
    #     dm  = np.zeros((n,N-1))                             
    #     for i in range(n):
    #         dm[i] = mag_ref[0,:]-mag_ref[i,:]              # dm from 0th frame to every other
    #     # Make trumpet plot:
    #     dm_med = np.median(dm,axis=0)
    #     # Using all stars:
    #     mag_all = mag_star + np.median(dm, axis=1)            # Using all stars

    #     # Which stars should be used to correct with:
    #     if dm_int==None:
    #         dm_int = [-0.01, 0.01, 3, 10]
    
    #     # Using interval:
    #     im_x = np.where((mag_ref[0,:]>dm_int[0])*(mag_ref[0,:]<dm_int[1]))[0]
    #     im_y = np.where((dm_med>dm_int[2])*(dm_med<dm_int[3]))[0]
    #     im_xy = list(set(im_x).intersection(im_y))
    #     mag_int = mag_star + np.median(dm[:,im_xy], axis=1)

    # # Plot results:
    # if coor_stars is not None:
    #     if plot==1 or plot==2:
    #         plots.plot_trumpet(mag_ref, dm_med, dm_int, im_xy, '$m_i$', '$m_0\, -\, m_i$',\
    #                            'Trumpet - Seeing Corretion')
    #         plt.show()
    #         plots.MAG(time, [mag_star, mag_all, mag_int], ['.', 'r.', 'g*'], '$t$ (hours) - UTC JD', '$m$')
    #         plt.show()

    if coor_stars is None and plot==1 or plot==2:
        plt.plot(time, mag_star, '+')#, '$t$ (hours) - UTC JD', '$m$')
        plt.show()
    return



def aperture(LF, x, y, R, plot):
    """
    This function find stellar flux with a circular aperture.
    ----------INPUT:
    LF             : Light Frame (LF)
    x, y           : Stellar coordinate
    R              : [r, r_min, r_max]: aperture, inner and outer radius for background
    plot           : Plot if you like. 
    ----------------
    """
    
    # from numpy import meshgrid, ogrid, sqrt, mean, median, argmax, sum, cos, sin, nonzero
    # Split varibles:
    r, r_min, r_max = R[0], R[1], R[2]
    x, y = int(x), int(y)
    # Make meshgrid maks:
    x_grid, y_grid = np.meshgrid(range(x-r_max,x+r_max), range(y-r_max,y+r_max))
    square = LF[y_grid, x_grid]                       # Square frame of each image
    c      = np.sqrt((x_grid-x)**2 + (y_grid-y)**2)   # Equation of a circle
    # Background:
    c_sky = ((c>=r_min)*(c<=r_max))*square            # r_max^2 dimensional array
    sky_mean = np.mean(c_sky[np.nonzero(c_sky)])      # Mean sky flux
    sky_medi = np.median(c_sky[np.nonzero(c_sky)])    # Median sky flux
    flux_sky = 3*sky_medi - 2*sky_mean                # Sky flux
    # Star:
    c_star     = (c<=r)*square
    star       = c_star[np.nonzero(c_star)]-flux_sky  # Stellar corrected pixels
    n_pix_star = np.sum(c<=r)                         # Number of used star pixels
    flux_star  = np.sum(star)                         # Flux from star

    # Plot zoom-in on star with aperture: 
    if plot==3:
        plots.FITS(square, 'linear')
        t = np.linspace(0,2*np.pi)
        x_off = (2*r_max+1)/2
        y_off = (2*r_max+1)/2
        plt.plot(x_off,y_off,'+r')
        plt.plot(r*np.cos(t)+x_off, r*np.sin(t)+y_off, 'g-')
        plt.plot(r_min*np.cos(t)+x_off, r_min*np.sin(t)+y_off, 'b-')
        plt.plot(r_max*np.cos(t)+x_off, r_max*np.sin(t)+y_off, 'b-')
        plt.show() # Shows only figure ones in loop
        sys.exit()
        
    return flux_sky, n_pix_star, flux_star
  


def SNR(flux_sky, n_pix_star, flux_star, gain, ron):
    """
    This function calculates the Signal-to-Noise Ratio (SNR). Here: 
    flux_sky       : Flux from each sky pixel
    n_pix_star     : Number of pixels used to find star flux
    flux_star      : Flux from stellar object
    gain           : Gain CCD (e-/ADU)
    ron            : Read out noise (ron) (e-)
    """
    SNR = (gain*flux_star/np.sqrt(gain*flux_star + n_pix_star*gain*flux_sky + n_pix_star*ron**2))
    return SNR



def optimal_aperture(LF, x, y, R, gain, ron):
    """
    This function can be used to find the aperture where the SNR is highest, hence, the best aperture. 
    """
    # Constants:
    r_min, r_max = R[1], R[2]
    r_i = range(r_max)
    # Find SNR as a funcion of aperture radius:
    SNR_r = np.zeros(r_max)
    for r in r_i:
        x, y = int(x), int(y)
        # Make meshgrid maks:
        x_grid, y_grid = np.meshgrid(range(x-r_max-2,x+r_max+2), range(y-r_max-2,y+r_max+2))
        square = LF[y_grid, x_grid]                    # Square frame of each image
        c      = np.sqrt((x_grid-x)**2 + (y_grid-y)**2)   # Equation of a circle
        # Background:
        c_sky = ((c>=r_min)*(c<=r_max))*square         # r_max^2 dimensional array
        sky_mean = np.mean(c_sky[np.nonzero(c_sky)])         # Mean sky flux
        sky_medi = np.median(c_sky[np.nonzero(c_sky)])       # Median sky flux
        flux_sky = 3*sky_medi - 2*sky_mean             # Sky flux
        # Star:
        c_star     = (c<=r)*square
        star       = c_star[np.nonzero(c_star)]-flux_sky  # Stellar corrected pixels
        n_pix_star = np.sum(c<=r)                         # Number of used star pixels
        flux_star  = np.sum(star)                         # Flux from star
        # SNR:
        SNR_r[r] = SNR(flux_sky, n_pix_star, flux_star, gain, ron)
    # Plotting:
    plt.plot(range(r_max), SNR_r, '*')
    plt.xlabel('Aperture Radius')
    plt.ylabel('SNR')
    plt.show()   
    return SNR_r


def center_of_flux(LF, x, y, n_pix, R):
    """
    This function finds the center of flux for all desired stellar object. Here LF is the masked image
    thus every pixel is set to zero except the star) and n_pix is the number of pixels one wish to use
    in order to find the COF.
    """
    # Cut out part of image that i relevant:
    x_grid, y_grid = np.meshgrid(range(x-R,x+R), range(y-R,y+R))
    LF = LF[y_grid, x_grid]      # Square frame of each image

    # plots.FITS(LF, 'linear', 2)
    # plt.show()
    # sys.exit()

    # Loops over all pixels:
    LF_copy  = np.copy(LF)     # Copy to avoid overwriting
    flux_max = np.zeros(n_pix)
    x_max = np.zeros(n_pix)
    y_max = np.zeros(n_pix)
    pixel = np.zeros(n_pix)
    for j in range(n_pix):
        flux_max[j] = np.max(LF_copy)               # Maximum value for array
        max_dex = np.where(LF_copy == flux_max[j])  # Find row, column for min value
        x_max[j] = max_dex[0][0]                    # max for x coordinate
        y_max[j] = max_dex[1][0]                    # max for y coordinate
        pixel[j] = j
        # Min pixel is et to max in order to find the next min:
        LF_copy[int(y_max[j]), int(x_max[j])] = 0 

    # Flux center is found:
    flux = sum(pixel)
    cen_x = 1/flux*np.dot(pixel,x_max)              # Flux-center in x
    cen_y = 1/flux*np.dot(pixel,y_max)              # Flux-center in y

    return int(cen_x)+x-2*R, int(cen_y)+y-2*R

# def center_of_flux_old(LF, x, y, R):
#     """This function finds the center of flux for all desired stellar object"""
#     # Packages:
#     from numpy import zeros, max, sum, round, meshgrid, argmax
#     # Extract data:
#     r, r_max = R[0], R[2]
#     x, y     = x.astype(int), y.astype(int)
#     # Outer: loops over all stars:
#     center_x = zeros(len(x))
#     center_y = zeros(len(y))
#     for i in range(len(x)):
#         x_grid, y_grid = meshgrid(range(x[i]-r_max,x[i]+r_max), range(y[i]-r_max,y[i]+r_max))
#         square = LF[y_grid, x_grid]      # Square frame of each image
#         # Test:
#         FITS(square,'linear'); plt.plot(x_grid-x[i]+r_max, y_grid-y[i]+r_max, 'b*'); plt.show()
#         # Inner: loops over all pixels: 
#         for j in range(2*r):
#             x_max_i = max(square)                       # Row and value of max pixel
#             pixel_i = square.index(x_max_i)             # Index of max value pixel
#             _, y_max_i       = max(square[x_max_i, :])   # Column for max pixel
#             x_max[j]         = x_max_i                   # Rows gathers in a list
#             y_max[j]         = y_max_i                   # Columns gathers in a list
#             pixel[j]         = pixel_i                   # Max pixels gathers in a list
#             square[x_max_i, y_max_i] = 0                 # Max pixel is set to zero to find the nextr
#         # Flux center is found:
#         flux      = sum(pixel)                           # Total flux from pixels
#         center_x  = 1/flux*dot(pixel,x_max)              # Flux-center in x
#         center_y  = 1/flux*dot(pixel,y_max)              # Flux-center in y
#     # Center of fluc for all stars:
#     cen = [center_x-r, center_y-r]
#     return cen


# def center_of_flux(LF, x, y, R):
#     """This function finds the center of flux for all desired stellar object"""
#     # Packages:
#     from numpy import zeros, max, sum, round, meshgrid, argmax
#     # Extract data:
#     r, r_max = R[0], R[2]
#     # x, y     = x.astype(int), y.astype(int)
#     # Outer: loops over all stars:
#     center_x = zeros(len(x))
#     center_y = zeros(len(y))
#     for i in range(len(x)):
#         x_grid, y_grid = meshgrid(range(x[i]-r_max,x[i]+r_max), range(y[i]-r_max,y[i]+r_max))
#         square = LF[y_grid, x_grid]      # Square frame of each image
#         # Test:
#         # FITS(square,'linear'); plt.plot(x_grid-x[i]+r_max, y_grid-y[i]+r_max, 'b*'); plt.show()
#         # Inner: loops over all pixels: 
#         for j in range(2*r):
#             x_max_i = max(square)                       # Row and value of max pixel
#             pixel_i = square.index(x_max_i)             # Index of max value pixel
#             _, y_max_i       = max(square[x_max_i, :])   # Column for max pixel                   
#             x_max[j]         = x_max_i                   # Rows gathers in a list
#             y_max[j]         = y_max_i                   # Columns gathers in a list  
#             pixel[j]         = pixel_i                   # Max pixels gathers in a list
#             square[x_max_i, y_max_i] = 0                 # Max pixel is set to zero to find the nextr
#         # Flux center is found:
#         flux      = sum(pixel)                           # Total flux from pixels
#         center_x  = 1/flux*dot(pixel,x_max)              # Flux-center in x
#         center_y  = 1/flux*dot(pixel,y_max)              # Flux-center in y
#     # Center of fluc for all stars:
#     cen = [center_x-r, center_y-r]
#     return cen


# def match_coordinates(array1, array2, threshold=5, plot=1):
#     """
#     This function match two set of coordinates. This is done by a purely geometrical technique and looking 
#     at the histogram. Find the minimum distance from i'th cat0 star to every other cat1 star.
#     ----------INPUT:
#     cat0           : first catalog : Usually be catalog from the observation.
#     cat1           : second catalog: Usually be the CDS catalog.
#     TOL            : Tolerance of distance between catalogs.
#     plot           : plot=1 display the result.
#     ---------OUTPUT:
#     indices        : Two column array with all the indices of matching stars.
#     """
#     # Constants:
#     aconst    = len(array1)
#     d         = np.zeros(aconst)
#     value_min = np.zeros(aconst)
#     index_min = {}
#     # FIND MINIMUM DANCE WITH PYTHAGOREAN GEOMETRY: 
#     for i in range(aconst):
#         d[i] = np.min(np.sqrt( (array2[:,0]-array1[:,0][i])**2 + (array2[:,1]-array1[:,1][i])**2) )
#         index_min[i] = np.argmin(d)
#         value_min[i] = d[index_min[i]]
#     # FIND ALL CAT0 OBJECTS THAT MATCH CAT1 OBECTS WITHIN THRESHOLD:
#     # When more stars are found within TOL, the one with min distance is used.
#     # cat1 indices of all stars:
     
#     index_min = list(index_min.values())
#     # find cat1 stars within threshold:
#     index_cat1 = np.where(value_min<threshold)[0]
#     # Final list of matching cat1 stars:
#     index_cat2 = [index_min[i] for i in index_cat1]

#     #-------------------- PLOT ILLUSTRATION --------------------#
#     if plot==1:
#         # Plotting a histogram showing the best threshold:
#         #pt.HIST(dis, 2000, 'Bins', '$\Delta N$', [0.0, 0.0025], [0, 1e3])

#         # Plot coordinates and match:
#         plt.scatter(array1[:,1], array1[:,0], marker='o', facecolors='none', edgecolors='r')
#         plt.scatter(array2[:,1], array2[:,0], marker='o', facecolors='none', edgecolors='b')
#         plt.scatter(array1[:,1][index_cat2], array1[:,0][index_cat2], marker='x', facecolors='k')
#         plt.title('{} stars in common out of ({}, {})'.format(len(index_cat2), aconst, len(array2)))
#         plt.xlabel('X (pixels)'); plt.ylabel('Y (pixels)')
#         plt.show()
#     #-----------------------------------------------------------#
#     sys.exit()
#     return index_cat2



    def image_reduction_fits(self, LF_name, FF_name, DF_name=None, BF_name=None, smooth=None, align=None, \
                         plot=None, save=None):
        """ This function takes all loaded flat-frames, dark-frames, bias-frames and combine them to a one
        master-flat-field image.
        ----------INPUT:
        path           : Directory path to data.
        LF_name        : Name of Light Frames (LF) except end number.
        FF_name        : Name of Flat  Frames (FF) except end number.
        DF_name        : Name of Dark  Frames (DF) except end number.
        BF_name        : Name of Bias  Frames (BF) except end number.
        N              : Number of [LF, FF, DF, BF] frames.
        plot           : If plot=1: plots all relevant frames..
        save           : If save=1: saves all corrected light frames in seperate files.
        ---------OUTPUT:
        CF_i           : Cube of Corrected Frames (CF)."""
        print '--------------------------------------------------------------------- image_reduction_fits'
        start_time = time.time()  # Take time
        
        # Construct 3D cube of light-frames and flat-frames with dim(number, x-dim, y-dim):
        LF_i = array([pyfits.getdata('{}{}_%03d.fits'.format(path, LF_name) %i) for i in range(N[0])])
        FF_i = array([pyfits.getdata('{}{}_%03d.fits'.format(path, FF_name) %i) for i in range(N[1])])
        FF   = median(FF_i, axis=0)    # Median correction image:
        
        # Uses Flats, Darks and Bias':
        if DF_name!=None and BF_name!=None:
            DF_i = array([pyfits.getdata('{}{}_%03d.fits'.format(path, DF_name) %i) for i in range(N[2])])
            BF_i = array([pyfits.getdata('{}{}_%03d.fits'.format(path, BF_name) %i) for i in range(N[3])])
            # Median correction images:
            DF   = median(DF_i, axis=0)
            BF   = median(BF_i, axis=0)
            # Perform correction: (Image Arithmetric Operation)
            CF_i = (LF_i - DF)/(FF - BF)*mean(FF - BF) 

        # Uses Flats and Darks:
        elif DF_name!=None:
            DF_i = array([pyfits.getdata('{}{}_%03d.fits'.format(path, DF_name) %i) for i in range(1, N[2]+1)])
            DF   = median(DF_i, axis=0)  # Median correction image:
            CF_i = (LF_i - DF)/FF        # Perform correction:
        # Uses Flats only:
        else:
            CF_i = LF_i/FF               # Perform correction:
        

        def smooth_image(self, data, smooth):
        """ This function takes all loaded flat-frames, dark-frames, bias-frames and combine them to a one
        master-flat-field image.
        ----------INPUT:
        path           : Directory path to data.
        LF_name        : Name of Light Frames (LF) except end number.
        FF_name        : Name of Flat  Frames (FF) except end number.
        DF_name        : Name of Dark  Frames (DF) except end number.
        BF_name        : Name of Bias  Frames (BF) except end number.
        N              : Number of [LF, FF, DF, BF] frames.
        plot           : If plot=1: plots all relevant frames..
        save           : If save=1: saves all corrected light frames in seperate files.
        ---------OUTPUT:
        CF_i           : Cube of Corrected Frames (CF)."""
        print '--------------------------------------------------------------------- smooth_image'

        # Data definition:
        CF_i = data
        
        # Clean data for hot pixels:
        if smooth=='median':
            HF_i = copy(CF_i)
            HF_i = median_filter(HF_i, size=3)

        # Replace hot and dead pixels only:
        if smooth=='hot':
        HF_i = copy(CF_i)
        n, h, w = shape(CF_i) # Number of frames, Height and Width
        for i in range(n):
            # Find good threshold to locate pixels:
            blur_i = median_filter(CF_i[i], size=2)
            diff_i = CF_i[i] - blur_i
            threshold = 2*std(diff_i)
            # Find the hot pixels (ignoring the edges):
            hot_pixels = nonzero((abs(diff_i[1:-1, 1:-1])>threshold))
            hot_pixels = array(hot_pixels) + 1  # +1 because 1. row and column was ignored
            print 'Frame {}: Hot/Dead pixels: {}'.format(i, len(hot_pixels[1]))
            for y,x in zip(hot_pixels[0], hot_pixels[1]):
                HF_i[i,y,x] = blur_i[y,x]
            # Now get the pixels on the edges (but not the corners):
            # Left and right sides:
            for index in range(1, h-1):
                # Left side:
                med  = median(CF_i[i, index-1:index+2, 0:2])
                diff = abs(CF_i[i, index, 0] - med)
                if diff>threshold: 
                    hot_pixels = hstack((hot_pixels, [[index],[0]]))
                    HF_i[i, index, 0] = med
                # Right side:
                med  = median(CF_i[i, index-1:index+2, -2:])
                diff = abs(CF_i[i, index,-1] - med)
                if diff>threshold: 
                    hot_pixels = hstack((hot_pixels, [[index], [w-1]] ))
                    HF_i[i, index, -1] = med
            # Then the top and bottom:
            for index in range(1, w-1):
                # Bottom:
                med  = median(CF_i[i, 0:2, index-1:index+2])
                diff = abs(CF_i[i, 0, index] - med)
                if diff>threshold: 
                    hot_pixels = hstack((hot_pixels, [[0], [index]] ))
                    HF_i[i, 0, index] = med
                # Top:
                med  = median(CF_i[i, -2:, index-1:index+2])
                diff = abs(CF_i[i, -1, index] - med)
                if diff>threshold: 
                    hot_pixels = hstack((hot_pixels, [[h-1], [index]] ))
                    HF_i[i, -1, index] = med
            # Then the corners:
            # Bottom left:
            med  = median(CF_i[i, 0:2, 0:2])
            diff = abs(CF_i[i, 0, 0] - med)
            if diff>threshold: 
                hot_pixels = hstack((hot_pixels, [[0], [0]] ))
                HF_i[i, 0, 0] = med
            # Bottom right:
            med  = median(CF_i[i, 0:2, -2:])
            diff = abs(CF_i[i, 0, -1] - med)
            if diff>threshold: 
                hot_pixels = hstack((hot_pixels, [[0], [w-1]] ))
                HF_i[i, 0, -1] = med
            # Top left:
            med  = median(CF_i[i, -2:, 0:2])
            diff = abs(CF_i[i, -1, 0] - med)
            if diff>threshold: 
                hot_pixels = hstack((hot_pixels, [[h-1], [0]] ))
                HF_i[i, -1, 0] = med
            # Top right:
            med  = median(CF_i[i, -2:, -2:])
            diff = abs(CF_i[i, -1, -1] - med)
            if diff>threshold: 
                hot_pixels = hstack((hot_pixels, [[h-1], [w-1]] ))
                HF_i[i, -1, -1] = med

    # Align frames:
    if align!=None:
        print 'hej'
                
    # Plot if you like:
    if plot==1:
        if FF_name!=None: FITS(FF, 'linear', 5) 
        if DF_name!=None: FITS(DF, 'linear')  
        if BF_name!=None: FITS(BF, 'linear', 5)    
        FITS(LF_i[0], 'linear', 2.5)
        FITS(CF_i[0], 'linear', 2)
        if smooth=='median': FITS(HF_i[0], 'linear', 2.5)
        if smooth=='hot':    FITS(HF_i[0], 'linear', 2)
  
    # Save frames in new fits files:
    if save==1:
        if smooth!=None:
            for i in range(1, N[0]+1):
                pyfits.writeto(('{}{}_%03d.fits'.format(path, 'CF') %i), HF_i[i-1], clobber=True)
        else:
            for i in range(1, N[0]+1):
                pyfits.writeto(('{}{}_%03d.fits'.format(path, 'CF') %i), CF_i[i-1], clobber=True)
                
    print ('Filter done in time: %s s' % (time.time() - start_time))
    return CF_i 




def image_reduction_jpg(path, N, LF_name, FF_name=None, DF_name=None, BF_name=None, plot=None, save=None):
    """------------------------------------------ FUNCTION ----------------------------------------------:
    This function takes all loaded flat-frames, dark-frames, bias-frames and combine them to a one
    master-flat-field image.
    ----------INPUT:
    path           : Directory path to data.
    LF_name        : Name of Light Frames (LF) except end number.
    FF_name        : Name of Flat  Frames (FF) except end number.
    DF_name        : Name of Dark  Frames (DF) except end number.
    BF_name        : Name of Bias  Frames (BF) except end number.
    N              : Number of [LF, FF, DF, BF] frames.
    plot           : If plot=1: plots all relevant frames..
    save           : If save=1: saves all corrected light frames in seperate files.
    ---------OUTPUT:
    CF_i           : Cube of Corrected Frames (CF)."""
    print '--------------------------------------------------------------------- image_reduction_jpg'
    start_time = time.time()  # Take time

    # Load data:
    LF_i = array([asarray(Image.open('{}{}_%02d.jpg'.format(path, LF_name) %i)) for i in range(N[0])])
    if FF_name!=None:
        FF_i = array([asarray(Image.open('{}{}_%02d.jpg'.format(path, FF_name) %i)) for i in range(N[1])])
    if DF_name!=None:
        DF_i = array([asarray(Image.open('{}{}_%02d.jpg'.format(path, DF_name) %i)) for i in range(N[2])])
    if BF_name!=None:
        BF_i = array([asarray(Image.open('{}{}_%02d.jpg'.format(path, BF_name) %i)) for i in range(N[3])])
        
    # Different statements:
    if FF_name!=None:
        CF_i = LF_i/mean(FF_i, axis=0)
        
    if DF_name!=None:
        CF_i = LF_i - mean(DF_i, axis=0)
        
    if BF_name!=None:
        CF_i = LF_i - mean(BF_i, axis=0)

    if FF_name!=None and DF_name!=None:
        print 'Using FLATS and DARKS'
        # Darks:
        DF_master = mean(DF_i, axis=0)
        # Flats:
        FF_i = FF_i/mean(FF_i, axis=0)   # Normalize
        FF_master = mean(FF_i, axis=0)
        # Corrections:
        CF_i = (LF_i - BF_master)/FF_master

    
    if FF_name!=None and BF_name!=None:
        print 'Using FLATS and BIAS'
        # Bias:
        BF_master = median(BF_i, axis=0)
        # Flats:
        FF_bias = FF_i - BF_master          # Bias correction
        FF_medi = median(FF_bias, axis=0)
        FF_master = FF_medi/median(FF_i, axis=0)   # Normalize
        # Corrections:
        CF_i = LF_i - BF_master - FF_master
        
        
    # Print variables:
    if FF_name!=None: print('FF (mean: %.2f) (std: %.2f)' % (mean(FF_i), std(FF_i)))
    if DF_name!=None: print('DF (mean: %.2f) (std: %.2f)' % (mean(DF_i), std(DF_i)))
    if BF_name!=None: print('BF (mean: %.2f) (std: %.2f)' % (mean(BF_i), std(BF_i)))

    # Plot if you like:
    if plot==1:
        
        # Plot Bias:
        if BF_name!=None:
            # Join axes:
            fig, axes = plt.subplots(1,2)
            # Plot data:
            for ax, data in zip(axes.flat, [BF_i[0], BF_master]):
                img = ax.imshow(data)
            #fig.subplots_adjust(right=1.2)
            cbar_ax = fig.add_axes([0.1, 0.05, 0.8, 0.05])
            fig.colorbar(img, cbar_ax, orientation='horizontal')
            plt.show()
            sys.exit()
        # Plot flats:
        if FF_name!=None:
            # Join axes:
            fig, axes = plt.subplots(1,2)
            # Plot data:
            for ax, data in zip(axes.flat, [FF_i[0], FF_master]):
                img = ax.imshow(data)
            #fig.subplots_adjust(right=0.8)
            cbar_ax = fig.add_axes([0.7, 0.15, 0.05, 0.7])
            fig.colorbar(img, cbar_ax)
            plt.show()
            
        plt.imshow(FF_i[0]); plt.show()
        plt.imshow(FF_master); plt.show()
        plt.imshow(LF_i[0]); plt.show()
        plt.imshow(CF_i[0]); plt.show()

    if save==1:
        for i in range(N[0]):
            pyfits.writeto(('{}{}_%03d.fits'.format(path, 'CF') %i), CF_i[i], clobber=True)


