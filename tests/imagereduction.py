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


