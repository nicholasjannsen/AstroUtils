"""
SOFTWARE DESCRIPTION:
---------------------

Written April 2018 -- Nicholas Jannsen
Typeset in Python 3

This routine is developed such with all the necessary calibration images it returns an LRGB image, and also a master image for each filter.
"""
import sys, time
import numpy as np

# Import reduction routine:
sys.path.append('../')
import Astro_Filters as astro

##########################################################################################################
#                                               PIPELINE                                                 #
##########################################################################################################

# def rgb_combine(path, target, plot=None, save=None):

#     #-----------------
#     # IMAGE REDUCTION:
#     #-----------------
#     print('1/3: Image reduction')

#     # Calibrate each filter and save images:
#     astro.image_reduction(path, 'LF_C', 'FF_C', 'DF', 'BF', plot, save)
#     astro.image_reduction(path, 'LF_R', 'FF_R', 'DF', 'BF', plot, save)
#     astro.image_reduction(path, 'LF_G', 'FF_G', 'DF', 'BF', plot, save)
#     astro.image_reduction(path, 'LF_B', 'FF_B', 'DF', 'BF', plot, save) 

#     #------------------------
#     # HOT/DEAD PIXEL REMOVAL:
#     #------------------------
#     print('2/3: Hot/dead pixel removal')

#     # Removing hot/dead pixels and save images to calibration files:
#     astro.pixel_outliers(path, 'CLF_C', 'hot', None, plot, save)
#     astro.pixel_outliers(path, 'CLF_R', 'hot', None, plot, save)
#     astro.pixel_outliers(path, 'CLF_G', 'hot', None, plot, save)
#     astro.pixel_outliers(path, 'CLF_B', 'hot', None, plot, save)
    
#     #-----------------
#     # IMAGE ALIGNMENT:
#     #-----------------
#     print('3/3: Image alignment')
    
#     # Aligning each filter and save shifted images:
#     astro.image_alignment(path, 'CLF_C', 'astroalign', None, plot, save)
#     astro.image_alignment(path, 'CLF_R', 'astroalign', None, plot, save)
#     astro.image_alignment(path, 'CLF_G', 'astroalign', None, plot, save)
#     astro.image_alignment(path, 'CLF_B', 'astroalign', None, plot, save)

#     # Align each master image:
#     print('Align each filter')
#     astro.image_alignment(path, 'ACLF', 'astroalign', 'LRGB', plot, save)
    
#     # Preview result: (Use Star Tools for further analysis)
#     astro.LRGB_combine(path, 'ACLF', 'linear', target, plot, save)



##########################################################################################################
#                                     Different tasks                                                    #
##########################################################################################################

def test(name):

     if name=='corot1b': # This is an exoplanet I observed with the IAC80 at Tenerife:
          path = '/home/nicholas/Data/IAC80/2015-01-06/corot1b/'
          #star_finder(path, 'red', 5, 1, 1)
          astro.aperture_photometry(path, 'red', 'star_coor', [1010, 981], ['circle', 10, 15, 18],\
                                    [5, 12, -0.01, 0.01], 1, 0)
     
     if name=='34leo': # Stellar magnetic active stars from the HK survey from Mount Wilson:
          path = '/home/nicholas/Data/ORO/28cm/34leo/'
          # phot.image_reduction(path, '34leo', 'flat', 'dark', 'bias', 1, 1)
          # phot.align_images(path, '34leo0', 1, 1)
          # phot.star_finder(path, '34leo', 100, 1, 1)
          astro.aperture_photometry(path, ['34leo0', '34leo0'], 'star_coor', \
                                   [1364, 1607], [17, 30, 40], [0, 20, -10, 10], 0, 0)
          
     if name=='ebs':
          """ EB: HD187276 """
          path = '/home/nicholas/Data/ORO/28cm/HD187276/'
          astro.image_reduction(path, 'LF', 'FF', 'DF', 'BF', 0, 0)
          # phot.align_images(path, '34leo0', 1, 1)
          # phot.star_finder(path, '34leo', 100, 1, 1)
          # phot.aperture_photometry(path, ['34leo0', '34leo0'], 'star_coor', \
          #                          [1364, 1607], [17, 30, 40], [0, 20, -10, 10], 1, 0)
     
     if name=='bet_Cas': # Delphini-1 camera test 2018
          path = '/home/nicholas/Data/ORO/28cm/bet_Cas_test/'
          # star_coor =
          astro.image_alignment(path, 'LF', 'astroalign', None, 1, 0)
          # astro.aperture_photometry(path, 'ALF', None, [948, 845], [100,160,170], None, 1, 0)
     
          
     # if name=='photometry_class':
     #      from Photometry_Class import Photometry_Class
     #      path = '/home/nicholas/Data/IAC80/2015-01-06/corot1b/'
     #      XX = Photometry_Class(path, 'red', 3, 0)
     #      XX.aperture_photometry('star_coor', [1010, 981], [10,15,18], [5,15,-0.01,0.01], 2, 0)

          
if __name__ == '__main__': 
     # FUNCTION CALL:
     test(name='corot1b')    
