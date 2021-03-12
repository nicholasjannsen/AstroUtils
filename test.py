#o!/Usr/Bin/Python
#-*- coding: utf-8 -*-

# Numpy:
import numpy as np
from numpy import inf, nan, sin, cos, pi, sqrt, diff, std, diag, log, mean, std, asarray
from numpy import mean, median, nanargmax, zeros, ones, ceil, delete
from numpy import arange, array, size, vstack, copy, loadtxt, where, savetxt
from numpy import min, max, sum, float, round, int
# Others:
import math, os, sys, time, random
# Plotting:
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
import os
from Plot_Tools import FITS
from PIL import Image

##################################################################################################

def test(name):
#------------------------------------------- FUNCTION -------------------------------------------:
# This function calls a subfunction by name and execute it.
#------------INPUT:
# name            : Name of called function.
#------------------------------------------------------------------------------------------------:      

     if name=='photometry':
          from Photometry import image_reduction, star_finder, aperture_photometry
          """ 34 Leo """
          path = '/home/nicholas/Data/ORO/28cm/34leo/'
          image_reduction(path, '34leo', 'flat', 'dark', 'bias', 1, 0)
          #star_finder(path, '34leo', 100, 1, 1)
          aperture_photometry(path, ['34leo', '34leo_align'], 'star_coor', \
                              [1607, 1364], [20,30,40], [5,12,-0.01,0.01], 2, 0)
          sys.exit()
          """ Corot-1b """
          path = '/home/nicholas/Data/IAC80/2015-01-06/corot1b/'
          #star_finder(path, 'red', 5, 1, 1)
          aperture_photometry(path, ['red', 'red_align'], 'star_coor', \
                              [1010, 981], [10,15,18], [5,12,-0.01,0.01], 1, 0)
          #----- Test class Phtometry:
          """ HD 27130 """
          #path = '/home/nicholas/Data/ORO/50cm/HD27130/2018-02-08/'
          #aperture_photometry(path, 32, 'LF', 'coor', [2018, 1260], [20,25,30], [5,8,-0.02,0.02], 3, 0)

          
     if name=='photometry_class':
          from Photometry_Class import Photometry_Class
          path = '/home/nicholas/Data/IAC80/2015-01-06/corot1b/'
          XX = Photometry_Class(path, 'red', 3, 0)
          XX.aperture_photometry('star_coor', [1010, 981], [10,15,18], [5,15,-0.01,0.01], 2, 0)
          
if __name__ == '__main__': #---------------------------- Main function --------------------------#
     # FUNCTION CALL:
     test(name='photometry')    
