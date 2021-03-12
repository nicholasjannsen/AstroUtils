# AstroUtils

This is a personal libriary of self-written python modules for astronomy. 


## Photometry

Written April 2018 -- Nicholas Jannsen
Typeset in Python 2.7

Notice! 
---
At the moment only Photometry.py works as the class enviroment is under development. 

Usage:
---
In short this software performs aperture photometry for every star that have a image coordinate asigned. The first coordinate needs to belong to the target star for which a final light curve is desired. However, the flux and magnitude is calculated for every star, and all stars can be used as reference to correct for seeing variations during the total observation. Such a corrections is often refered to 'field (aperture) photometry' as the combined field of stars is used to correct for weather conditions during the night. The main function that needs to be called is 'aperture_photometry'.

Input paramters:
---
```
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
```
How to run the code:
---
Name of stellar coordinates. With this string the software look for a textfile where all the coordinates is within (x, y). The software `starfinder` within this software can be used to find and save all stellar coordinates. 

Aperture size of stellar aperture `r`, and minimum and maximum sky aperture radius `r_min` and `r_max`, respectively.

Name of the image data/fits files. The software is written such that if your images is named e.g. `image2018-XX-XXTXX-XX-XX.fits` you can simply use `image` (if no other data in the same folder is named exactly this).

The software only uses circular apertures and to select a good aperture size for the target star just set `plot=3` to plot the Signal-to-Noise Ratio (SNR) as a function of aperture radius. The best aperture is the usually the local maximum of this plot. This aperture is then also used for every reference star.
