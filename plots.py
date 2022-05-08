#!/usr/bin/env python3

import sys
import glob
import h5py
import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

import astropy.units as u
from astropy.coordinates import SkyCoord

from pylab import MaxNLocator

from colorama import Fore, Style

from scipy.ndimage import median_filter

from PIL import Image

#==============================================================#
#                    GLOBAL HARDCODE SETTINGS                  #
#==============================================================#

# Matplotlib.savefig resolution setting (default is 95)
dpi = 300

sp = 3
sms = 5
lw = 0.3
al = 0.5
dx = 150
cm = 'plasma'


#==============================================================#
#                           FUNCTIONS                          #
#==============================================================#


def plot(data, mark, xlab, ylab, title=None, subplot=0, legend=1,  axis=[1,1]):
    """
    General function to make fast plots in one command line:
    --------INPUT:
    data       (array)  : Data structure e.g. [data0, data1];  data0 and data1 have a x, y coloumn.
    mark       (list)   : If one have 2 datasets use e.g. ['b-', 'k.']
    xlab, ylab (string) : Labels on x and y
    title      (string) : Title
    legpos     (float)  : This can be 1, 2, 3, and 4 corresponding to each quadrant.
    subplot    (float)  : Different types of subplots.
    axis       (list)   : Procentage edge-space in x and y. E.g. [1, 5] to 1% in x and 5% in y. 
    """
    # Type of subplot:
    if subplot != 0: plot_subplot(subplot)

    # Plot data:
    if legend == 1:
        for i in range(len(data)):
            plt.plot(data[i][:,0], data[i][:,1], mark[i])
        plot_settings(xlab, ylab, title)
    if legend != 1:
        for i in range(len(data)):
            plt.plot(data[i][:,0], data[i][:,1], mark[i], label=legend[i+1])
        plot_settings(xlab, ylab, title, legend[0])

    # Axes setting:
    plot_axis(data[0][:,0], data[0][:,1], axis[0], axis[1])
    plt.show()






def SURF(x, y, z, xlab, ylab, zlab, title):
    # Find (x, y) value for maximum peak:
    z_max   = max(z)
    z_max_i = np.where(z==z_max)
    print('Best Period: {:.6f} days'.format(x[z_max_i[0][0]]))
    print('Best Phase : {:.6f} days'.format(y[z_max_i[1][0]]))
    # 3D plot:
    y, x = meshgrid(y, x)
    fig  = plt.figure()
    ax   = fig.gca(projection='3d')
    surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0)
    # Axes labels and title:
    ax.set_xlabel(xlab, fontsize=12); ax.tick_params(axis='x', labelsize=10)
    ax.set_ylabel(ylab, fontsize=12); ax.tick_params(axis='y', labelsize=10)
    ax.set_zlabel(zlab, fontsize=12); ax.tick_params(axis='z', labelsize=10)
    plt.title(title,    fontsize=15)
    # Extra settings:
    # ax.invert_xaxis()                          # Invert x-axis
    # ax.view_init(30, 45)                       # Viewing angle 
    # fig.colorbar(surf, shrink=0.5, aspect=8)   # Colorbar
    plt.show()









def plotTransitLightcurve(data, axes, zoom, xlabel, ylabel): # TODO do not work yet
    """
    This module make a nice plot of a exoplanet phase curve
    """
    # Subplot first:
    g  = gs.GridSpec(4, 1)
    ax = plt.subplot(g[0:2])
    # Plot data:
    plt.plot(data[:,0], data[:,1], 'r-')
    plot_settings(xlabel, ylabel)
    plot_axis(data[:,0], data[:,1], 1, 2)
    # Fancy subplot:
    #bbox_to_anchor=(380, 345)) # customized position 
    axins = zoomed_inset_axes(ax, zoom, loc=3)
    axins.plot(data[:,0], data[:,1], 'r-')
    x1, x2, y1, y2 = axes[0], axes[1], axes[2], axes[3] # specify the limits
    axins.set_xlim(x1, x2); plt.xticks(visible=False)
    axins.set_ylim(y1, y2); plt.yticks(visible=False)
    mark_inset(ax, axins, loc1=2, loc2=1, fc="none", ec="0.5")
    plt.show()





def plotPhasefoldLightcurve(data, P, phi, dT, depth, mark='o'):
    # Time, P, phi, dT is given in [days]
    t     = data[:,0]
    S     = data[:,1]
    phase = mod(t, P)
    # Whole phase plot:
    plt.plot(phase, S, mark)
    xpos = min(phase)+(max(phase)-min(phase))*0.70
    ypos = max(S)-(max(S)-min(S))*0.98
    plt.text(xpos, ypos,'$P={:.4f}$ days'.format(P), fontsize=18)
    plt.title('Phase folded', fontsize=18)
    plot_axes2([min(phase), max(phase)], S)
    plot_settings('$t$ $mod$ $P$ [days]', 'Signal')
    plt.show()
    # Zoom-in on transit:
    plt.plot(phase, S, mark)
    plot_axes2([phi-dT, phi+dT], [depth, max(S)])
    xpos  = phi-dT+((phi+dT)-(phi-dT))*0.70
    ypos1 = depth+0.0000
    ypos2 = depth+0.0002 
    plt.text(xpos, ypos2, '$P={:.4f}$ days'.format(P), fontsize=18)
    plt.text(xpos, ypos1, '$\phi={:.4f}$ days'.format(phi), fontsize=18)
    plt.title('Phase  folded', fontsize=18)
    plot_settings('$t$ $mod$ $P$ [days]', 'Signal')
    plt.show()





def drawStarsInSkyAitoff(fig, raStars, decStars, magStars, skymap=None, cbarOrientation=None):
    """
    Project and plot a catalog of stars on the sky in a Aitoff Galactic projection.
    This plot uses the astropy library to make the ICRS to Galactic coordinate
    transformation and the file gaiaDR3.png to set as background reference image.

    INPUT:
    @param array raStars:           Right ascension of stars    [deg]
    @param array decStars:          Declination of stars        [deg]
    @param array magStars:          Magnitudes of stars
    @param str   cbarOrientation:   Colorbar orientation. Default 'horizontal' else 'vertical'

    RETURN:
    None
    """
    # Convert coordinates from ICRS to Galactic using astropy

    gal = SkyCoord(raStars, decStars, frame='icrs', unit=u.deg)
    gal = gal.galactic

    # Plot Aitoff projection in Galactic coordinates

    plt.title('Aitoff projection in Galactic coordinates', fontsize=18, y=1.02)
    fig, ax = fig
    fs = 16
    if len(raStars) <= 1e2: ms = 3
    if len(raStars) >= 1e2 and len(raStars) < 1e3: ms = 1
    if len(raStars) >= 1e3 and len(raStars) < 1e4: ms = 0.5
    if len(raStars) >= 1e4: ms = 0.1

    # Plot Galactic map as background (e.g. Gaia DR3)
    # E.g.: skymap = plt.imread('skymap.png')

    if skymap is not None:
        ax.imshow(skymap)

    # Add the sky projection ontop as transparent layer

    axes = fig.add_subplot(111, projection='aitoff', facecolor='none')

    # Plot the targets on the sky

    im = plt.scatter(-gal.l.wrap_at('180d').radian, gal.b.radian, c=magStars, s=ms, cmap='autumn_r', zorder=3)

    # Vertical or horizontal colorbar showing magnitudes

    if cbarOrientation == 'vertical':
        cbarax = fig.add_axes([0.905, 0.2, 0.02, 0.57])
        cbar = plt.colorbar(im, orientation='vertical', cax=cbarax, extend='both')
        cbar.set_label('Bessel V Magnitude', fontsize=fs)
        cbar.ax.tick_params(labelsize=fs)
    else:
        cbarax = fig.add_axes([0.25, 0.06, 0.525, 0.03])
        cbar = plt.colorbar(im, orientation='horizontal', cax=cbarax, extend='both')
        cbar.set_label('Gaia V Magnitude', fontsize=fs)
        cbar.ax.tick_params(labelsize=fs)

    # Change the tick labels so that they are 0->360, rather than -180->+180

    tickLabels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
    tickLabels = np.remainder(tickLabels+360, 360)
    axes.set_xticklabels(tickLabels)

    # Change y ticks and remove last to make space for title

    tickLabels = np.array([-75, -60, -45, -30, -15, 0, 15, 30, 45, 60, ''])
    axes.set_yticklabels(tickLabels)

    # Change color of x tick labels

    axes.tick_params(axis='x', colors='w')

    # Increase x and y tick labels

    axes.xaxis.set_tick_params(labelsize=fs+1)
    axes.yaxis.set_tick_params(labelsize=fs)

    # Add axis labels

    axes.set_xlabel('Longitude [deg]', fontsize=fs)
    axes.set_ylabel('Latitude [deg]', fontsize=fs)

    # Set grid and remore outer ticks (if set by default)

    axes.grid(True, alpha=0.3)
    ax.axis('off')

    # TODO plot PLATO FOV as alpha channel

    # Convert Galactic coordinates to RA/Dec (Internation Celestial Reference System)
    #galPointing = SkyCoord(pointingField[0], pointingField[1], unit=u.deg, frame='galactic')
    #equPointing = galPointing.icrs

    # tilt = np.deg2rad(9.2)
    # raPlatform  = np.deg2rad(-85)
    # decPlatform = np.deg2rad(-50)
    # rFOV = np.deg2rad(18)

    # coords = np.array([[raPlatform-tilt, decPlatform+tilt], [raPlatform-tilt, decPlatform-tilt],
    #                    [raPlatform+tilt, decPlatform-tilt], [raPlatform+tilt, decPlatform+tilt]])

    # px1, py1 = circle(coords[0,0], coords[0,1], rFOV)
    # px2, py2 = circle(coords[1,0], coords[1,1], rFOV)
    # px3, py3 = circle(coords[2,0], coords[2,1], rFOV)
    # px4, py4 = circle(coords[3,0], coords[3,1], rFOV)

    # Plot PLATO goup-FOV circles
    # g1 = plt.Circle(coords[0], rFOV, fc='c', fill=True, alpha=0.2, zorder=1)
    # g2 = plt.Circle(coords[1], rFOV, fc='c', fill=True, alpha=0.2, zorder=1)
    # g3 = plt.Circle(coords[2], rFOV, fc='c', fill=True, alpha=0.2, zorder=1)
    # g4 = plt.Circle(coords[3], rFOV, fc='c', fill=True, alpha=0.2, zorder=1)
    # axes.add_artist(g1)
    # axes.add_artist(g2)
    # axes.add_artist(g3)
    # axes.add_artist(g4)
    # axes.plot(px1, py1, 'c', lw=0.7, zorder=2)
    # axes.plot(px2, py2, 'c', lw=0.7, zorder=2)
    # axes.plot(px3, py3, 'c', lw=0.7, zorder=2)
    # axes.plot(px4, py4, 'c', lw=0.7, zorder=2)

    # TODO plot a cartesian view of all PIC targets visible by different CCDs
    # coords = SkyCoord(ra, dec, unit=u.deg, frame='icrs')
    # from astropy.coordinates import EarthLocation
    # location = EarthLocation.of_site('greenwich')
    # from astropy.time import Time
    # from astropy.coordinates import AltAz
    # local_frame = AltAz(location=location, obstime=Time.now())
    # local_coords = coords.transform_to(local_frame)
    # coords_catesian = local_coords.cartesian
    # plt.figure()
    # plt.plot(ra, dec, 'b.', alpha=0.3)
    # plt.show()





def plotStellarSampleDistributions(fig, mag, magCon, magRange, numConPerTar, distCon):

    fig1, ax = fig

    magbinTar  = 0.1
    if magRange[1]-magRange[0] > 10: magbinTar = 0.2
    binsizeTar = int((magRange[1] - magRange[0]) / magbinTar) + 1
    binlistTar = np.linspace(magRange[0], magRange[1], binsizeTar)

    ax[0,0].hist(mag, binlistTar, facecolor='b', edgecolor='b', fill=True, alpha=0.3)
    ax[0,0].set_title('Magnitude distribution of PIC targets')
    ax[0,0].set_xlabel('Gaia V Magnitude')
    ax[0,0].set_ylabel('Number of stars')
    ax[0,0].locator_params(axis='y', integer=True)
    ax[0,0].tick_params(axis='x', which='minor', bottom=True, top=False)
    ax[0,0].tick_params(axis='x', which='major', bottom=True, top=False)
    ax[0,0].tick_params(axis='y', which='minor', left=False, right=False)
    ax[0,0].tick_params(axis='y', which='major', left=True, right=False)
    ax[0,0].grid(axis='y', color='gray', alpha=0.3)

    magbinCon  = 0.2
    binsizeCon = int((np.max(magCon) - np.min(magCon)) / magbinCon) + 1
    binlistCon = np.linspace(round(np.min(magCon)), round(np.max(magCon)), binsizeCon)

    ax[0,1].hist(magCon, binlistCon, facecolor='m', edgecolor='m', fill=True, alpha=0.3)
    ax[0,1].set_title('Magnitude distribution of PIC contaminants')
    ax[0,1].set_xlabel('Gaia V Magnitude')
    ax[0,1].set_ylabel('Number of stars')
    ax[0,1].tick_params(axis='x', which='minor', bottom=True, top=False)
    ax[0,1].tick_params(axis='x', which='major', bottom=True, top=False)
    ax[0,1].tick_params(axis='y', which='minor', left=False, right=False)
    ax[0,1].tick_params(axis='y', which='major', left=True, right=False)
    ax[0,1].grid(axis='y', color='gray', alpha=0.3)

    numbinCon  = 1
    binsizeNum = int((np.max(numConPerTar) - 0) / numbinCon) + 2
    binlistNum = np.linspace(-0.5, np.max(numConPerTar)+0.5, binsizeNum)  # -0.5 because num x-axis

    ax[1,0].hist(numConPerTar, binlistNum, facecolor='g', edgecolor='g', fill=True, log=True, alpha=0.3)
    ax[1,0].set_title('Number distribution of contaminants per target')
    ax[1,0].set_xlabel('Number of contaminants')
    ax[1,0].set_ylabel('Number of targets')
    ax[1,0].tick_params(axis='x', which='minor', bottom=False, top=False)
    ax[1,0].tick_params(axis='x', which='major', bottom=True, top=False)
    ax[1,0].tick_params(axis='y', which='minor', left=False, right=False)
    ax[1,0].tick_params(axis='y', which='major', left=True, right=False)
    # ax[1,0].xaxis.set_major_locator(MaxNLocator(integer=True))
    # ax[1,0].locator_params(axis='both', integer=True)
    ax[1,0].grid(axis='y', color='gray', alpha=0.3)

    distbinCon  = 1.0
    binsizeDist = int((np.max(distCon) - np.min(distCon)) / distbinCon) + 2  # +1 extra because zero is rare
    binlistDist = np.linspace(round(np.min(distCon)), round(np.max(distCon)), binsizeDist)

    ax[1,1].hist(distCon, binlistDist, facecolor='orange', edgecolor='orange', fill=True, alpha=0.4)
    ax[1,1].set_title('Distance distribution of contaminants')
    ax[1,1].set_xlabel('Distances [arcsec]')
    ax[1,1].set_ylabel('Number of stars')
    ax[1,1].locator_params(axis='y', integer=True)
    ax[1,1].tick_params(axis='x', which='minor', bottom=True, top=False)
    ax[1,1].tick_params(axis='x', which='major', bottom=True, top=False)
    ax[1,1].tick_params(axis='y', which='minor', left=False, right=False)
    ax[1,1].tick_params(axis='y', which='major', left=True, right=False)
    ax[1,1].grid(axis='y', color='gray', alpha=0.3)

    plt.tight_layout()






def plotTimeserie(time, data, title, labels, rms=False, ylims=False, save=False):
    """
    TIMESERIES

    INPUT:
    @param time     : Array of times.
    @param datasets : List or array of individual signal arrays.
    @param title    : Title in a string.
    @param labels   : String of labels where the first is the xlabel and the rest is ylabels.
                      Last entry is the physical unit of the signal.
    @param rms      : Array of Root-Mean-Square (RMS) values for each input signal.

    OUTPUT:
    Plot or/and saved plot to PNG.
    """
    # Hardcode parameters
    colors = ['royalblue', 'lightseagreen', 'limegreen', 'y', 'r', 'k']

    # Datasets to loop over
    numData = len(data)

    # Handle yaxis limits
    if ylims is False:
        lim = 1.2*np.max(np.abs(data))

    # Adjust linewidth after data
    if len(time) < 1e3: lw = 1
    else: lw = 0.5

    # Make plot
    fig, ax = plt.subplots(numData, 1, figsize=(5, 2*numData))

    for plot in range(numData):

        # Plot timeseries
        ax[plot].plot(time, data[plot], '-', c=colors[plot], lw=lw)

        # Add RMS lines
        if rms is not False:
            ax[plot].axhline(+rms[plot], c='k', ls='--', lw=0.7, label='RMS = {0:.3f} {1}'.format(rms[plot], labels[-1]))
            ax[plot].axhline(-rms[plot], c='k', ls='--', lw=0.7)
            ax[plot].legend(loc='upper right')

        # Latter settings
        ax[plot].set_ylabel('{0} [{1}]'.format(labels[plot+1], labels[-1]))
        ax[plot].set_xlim(np.min(time), np.max(time))
        ax[plot].set_ylim(-lim, +lim)

        # Remove tick labels on x axis except for last plot
        if plot < numData-1:
            ax[plot].tick_params(labelbottom=False)
            ax[plot].set_ylim(-lim, +lim)

    # Remaining
    ax[0].set_title(title)
    ax[-1].set_xlabel(labels[0])
    plt.tight_layout()
    plt.subplots_adjust(hspace = .001)
    plt.show()

    # Save if string prefix is given
    if save is not False:
        outputDir  = save[0]
        prefixName = save[1]
        filename = outputDir + '/plot{}Timeseries.png'.format(prefixName)
        fig.savefig(filename, bbox_inches='tight', dpi=dpi)







def plotPSD(time, datasets, title, labels, freqlim=False, save=False):
    """
    POWER DENSITY SPECTRE

    INPUT:
    @param time : array of times
    @param datasets : list or array of individual signal arrays
    @param title : title in a string
    @param labels : String of labels where the first is the xlabel and the rest is ylabels
    @param rms : array of Root-Mean_square values for each signal variation

    OUTPUT:
    Plot or/and saved plot to PNG.
    """
    # Hardcode parameters
    colors = ['tomato', 'darkorange', 'gold', 'k', 'k', 'k']

    # Datasets to loop over
    numData = len(datasets)

    # Adjust linewidth after data
    if len(time) < 1e3: lw = 1
    else: lw = 0.5

    # Prepare Median convolve data (control boundaries)
    filt = int(len(time) / 10.)
    if filt > 10000: filt = 3000

    # Compute frequencies uptil the Nyquist frequency
    frequencies = np.fft.fftfreq(len(time), d=np.diff(time)[0])

    # Make plot
    fig, ax = plt.subplots(numData, 1, figsize=(5, 2*numData))

    # Special care for the ylabel
    if (numData % 2) == 0:
        fig.text(0.0, 0.5, r'Amplitude [ppm$^2$ $\mu$Hz$^{-1}$]', va='center', rotation='vertical')
        plt.tight_layout(pad=3)
    else:
        ax[int(numData/2)].set_ylabel(r'Amplitude [arcsec$^2$ $\mu$Hz$^{-1}$]')
        plt.tight_layout(pad=2)

    # Common plot items to loop over
    for plot in range(numData):

        # Require even number of data points
        if (len(time) % 2) != 0:
            dx = len(time) - 1
        else:
            dx = len(time)
        time = time[:dx]
        data = datasets[:dx][plot]

        # Start with wavelength varying signal; units are DN/s, values are real
        sp_model    = np.fft.fft(data)
        power_model = np.sqrt(sp_model.real**2 + sp_model.imag**2)
        sp_med      = median_filter(power_model, filt)

        # Sort away the negative reflection image
        positives = np.where(frequencies > 0)
        frequencies = frequencies[positives]
        power_model = power_model[positives]
        sp_med      = sp_med[positives]

        # Plot PSD
        ax[plot].plot(frequencies, power_model, '-', c=colors[plot], lw=lw, label=labels[plot])
        ax[plot].plot(frequencies, sp_med, 'k-', lw=lw+0.5, label='Median filter')

        # Log scaling
        ax[plot].set_xscale("log")
        ax[plot].set_yscale("log")
        ax[plot].legend(loc='lower left')

        # Remove tick labels on x axis except for last plot
        if plot < numData-1: ax[plot].tick_params(labelbottom=False)

        # xlimits
        if freqlim is not False:
            ax[plot].set_xlim(freqlim, np.max(frequencies))
        else:
            ax[plot].set_xlim(np.min(frequencies), np.max(frequencies))

    # Special settings when plotting from Prime
    ax[0].set_ylim(1.5e-4, 1.5e3)
    ax[1].set_ylim(1.5e-4, 1.5e3)
    ax[2].set_ylim(1.5e-4, 1.5e3)
    ax[0].set_title(title)

    # Perform plot after loop
    ax[-1].set_xlabel(r'Frequency [$\mu$Hz]')
    plt.tight_layout()
    plt.subplots_adjust(hspace = .001)
    plt.show()

    # Save plot if specified
    if save is not False:
        outputDir  = save[0]
        prefixName = save[1]
        filename = outputDir + '/plot{}PSD.png'.format(prefixName)
        fig.savefig(filename, bbox_inches='tight', dpi=dpi)







def plotYawPitchRoll(time, data, rms, clabel, save=False, longrun=False):
    """
    PLOT PARAMETER CORRELATIONS

    INPUT:
    @param time    : Array of times.
    @param data    : List or array of individual data arrays.
    @param rms     : List or array of float Root-Mean-square (RMS) values for each time varying dataset.
    @param title   : String of the main title.
    @param clabel  : String with the color-bar time-label.
    @param longrun : Flag to plot the correlations for the entire input timeseries if True.

    OUTPUT:
    Plot or/and saved plot to PNG.
    """
    # Hardcode values
    labels = ['Yaw [arcsec]', 'Pitch [arcsec]', 'Roll [arcsec]']
 
    # A single plot

    if len(data) == 2:
        pass
        # fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        # im = ax.scatter(data[0], data[1], c=time, s=3, cmap='plasma')
        # # Adjust and customize colorbars
        # div = make_axes_locatable(ax)
        # cax = div.append_axes("right", size="5%", pad=0.015)
        # cbar = plt.colorbar(im, ax=ax, cax=cax, extend='max')
        # cbar.set_label(labels[0])
        # # First set grid
        # ax.grid(c='gray', ls='-', lw=0.5, alpha=0.5)
        # # Labels
        # ax.set_xlabel(labels[1])
        # ax.set_ylabel(labels[2])
        # # Limits
        # lim = 4*np.mean(rms)
        # ax.set_xlim(-lim, lim)
        # ax.set_ylim(-lim, lim)
        # # Secure equal axes
        # ax.set_aspect('equal', 'box')

    # PLOT CORRELATIONS ON SHORT TIME SCALES

    else:
        # Limits and grid
        lim = 0.125  #np.max(np.abs(data))
        nticks = 5

        time = time - np.min(time)
        fig, ax = plt.subplots(sp, sp, figsize=(10, 7))

        # ax[0, 0].set_title(titles[0])
        # ax[0, 1].set_title(titles[1])
        # ax[0, 2].set_title(titles[2])

        for row in range(sp):

            # Plots
            ax[row, 0].plot(data[1][dx*row:dx*(row+1)], data[0][dx*row:dx*(row+1)], 'k-', alpha=al, lw=lw, zorder=1)
            ax[row, 1].plot(data[2][dx*row:dx*(row+1)], data[0][dx*row:dx*(row+1)], 'k-', alpha=al, lw=lw, zorder=1)
            ax[row, 2].plot(data[2][dx*row:dx*(row+1)], data[1][dx*row:dx*(row+1)], 'k-', alpha=al, lw=lw, zorder=1)
            im0 = ax[row, 0].scatter(data[1][dx*row:dx*(row+1)], data[0][dx*row:dx*(row+1)], c=time[dx*row:dx*(row+1)], s=sms, cmap=cm, zorder=2)
            im1 = ax[row, 1].scatter(data[2][dx*row:dx*(row+1)], data[0][dx*row:dx*(row+1)], c=time[dx*row:dx*(row+1)], s=sms, cmap=cm, zorder=2)
            im2 = ax[row, 2].scatter(data[2][dx*row:dx*(row+1)], data[1][dx*row:dx*(row+1)], c=time[dx*row:dx*(row+1)], s=sms, cmap=cm, zorder=2)

            # Labels
            ax[2, 0].set_xlabel(labels[1])
            ax[2, 1].set_xlabel(labels[2])
            ax[2, 2].set_xlabel(labels[2])
            ax[row, 0].set_ylabel(labels[0])
            ax[row, 1].set_ylabel(labels[0])
            ax[row, 2].set_ylabel(labels[1])

            # Remove tick labels on x axis except for last plot
            if row < 2:
                ax[row, 0].tick_params(labelbottom=False)
                ax[row, 1].tick_params(labelbottom=False)
                ax[row, 2].tick_params(labelbottom=False)

            # Duplicate settings for each plotted row
            for col, im in zip(range(sp), [im0, im1, im2]):

                # Axes limits
                ax[row, col].set_xlim(-lim, lim)
                ax[row, col].set_ylim(-lim, lim)
                ax[row, col].set_aspect('equal', 'box')

                # Force the same number of ticks
                ax[row, col].xaxis.set_major_locator(MaxNLocator(nticks))
                ax[row, col].yaxis.set_major_locator(MaxNLocator(nticks))

                # Set grid
                ax[row, col].grid(c='gray', ls='-', lw=lw, alpha=al)

                # Color bars
                div = make_axes_locatable(ax[row, col])
                cax = div.append_axes('right', size='10%', pad=0.1)
                cbar = plt.colorbar(im, ax=ax[row, col], cax=cax, extend='max')
                cbar.ax.invert_yaxis()
                if im == im2:
                    cbar.set_label(clabel)
                else:
                    cbar.remove()

        # Adjust subplot spacing
        fig.tight_layout()
        fig.subplots_adjust(hspace = .001)
        plt.show()


    # PLOT CORRELATIONS FOR ENTIRE TIMESERIES

    if longrun is True:

        # Limits and grid
        lim = 0.28  #np.max(np.abs(data))
        nticks = 6
        time = time/(60*60)

        fig1, ax1 = plt.subplots(1, 3, figsize=(10, 2.8))

        # Plot
        ax1[0].plot(data[1], data[0], 'k-', alpha=al, lw=lw, zorder=1)
        ax1[1].plot(data[2], data[0], 'k-', alpha=al, lw=lw, zorder=1)
        ax1[2].plot(data[2], data[1], 'k-', alpha=al, lw=lw, zorder=1)
        im0 = ax1[0].scatter(data[1], data[0], c=time, s=2, cmap='magma', zorder=2)
        im1 = ax1[1].scatter(data[2], data[0], c=time, s=2, cmap='magma', zorder=2)
        im2 = ax1[2].scatter(data[2], data[1], c=time, s=2, cmap='magma', zorder=2)

        # Labels
        ax1[0].set_xlabel(labels[1])
        ax1[0].set_ylabel(labels[0])
        ax1[1].set_xlabel(labels[2])
        ax1[1].set_ylabel(labels[0])
        ax1[2].set_xlabel(labels[2])
        ax1[2].set_ylabel(labels[1])

        # Duplicate settings
        for plot in range(3):

            # Adjust axes
            ax1[plot].set_xlim(-lim, lim)
            ax1[plot].set_ylim(-lim, lim)
            ax1[plot].set_aspect('equal', 'box')

            # Force the same number of ticks
            ax1[plot].xaxis.set_major_locator(MaxNLocator(nticks))
            ax1[plot].yaxis.set_major_locator(MaxNLocator(nticks))

            # Plot grid
            ax1[plot].set_axisbelow(False)
            ax1[plot].grid(c='gray', ls='-', lw=lw, alpha=al)

            # We put the colorbar for then remove
            # This is need to keep the same size of each subplot..
            div1 = make_axes_locatable(ax1[2])
            cax1 = div1.append_axes("right", size="10%", pad=0.1)
            cbar1 = plt.colorbar(im2, ax=ax1[2], cax=cax1, extend='max')
            cbar1.ax.invert_yaxis()
            if plot == 2:
                cbar1.set_label('Time [hours]')
            else:
                cbar1.remove()

        # Plot this figure separately
        fig1.tight_layout()
        fig1.subplots_adjust(hspace = .001)
        plt.show()

    # Save if specified
    if save is not False:
        outputDir  = save[0]
        prefixName = save[1]
        filename = outputDir + '/plot{}Correlations.png'.format(prefixName)
        fig.savefig(filename, bbox_inches='tight', dpi=dpi)

        # Save the long timeseries
        if longrun is True:
            filename = outputDir + '/plot{}CorrelationsLong.png'.format(prefixName)
            fig1.savefig(filename, bbox_inches='tight', dpi=dpi)





def skycoor(cat0, cat1):
    plt.plot(cat0[:,0], cat0[:,1], 'wo', markersize=7, label='$Cat_0$', markeredgecolor='k')
    plt.plot(cat1[:,0], cat1[:,1], 'bo', markersize=3, label='$Cat_1$', markeredgecolor='none')
    plt.legend(bbox_to_anchor=(0.25, 1.02), numpoints=1, frameon=False, fontsize=15)
    # Other settings:
    plot_settings(r'$\alpha$ [deg]', r'$\delta$ [deg]')
    # plt.ticklabel_format(useOffset=False)
    plt.axis('equal')
    #plt.ylim(11.59, 12.01)
    plt.minorticks_on()
    # plt.legend(bbox_to_anchor=(1.01, 1.02),numpoints=1,frameon=False,fontsize=FS)
    # plt.xticks(np.arange(132.65,133.05,0.1))
    # plt.legend(loc='upper right',numpoints=1,frameon=False,fontsize=FS)
    plt.show()





def linear(inputArray, scale_min=None, scale_max=None):
    """
    Performs linear scaling of the input np array.

	@type inputArray: np array
	@param inputArray: image data array
	@type scale_min: float
	@param scale_min: minimum data value
	@type scale_max: float
	@param scale_max: maximum data value
	@rtype: np array
	@return: image data array
	"""
    imageData = np.array(inputArray, copy=True)

    if scale_min == None:
        scale_min = imageData.min()
    if scale_max == None:
        scale_max = imageData.max()

    imageData = imageData.clip(min=scale_min, max=scale_max)
    imageData = (imageData - scale_min) / (scale_max - scale_min)
    indices = np.where(imageData < 0)
    imageData[indices] = 0.0
    indices = np.where(imageData > 1)
    imageData[indices] = 1.0

    return imageData


# def sqrt(inputArray, scale_min=None, scale_max=None):
#     """
#     Performs sqrt scaling of the input np array.

# 	@type inputArray: np array
# 	@param inputArray: image data array
# 	@type scale_min: float
# 	@param scale_min: minimum data value
# 	@type scale_max: float
# 	@param scale_max: maximum data value
# 	@rtype: np array
# 	@return: image data array
# 	"""

# 	imageData = np.array(inputArray, copy=True)

# 	if scale_min == None:
# 		scale_min = imageData.min()
# 	if scale_max == None:
# 		scale_max = imageData.max()

# 	imageData = imageData.clip(min=scale_min, max=scale_max)
# 	imageData = imageData - scale_min
# 	indices = np.where(imageData < 0)
# 	imageData[indices] = 0.0
# 	imageData = np.sqrt(imageData)
# 	imageData = imageData / math.sqrt(scale_max - scale_min)

# 	return imageData


# def log(inputArray, scale_min=None, scale_max=None):
#     """
#     Performs log10 scaling of the input np array.

# 	@type inputArray: np array
# 	@param inputArray: image data array
# 	@type scale_min: float
# 	@param scale_min: minimum data value
# 	@type scale_max: float
# 	@param scale_max: maximum data value
# 	@rtype: np array
# 	@return: image data array
# 	"""
# 	imageData=np.array(inputArray, copy=True)

# 	if scale_min == None:
# 		scale_min = imageData.min()

# 	if scale_max == None:
# 		scale_max = imageData.max()

#     factor = math.log10(scale_max - scale_min)
# 	indices0 = np.where(imageData < scale_min)
# 	indices1 = np.where((imageData >= scale_min) & (imageData <= scale_max))
# 	indices2 = np.where(imageData > scale_max)
# 	imageData[indices0] = 0.0
# 	imageData[indices2] = 1.0
# 	try:
#         imageData[indices1] = np.log10(imageData[indices1])/factor
# 	except:
#         print "Error on math.log10 for ", (imageData[i][j] - scale_min)

#     return imageData


# def asinh(inputArray, scale_min=None, scale_max=None, non_linear=2.0):
# 	"""
#     Performs asinh scaling of the input np array.

# 	@type inputArray: np array
# 	@param inputArray: image data array
# 	@type scale_min: float
# 	@param scale_min: minimum data value
# 	@type scale_max: float
# 	@param scale_max: maximum data value
# 	@type non_linear: float
# 	@param non_linear: non-linearity factor
# 	@rtype: np array
# 	@return: image data array
# 	"""

# 	imageData=np.array(inputArray, copy=True)

# 	if scale_min == None:
# 		scale_min = imageData.min()
# 	if scale_max == None:
# 		scale_max = imageData.max()

# 	factor = np.arcsinh((scale_max - scale_min)/non_linear)
# 	indices0 = np.where(imageData < scale_min)
# 	indices1 = np.where((imageData >= scale_min) & (imageData <= scale_max))
# 	indices2 = np.where(imageData > scale_max)
# 	imageData[indices0] = 0.0
# 	imageData[indices2] = 1.0
# 	imageData[indices1] = np.arcsinh((imageData[indices1] - \
# 	scale_min)/non_linear)/factor

# 	return imageData


def FITS(img, scale='linear', ds=2, cmap='gray', colorbar=True, rgb=False):

    if scale == 'linear': scale = linear
    if scale == 'sqrt'  : scale = sqrt
    if scale == 'log'   : scale = log
    if scale == 'asinh' : scale = asinh

    # if rgb is True:
    #     pylab.clf()
    #     pylab.imshow(img, aspect='equal')
    # else:
    img_min, img_max = img.mean()-ds*img.std(), img.mean()+ds*img.std()

    plt.imshow(scale(img, img_min, img_max), cmap=cmap, origin='lower')

    if (colorbar is True) and (rgb is False):
        # Prepare for scaled colorbar:
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        ax = plt.gca()
        im = ax.imshow(scale(img, img_min, img_max), cmap=cmap, origin='lower')

        # create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        divider = make_axes_locatable(ax)
        cax  = divider.append_axes('right', size='5%', pad=0.05)  # right, left, top, bottom

        # Make colorbar and append ylabel and axis labels:
        cbar = plt.colorbar(im, cax=cax)# orientation='horizontal')
        cbar.ax.set_ylabel('Normalized Counts')






def gif():

    # Initialize input parameters:
    argv = sys.argv      # Arguments parsed to pipeline
    argc = len(argv)     # Number of arguments parsed to pipeline

    # Help usage function:
    def help():
        print(Fore.BLUE+Style.BRIGHT+"""Usage : %s <gif-filename> <path/to/images>"""
              % argv[0][2:]+Style.RESET_ALL)

    # Print usage if no arguments are given:
    if argc==1: help(); sys.exit()

    # Parsing the path:
    files = argv[2:]

    # filepaths
    fp_in  = "/path/to/image_*.png"
    fp_out = "{}.gif".fotmat(argv[1])

    # https://pillow.readthedocs.io/en/stable/handbook/image-file-formats.html#gif
    img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]
    img.save(fp=fp_out, format='GIF', append_images=imgs, save_all=True, duration=200, loop=0)
