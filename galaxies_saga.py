import numpy as np
import matplotlib.pyplot as plt
from Plot_Tools import plot_settings
import uncertainties.unumpy as up
from uncertainties import ufloat
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
import sys

# saga1:
data1 = np.loadtxt('/home/nicholas/Dropbox/Galaxies/data/saga1.dat')
KIC  = data1[:,0]
FeH  = up.uarray(data1[:,5],  data1[:,6])
T    = up.uarray(data1[:,7],  data1[:,8])
R    = up.uarray(data1[:,9], data1[:,10])
M    = up.uarray(data1[:,12], data1[:,13])
logg = up.uarray(data1[:,15], data1[:,16])
rho  = up.uarray(data1[:,18], data1[:,19])
d    = up.uarray(data1[:,21], data1[:,22])
star   = data1[:,24]
FeH_qf = data1[:,25]
pho_qf = data1[:,26]
by     = data1[:,27]
m1     = data1[:,28]
c1     = data1[:,29]


############## Make error definition:

def val(x): return up.nominal_values(x)
def err(x): return up.std_devs(x)

###############

#print data1[0,6]; sys.exit()
# saga2:
data2 = np.loadtxt('/home/nicholas/Dropbox/Galaxies/data/saga2.dat')
age  = up.uarray(data2[:,1], data2[:,2])/1e3     # Gyr

# Extra:
data3 = np.loadtxt('/home/nicholas/Dropbox/Galaxies/data/saga3.dat')
y    = up.uarray(data3[:,6], data3[:,7])
E_BV = data3[:,8]

##############################################################################################################


# Correct for extinction:
y0 = y - 0.75*E_BV
# Absolute y magnitude:
My = y0 - 5*up.log10(d) + 5  

#------------------------
# SORT AFTER COLOR RANGE:
#------------------------

# Color selected:
i0 = (0.6<by)*(by<0.8)*(y0<13.5)   #, np.where(by<0.8)[0]

# y-magitude vs. reddening:

# fig, ax = plt.subplots(1,1)
# xlab, ylab = '$b-y$', '$y_0$'
# ax.errorbar(by,     val(y0),     err(y0),     0, 'k.', label='Raw',   alpha=.2)
# ax.errorbar(by[i0], val(y0[i0]), err(y0[i0]), 0, 'r.', label='Range', alpha=.3)
# plot_settings(fig, ax, xlab, ylab, 3)
# plt.gca().invert_yaxis()
# plt.show()

#sys.exit()

#--------------------------
# SORT AFTER UNCERTAINTIES:
#--------------------------

# Stars with an procentage wise uncertainty:
hist_M  = err(M) /val(M) * 100 
hist_y0 = err(y0)/val(y0)* 100
hist_d  = err(d) /val(d) * 100
n = 100

# fig, ax = plt.subplots(1,1)
# xlab, ylab = 'Procentage (\%)', '$N$'
# plt.hist(hist_M,  bins=n, histtype='step', color='k', label='$M$')
# plt.hist(hist_y0, bins=n, histtype='step', color='b', label='$y_0$')
# plt.hist(hist_d,  bins=n, histtype='step', color='r', label='$d$')
# plot_settings(fig, ax, xlab, ylab, 1)
# plt.show()

# From the sigma-mass plots we harcode the thresholdes:
n0, n1, n2  = 7.7, 0.15, 4
i1 = (hist_M<n0)*(hist_y0<n1)*(hist_d<n2)
i2 = i1*i0
#--------------------------
# IMF: MAGNITUDE VS. MASSS:
#--------------------------

# Plot in 3 seperate plots:

# fig, ax = plt.subplots(3,1)#[plt.subplots(1,3, ) for i in range(3)]
# xlab, ylab = '$M_y$', 'Mass ($M_{\odot}$)'
# MM, MM_err = up.nominal_values(M),  up.std_devs(M)
# MY, MY_err = up.nominal_values(My), up.std_devs(My) 
# #ax = axs[0]
# ax[0].errorbar(MY,     MM,     MM_err,     MY_err,     'k.', label='Raw',         alpha=.1)
# ax[1].errorbar(MY[i0], MM[i0], MM_err[i0], MY_err[i0], 'b.', label='Color',       alpha=.1)
# ax[2].errorbar(MY[i1], MM[i1], MM_err[i1], MY_err[i1], 'r.', label='Uncertainty', alpha=.5)
# plot_settings(fig, ax[0], xlab)
# plot_settings(fig, ax[1], xlab, ylab, 1)
# plot_settings(fig, ax[2])
# plt.subplots_adjust(wspace=0, hspace=0)
# plt.show()

# fig, ax = plt.subplots(1,1)
# xlab, ylab = '$M_y$', 'Mass ($M_{\odot}$)'
MM, MM_err = val(M),  err(M)
MY, MY_err = val(My), err(My) 
# ax.errorbar(MY,     MM,     MM_err,     MY_err,     'k.', label='Raw',           markersize=3, alpha=.2)
# ax.errorbar(MY[i0], MM[i0], MM_err[i0], MY_err[i0], 'b.', label='Color',         markersize=3, alpha=.2)
# ax.errorbar(MY[i1], MM[i1], MM_err[i1], MY_err[i1], 'r.', label='Color + Uncertainty', markersize=3, alpha=.4)
# plot_settings(fig, ax, xlab, ylab, 1)
# plt.subplots_adjust(wspace=0, hspace=0)
# plt.show()



# Overplot line: 
# from scipy.interpolate import interp1d
# x = [2.00, 2.50, 3.00, 3.50, 4.00, 5.00, 6.00]
# y = [1.80, 1.60, 1.40, 1.25, 1.15, 0.95, 0.85]
# xnew = np.linspace(2, 6, 100)
# f_inter = interp1d(x, y, kind='cubic')

# fig, ax = plt.subplots(1,1)
# xlab, ylab =  '$M_y$', '$M$ ($M_{\odot}$)'
# ax.errorbar(MY[i1], MM[i1], MM_err[i1], MY_err[i1], 'r.', label='Saga stars', markersize=3, alpha=.4)
# plt.plot(xnew, f_inter(xnew), 'k-', label='$M_V$ theo-curve', alpha=.8)
# plot_settings(fig, ax, xlab, ylab, 1)
# plt.show()


#------------------------
# REDDENING VS. DISTANCE:
#------------------------

# fig, ax = plt.subplots(1,1)
# xlab, ylab = '$d$ (kpc)', '$E(B-V)$'
# ax.errorbar(val(d)/1e3, E_BV, 0, err(d)/1e3,  'k.', alpha=.2)
# plot_settings(fig, ax, xlab, ylab, 1)
# plt.show()

# This plot shows 3 distinct reddening branches equal to the bulge, the thin, and the tick disk:

# Let crude sortation of uncertainties:
n0, n1, n2  = 15, 0.15, 8
i4 = (hist_M<n0)*(hist_y0<n1)*(hist_d<n2)

# The 3 branches:
b0 = (E_BV>0.050)*(E_BV<0.093)*(d>1400)*i4
b1 = (E_BV>0.093)*(E_BV<0.120)*(d>1600)*i4
b2 = (E_BV>0.120)             *(d>1700)*i4        

# Clump of stars:
i5 = (E_BV>0.137)*(E_BV<0.141)*(d>2300)*(d<2650)

# distance vs. reddening:

# fig, ax = plt.subplots(1,1)
# xlab, ylab = '$d$ (kpc)', '$E(B-V)$'
# ax.errorbar(val(d)/1e3, E_BV, 0, up.std_devs(d)/1e3,  'k.', alpha=.1)
# plt.plot(val(d[b0])/1e3, E_BV[b0], 'b.', label='Branch 1', alpha=.2)
# plt.plot(val(d[b1])/1e3, E_BV[b1], 'r.', label='Branch 2', alpha=.2)
# plt.plot(val(d[b2])/1e3, E_BV[b2], 'g.', label='Branch 3', alpha=.2)
# #plt.plot(val(d[i5])/1e3, E_BV[i5], 'm.', label='Branch 3', alpha=.2)
# plot_settings(fig, ax, xlab, ylab, 4)
# plt.show()

# y-magitude vs. reddening:

# fig, ax = plt.subplots(1,1)
# xlab, ylab = '$y_0$', '$E(B-V)$'
# ax.errorbar(val(y0), E_BV, 0, err(y), 'k.', alpha=.2)
# #plt.plot(up.nominal_values(y[b0]), E_BV[b0], 'b.', alpha=.2)
# #plt.plot(up.nominal_values(y[b1]), E_BV[b1], 'r.', alpha=.2)
# #plt.plot(up.nominal_values(y[b2]), E_BV[b2], 'g.', alpha=.2)
# plot_settings(fig, ax, xlab, ylab, 1)
# plt.show()

# IMF for the 3 branches:

# fig, ax = plt.subplots(1,1)
# xlab, ylab = '$M_y$', 'Mass ($M_{\odot}$)'
# ax.errorbar(MY[b0], MM[b0], MM_err[b0], MY_err[b0], 'k.', alpha=.1)
# ax.errorbar(MY[b1], MM[b1], MM_err[b1], MY_err[b1], 'k.', alpha=.1)
# ax.errorbar(MY[b2], MM[b2], MM_err[b2], MY_err[b2], 'k.', alpha=.1)
# plt.plot(MY[b0], MM[b0], 'b.', label='Branch 1', alpha=.3)
# plt.plot(MY[b1], MM[b1], 'r.', label='Branch 2', alpha=.3)
# plt.plot(MY[b2], MM[b2], 'g.', label='Branch 3', alpha=.3)
# #plt.plot(MY[i5], MM[i5], 'm.', label='Branch 3', alpha=1)
# plot_settings(fig, ax, xlab, ylab, 1)
# plt.show()



#--------------
# MASS VS. AGE:
#--------------

# fig, ax = plt.subplots(1,1)
# xlab, ylab =  'Age (Gyr)', 'Mass $(M_{\odot}$)'
# ax.errorbar(val(M), val(age), err(age), err(M), 'r.', ecolor='k', alpha=0.05)
# plot_settings(fig, ax, xlab, ylab, 1)
# plt.show()

#------------------
# AGE VS. DISTANCE:
#------------------

# fig, ax = plt.subplots(1,1)
# xlab, ylab =  'Age (Gyr)', '$d$ (pc)'
# ax.errorbar(val(d), val(age), err(age), err(d), 'b.', ecolor='k', alpha=0.1)
# plt.plot(val(d), val(age), 'b.', alpha=.2)
# plot_settings(fig, ax, xlab, ylab, 1)
# plt.show()

#------------------
# LUMI vs. MASS   :
#------------------




# Transform data into luminosity;
import scaling_relations as scaling
T, _, = scaling.T('torres', 'B-V', by, None, 'dwarf')

BC_V  = scaling.BC_V('torres', np.log10(T))
M_bol = 4.75
L = 10**((My+BC_V-M_bol)/(-2.5))

# Model fit:
m = np.linspace(0.1, 4, 500)
LL = (0.1400*m**2 + 46.81*m**9)/(2.395*10**(-10)*m**(-9) + 1 + 57.99*m**(4.5))

fig, ax = plt.subplots(1,1)
xlab, ylab =  '$\log(M/M_{\odot})$', '$\log(L/L_{\odot})$'
ax.errorbar(val(M[i2]), val(L[i2]), err(L[i2]), err(M[i2]), 'r.', alpha=0.3)
plt.loglog(m, LL, 'k--', alpha=.5)
plot_settings(fig, ax, xlab, ylab, 1)
plt.show()

