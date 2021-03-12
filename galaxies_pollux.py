""" Pollux """

# Package requirements:
import numpy as np
import uncertainties.unumpy as up
from uncertainties import ufloat    # (nominal_value, std_dev)
import sys

# Units in cgs:
c     = 2.99792458*1e10 # [cm/s]
sigma = 5.67051*1e-5    # []
G     = 6.67259*1e-8    # [cm^3/(g*s^2)]
# Solar units in cgs:
Msun  = 1.989*1e33      # [g]
Rsun  = 6.96392*1e10    # [cm]
Lsun  = 3.846*1e33      # [erg/s]
Tsun  = 5776.0          # [K]
gsun  = 27542.29        # [] log(g)=4.44
Msun_bol = 4.75
# Solar system units in cgs:
Mearth = 5.976*1e27     # [g]
Rearth = 6.397*1e8      # [cm]
Mjup = 1.899*1e30       # [g]
Rjup = 7.140*1e9        # [cm]
# Astronomical units in cgs:
AU_cm = 1.496*1e13         # [cm]
pc_cm = 3.086*1e18         # [cm]
ly_cm = 9.463*1e17         # [cm]
rad_as = 206265            # [as]

# Reddening corrected observables:
U = 3.0
B = 2.14
V = 1.14
R = 0.30
I = -0.11
J = -0.52
H = -1.00
K = -1.11

# Observables for pollux:
pi  = ufloat(96.54, 0.27)*1e-3           # [as] parallax
theta = ufloat(8.8, 0.1)*1e-3          # [as] angular diameter  
m   = np.mean((0.08, 0.19,-0.07))      # Metallicity mean value
g   = 10**ufloat(2.685, 0.09)          # Surface gravity

##############################################################################

#----- CHECK PARALLXE FROM: Nordgren et al. (2002)
R_n = ufloat(8.8, 0.1)*Rsun  # [cm]
d_n = 2*R_n/(theta/rad_as)     # [cm]
p_n = 1/(d_n/pc_cm)*1e3        # [mas]

#----- TEMPERATURE FROM SCALING RELATIONS:
from T_Scaling_Relation import T_scaling_relation
# Usinf Alonso:
# Ta1 = T_scaling_relation('alonso', 'U-V', U-V, m)
# Ta2 = T_scaling_relation('alonso', 'B-V', B-V, m)
# Ta3 = T_scaling_relation('alonso', 'V-R', V-R, m)
# Ta4 = T_scaling_relation('alonso', 'V-I', V-I, m)
# Ta5 = T_scaling_relation('alonso', 'R-I', R-I, m)
# Ta6 = T_scaling_relation('alonso', 'V-K', V-K, m)
# Ta7 = T_scaling_relation('alonso', 'J-H', J-H, m)
# Ta8 = T_scaling_relation('alonso', 'J-K', J-K, m)
# Ta9 = T_scaling_relation('alonso', 'I-K', I-K, m)
# # Using Remirez:
# Tm1 = T_scaling_relation('remirez', 'B-V', B-V, m, 'giant')
# Tm2 = T_scaling_relation('remirez', 'V-R', V-R, m, 'giant')
# Tm3 = T_scaling_relation('remirez', 'V-I', V-I, m, 'giant')
# Tm4 = T_scaling_relation('remirez', 'R-I', R-I, m, 'giant')
# Tm5 = T_scaling_relation('remirez', 'V-J', V-J, m, 'giant')
# Tm6 = T_scaling_relation('remirez', 'V-H', V-H, m, 'giant')
# Tm7 = T_scaling_relation('remirez', 'V-K', V-K, m, 'giant')
# Value we adopt:

T = (ufloat(4867, 80)+\
      ufloat(4855, 96)+\
      ufloat(4830, 55)+\
      ufloat(4559, 150)+\
      ufloat(3896, 41)+\
      ufloat(4790, 125)+\
      ufloat(4292, 40)+\
      ufloat(5199, 150)+\
      ufloat(4867, 38)+\
      ufloat(4851, 32)+\
      ufloat(4844, 25)+\
      ufloat(4816, 28)+\
      ufloat(4537, 170)+\
      ufloat(4709, 125)+\
      ufloat(4818, 130))/15


# Mass using R and g:
L = ufloat(43, 0.1)*Lsun
R = up.sqrt((L/Lsun) * (T/Tsun)**-4)*Rsun
M1 = g*R**2/G

# Bolmetric correction:
a = -0.190537291496456*1e5
b = 0.155144866764412*1e5
c = -0.42127881930171*1e4
d = 0.38147632842234*1e3
logT = up.log10(T)
BCv = a + b*logT + c*logT**2 + d*logT**3 

# mass:
M2 = 10**(up.log10(g/gsun) - 4*up.log10(T/Tsun) - 0.4*V - 0.4*BCv - 2*up.log10(pi) - 0.12)*Msun

# luminosity:
d     = 1/pi
M_V   = V - 5*up.log10(1/pi) + 5
M_bol = BCv + M_V
L = Lsun*10**((M_bol-Msun_bol)/(-2.5))

print T
print R/Rsun
print M1/Msun
print M2/Msun
print d
print M_V
print M_bol
print L/Lsun
