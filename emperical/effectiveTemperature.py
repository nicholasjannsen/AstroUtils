#!/usr/bin/env python3

"""
Written Marts 2018 - Nicholas Jannsen

This small piece of software can be used to quickly use imperical scaling relations for both observational and fundamental stellar parameters. These relations is directly from literature. An opdate for further use will follow, however, for no the software can be used to find effective temperatures from the observed colors and metallicities, and bolometric corrections from colors only.
"""

def T(method, color='B-V', C=1.0, M=0.0, star='dwarf'):

    if method=='Alonso1999':
        """
        Alonso et al. (1999) - Giant stars:
        These calibrations are only applicable for giant stars with spectral type F0-K5
        """

        if color=='U-V':
            # Eq.1:
            if +0.40<=C<=+1.20 and +2.20>=M>-0.50  or  +0.35<=C<=+1.20 and -0.50>=M>-1.50  or\
               +0.40<=C<=+1.20 and -1.50>=M>-2.50  or  +0.50<=C<=+1.20 and -2.50>=M>-3.00:
                a = [0.6388, 0.4065, -0.1117, 2.308e-3, -7.783e-2, -1.200e-2]
                sigma_Teff = 164
            # Eq. 2:
            if +1.50<=C<=+3.50 and +2.20>=M>-0.50  or  +0.35<=C<=+3.50 and -0.50>=M>-1.50  or\
               +1.50<=C<=+3.25 and -1.50>=M>-2.50:
                a = [0.8323, 9.374e-2, 1.184e-2, -2.351e-2, -0.1392, -1.944e-2]; sigma_Teff = 80

        if color=='B-V':
            # Eq. 3:
            if +0.20<=C<=+0.80 and +0.20>=M>-0.50  or  +0.35<=C<=+0.80 and -0.50>=M>-1.50  or\
               +0.35<=C<=+0.80 and -1.50>=M>-2.50  or  +0.50<=C<=+0.80 and -2.50>=M>-3.00:
                a = [0.5716, 0.5404, -6.126e-2, 4.862e-2, -1.777e-2, -7.969e-3]; sigma_Teff = 167
            # Eq. 4:
            if +0.70<=C<=+1.90 and +0.20>=M>-0.50  or  +0.70<=C<=+1.80 and -0.50>=M>-1.50  or\
               +0.70<=C<=+1.35 and -1.50>=M>-2.50  or  +0.70<=C<=+1.00 and -2.50>=M>-3.00:
                a = [0.6177, 0.4354, -4.025e-3, -5.204e-2, -0.1127, -1.385e-2]; sigma_Teff = 96

        if color=='V-R':
            # Eq. 5:
            if +0.15<=C<=+1.70 and +0.20>=M>-0.50  or  +0.45<=C<=+1.50 and -0.50>=M>-1.50  or\
               +0.50<=C<=+1.00 and -1.50>=M>-2.50  or  +0.55<=C<=+0.85 and -2.50>=M>-3.00:
                a = [0.4972, 0.8841, -0.1904, 1.197e-2, -1.025e-2, -5.500e-3]; sigma_Teff = 150

        if color=='V-I':
            # Eq. 6:
            if +0.20<=C<=+2.90 and +0.20>=M>-0.50  or  +0.80<=C<=+2.00 and -0.50>=M>-1.50  or\
               +0.85<=C<=+2.20 and -1.50>=M>-2.50  or  +1.00<=C<=+1.70 and -2.50>=M>-3.00:
                a = [0.5379, 0.3981, 4.432e-2, 0.0, 0.0, 0.0, -2.693e-2]; sigma_Teff = 125

        if color=='R-I':
            # Eq. 7:
            if +0.15<=C<=+1.40 and +0.20>=M>-0.50  or  +0.25<=C<=+0.80 and -0.50>=M>-1.50  or\
               +0.35<=C<=+0.70 and -1.50>=M>-2.50  or  +0.40<=C<=+0.65 and -2.50>=M>-3.00:
                a = [0.4974, 1.345, -0.5008, 8.134e-2, 3.705e-2, -6.184e-3]; sigma_Teff = 150

        if color=='V-K':
            # Eq. 8:
            if +0.20<=C<=+2.50 and +0.20>=M>-0.50  or  +1.00<=C<=+2.50 and -0.50>=M>-1.50  or\
               +1.20<=C<=+2.50 and -1.50>=M>-2.50  or  +1.70<=C<=+2.50 and -2.50>=M>-3.00:
                a = [0.5558, 0.2105, 1.981e-3, 9.965e-3, 1.325e-2, -2.726e-3]; sigma_Teff = 40
            # Eq. 9:
            if +2.00<=C<=+4.90 and +0.20>=M>-0.50  or  +2.00<=C<=+4.60 and -0.50>=M>-1.50  or\
               +2.00<=C<=+3.40 and -1.50>=M>-2.50  or  +2.00<=C<=+2.80 and -2.50>=M>-3.00:
                a = [0.3770, 0.3660, -3.170e-2, 3.074e-3, -2.765e-3, -2.973e-3]; sigma_Teff = 25

        if color=='J-H':
            # Eq. 10:
            if +0.00<=C<=+0.90 and +0.20>=M>-0.50  or  +0.20<=C<=+0.80 and -0.50>=M>-1.50  or\
               +0.30<=C<=+0.70 and -1.50>=M>-2.50  or  +0.35<=C<=+0.65 and -2.50>=M>-3.00:
                a = [0.5977, 1.015, 1.020e-1, 1.029e-2, 3.006e-2, 1.013e-2]; sigma_Teff = 170

        if color=='J-K':
            # Eq. 11:
            if +0.00<=C<=+1.10 and +0.20>=M>-0.50  or  +0.20<=C<=+1.00 and -0.50>=M>-1.50  or\
               +0.30<=C<=+0.90 and -1.50>=M>-2.50  or  +0.40<=C<=+0.80 and -2.50>=M>-3.00:
                a = [0.5816, 0.9134, -0.1443]; sigma_Teff = 125

        if color=='V-L':
            # Eq. 12:
            if +0.40<=C<=+5.00 and +0.20>=M>-0.50:
                a = [0.5641, 1882, 1.890e-2, 0.0, 0.0, 0.0, 4.651e-3]; sigma_Teff = 65

        if color=='I-K':
            # Eq. 13:
            if +0.00<=C<=+1.90 and +0.20>=M>-0.50  or  +0.50<=C<=+1.60 and -0.50>=M>-1.50  or\
               +0.50<=C<=+1.50 and -1.50>=M>-2.50  or  +0.80<=C<=+1.20 and -2.50>=M>-3.00:
                a = [0.5859, 0.4846, -2.457e-2]; sigma_Teff = 130

        if color=='b-y':
            # Eq. 14:
            if +0.00<=C<=+0.55 and +0.20>=M>-0.50  or  +0.30<=C<=+0.55 and -0.50>=M>-1.50  or\
               +0.35<=C<=+0.55 and -1.50>=M>-2.50  or  +0.40<=C<=+0.55 and -2.50>=M>-3.00:
                a = [0.5815, 0.7263, 6.856e-2, 6.832e-2, -1.062e-2, -1.079e-2]; sigma_Teff = 110
            # Eq. 15:
            if +0.50<=C<=+1.00 and +0.20>=M>-0.50  or  +0.50<=C<=+0.90 and -0.50>=M>-1.50  or\
               +0.50<=C<=+0.80 and -1.50>=M>-2.50  or  +0.50<=C<=+0.70 and -2.50>=M>-3.00:
                a = [0.4399, 1.209, -0.3541, -8.443e-2, -0.1063, -1.686e-2]; sigma_Teff = 70

        if color=='u-b':
            # Eq. 16:
            if +1.60<=C<=+4.00 and +0.20>=M>-0.50  or  +1.60<=C<=+3.70 and -0.50>=M>-1.50  or\
               +1.60<=C<=+3.40 and -1.50>=M>-2.50  or  +1.60<=C<=+2.60 and -2.50>=M>-3.00:
                a = [0.5883, 0.2008, -5.931e-3, -5.319e-3, -1.000e-1, -1.542e-2]; sigma_Teff = 110


    if method=='Remirez2005':
        """
        Remirez and Melendez (2005) - Giant stars:
        These calibrations are only applicable for giant stars with spectral type F0-K5
        """

        if color=='B-V':
            if star=='dwarf':
                a = [0.5002, 0.6440, -0.0690, -0.0230, -0.0566, -0.0170]; sigma_Teff = 88
                if 0.310<=C<=1.507 and -0.5<M<+0.5:  P = [-261.548, 684.977, -470.049, 79.8977]
                if 0.307<=C<=1.202 and -1.5<M<=-0.5: P = [-324.033, 1516.44, -2107.37, 852.150]
                if 0.335<=C<=1.030 and -2.5<M<=-1.5: P = [30.5985, -46.7882]
                if 0.343<=C<=0.976 and -4.0<M<=-2.5: P = [139.965, -292.329]
            if star=='giant':
                a = [0.5737, 0.4882, -0.0149, 0.0563, -0.1160, -0.0114]; sigma_Teff = 51
                if 0.144<=C<=1.668 and -0.5<M<+0.5:
                    P = [112.116, -372.622, 67.1254, 395.333, -203.471]
                if 0.664<=C<=1.558 and -1.5<M<=-0.5: P = [-12.9762]
                if 0.605<=C<=1.352 and -2.5<M<=-1.5: P = [606.032, -1248.79, 627.453]
                if 0.680<=C<=1.110 and -4.0<M<=-2.5: P = [-9.26209]

        if color=='b-y':
            if star=='dwarf':
                a = [0.4129, 1.2570, -0.2268, -0.0242, -0.0464, -0.0200]; sigma_Teff = 87
                if 0.248<=C<=0.824 and -0.5<M<+0.5:  P = [-1237.11, 6591.29, -11061.3, 5852.18]
                if 0.234<=C<=0.692 and -1.5<M<=-0.5:
                    P = [-2617.66, 22607.4, -68325.4, 86072.5, -38602.2]
                if 0.290<=C<=0.672 and -2.5<M<=-1.5: P = [103.927, -312.419, 225.430]
                if 0.270<=C<=0.479 and -4.0<M<=-2.5: P = [-294.106, 648.320]
            if star=='giant':
                a = [0.5515, 0.9085, -0.1494, 0.0616, -0.0668, -0.0083]; sigma_Teff = 68
                if 0.053<=C<=1.077 and -0.5<M<+0.5:  P = [-124.159, 553.827, -490.703]
                if 0.309<=C<=0.893 and -1.5<M<=-0.5: P = [888.088, -2879.23, 2097.89]
                if 0.388<=C<=0.702 and -2.5<M<=-1.5: P = [1867.63, -6657.49, 5784.81]
                if 0.404<=C<=0.683 and -4.0<M<=-2.5: P = [348.237, -659.093]

        if color=='V-R':
            if star=='dwarf':
                a = [0.433, 1.4399, -0.5419, -0.0481, -0.0239, -0.0125]; sigma_Teff = 84
                if 0.204<=C<=0.880 and -0.5<M<+0.5:
                    P = [-2666.55, 27264.5, -103923.0, 174663.0, -104940.0, -23249.4, 32644.9]
                if 0.284<=C<=0.546 and -1.5<M<=-0.5: P = [4.20153]
                if 0.264<=C<=0.532 and -2.5<M<=-1.5: P = [123.940, -342.217]
                if 0.240<=C<=0.336 and -4.0<M<=-2.5: P = [8.55498]
            if star=='giant':
                a = [0.3849, 1.6205, -0.6395, 0.1060, -0.0875, -0.0089]; sigma_Teff = 41
                if 0.299<=C<=1.106 and -0.5<M<+0.5:  P = [-8.51797, 15.6675]
                if 0.387<=C<=0.752 and -1.5<M<=-0.5: P = [-10.7764]
                if 0.429<=C<=0.598 and -2.5<M<=-1.5: P = [61.9821, -78.7382]
                if 0.394<=C<=0.550 and -4.0<M<=-2.5: P = [27.9886, -100.149]

        if color=='V-I':
            if star=='dwarf':
                a = [0.3295, 0.9516, -0.2290, -0.0316, 0.0003, -0.0081]; sigma_Teff = 68
                if 0.491<=C<=1.721 and -0.5<M<+0.5:
                    P = [-2757.79, 9961.33, -10546.6, -1746.05, 10512.3, -6653.57, 1301.21]
                if 0.597<=C<=1.052 and -1.5<M<=-0.5: P = [-22.9008, 40.2078]
                if 0.547<=C<=1.026 and -2.5<M<=-1.5: P = [-667.732, 1709.88, -1069.62]
            if star=='giant':
                a = [0.3575, 0.9069, -0.2025, 0.0395, -0.0551, -0.0061]; sigma_Teff = 40
                if 0.573<=C<=2.000 and -0.5<M<+0.5:  P = [0.42933]
                if 0.795<=C<=1.524 and -1.5<M<=-0.5: P = [-0.14180]
                if 0.870<=C<=1.303 and -2.5<M<=-1.5: P = [9.31011]
                if 0.812<=C<=1.095 and -4.0<M<=-2.5: P = [-23.0514]

        if color=='R-I':
            if star=='dwarf':
                a = [0.2919, 2.1141, -1.0723, -0.0756, 0.0267, -0.0041]; sigma_Teff = 76
                if 0.242<=C<=0.838 and -0.5<M<+0.5:
                    P = [-3326.97, 26263.8, -75355.8, 94246.5, -43334.8]
                if 0.300<=C<=0.718 and -1.5<M<=-0.5: P = [12.4740]
                if 0.283<=C<=0.551 and -2.5<M<=-1.5: P = [-5837.31, 41439.2, -94729.8, 69584.8]
                if 0.290<=C<=0.364 and -4.0<M<=-2.5: P = [32.1826]
            if star=='giant':
                a = [0.4351, 1.6549, -0.7215, -0.0610, 0.0332, -0.0023]; sigma_Teff = 62
                if 0.413<=C<=0.793 and -0.5<M<+0.5:  P = [61.3557, -116.711]
                if 0.383<=C<=0.771 and -1.5<M<=-0.5: P = [-16.8645]
                if 0.434<=C<=0.725 and -2.5<M<=-1.5: P = [32.0870]
                if 0.364<=C<=0.545 and -4.0<M<=-2.5: P = [-15.6318]

        if color=='V-J':
            if star=='dwarf':
                a = [0.4050, 0.4792, -0.0617, -0.0392, 0.0401, -0.0023]; sigma_Teff = 62
                if 0.815<=C<=2.608 and -0.5<M<+0.5:  P = [422.406, -910.603, 621.335, -132.566]
                if 0.860<=C<=2.087 and -1.5<M<=-0.5: P = [-466.616, 658.349, -220.454]
                if 0.927<=C<=1.983 and -2.5<M<=-1.5: P = [-862.072, 1236.84, -423.729]
                if 0.891<=C<=1.932 and -4.0<M<=-2.5: P = [-1046.10, 1652.06, -597.340]
            if star=='giant':
                a = [0.2943, 0.5604, -0.0677, 0.0179, -0.0532, -0.0088]; sigma_Teff = 38
                if 1.259<=C<=2.400 and -0.5<M<+0.5:  P = [-122.595, 76.4847]
                if 1.030<=C<=3.418 and -1.5<M<=-0.5: P = [-10.3848]
                if 1.033<=C<=2.679 and -2.5<M<=-1.5: P = [4.18695]
                if 0.977<=C<=2.048 and -4.0<M<=-2.5: P = [-67.7716, 28.9202]

        if color=='V-H':
            if star=='dwarf':
                a = [0.4931, 0.3056, -0.0241, -0.0396, 0.0678, 0.0020]; sigma_Teff = 57
                if 0.839<=C<=3.215 and -0.5<M<+0.5:  P = [-53.5574, 36.0990, 15.6878, -8.84468]
                if 1.032<=C<=2.532 and -1.5<M<=-0.5: P = [1.60629]
                if 1.070<=C<=2.535 and -2.5<M<=-1.5: P = [506.559, -1277.52, 939.519, -208.621]
                if 1.093<=C<=2.388 and -4.0<M<=-2.5: P = [-471.588, 643.972, -199.639]
            if star=='giant':
                a = [0.4354, 0.3405, -0.0263, -0.0012, -0.0049, -0.0027]; sigma_Teff = 32
                if 1.194<=C<=3.059 and -0.5<M<+0.5:  P = [-377.022, 334.733, -69.8093]
                if 1.293<=C<=4.263 and -1.5<M<=-0.5: P = [71.7949, -55.5383, 9.61821]
                if 1.273<=C<=3.416 and -2.5<M<=-1.5: P = [-27.4190, 20.7082]
                if 1.232<=C<=2.625 and -4.0<M<=-2.5: P = [-46.2946, 20.1061]

        if color=='V-K':
            if star=='dwarf':
                a = [0.4942, 0.2809, -0.0180, -0.0294, 0.0444, -0.0008]; sigma_Teff = 50
                if 0.895<=C<=3.360 and -0.5<M<+0.5:
                    P = [-1425.36, 3218.36, -2566.54, 859.644, -102.554]
                if 1.060<=C<=2.665 and -1.5<M<=-0.5: P = [2.35133]
                if 1.101<=C<=2.670 and -2.5<M<=-1.5:
                    P = [-1849.46, 4577.00, -4284.02, 1770.38, -268.589]
                if 1.126<=C<=2.596 and -4.0<M<=-2.5: P = [215.721, -796.519, 714.423, -175.678]
            if star=='giant':
                a = [0.4405, 0.3272, -0.0252, -0.0016, -0.0053, -0.0040]; sigma_Teff = 28
                if 1.244<=C<=3.286 and -0.5<M<+0.5:  P = [-72.6664, 36.5361]
                if 1.366<=C<=4.474 and -1.5<M<=-0.5: P = [86.0358, -65.4928, 10.8901]
                if 1.334<=C<=3.549 and -2.5<M<=-1.5: P = [-6.96153, 14.3298]
                if 1.258<=C<=2.768 and -4.0<M<=-2.5: P = [-943.925, 1497.64, -795.867, 138.965]


    if method=='Torres2010':
        """
        Guillermo Torees (2010) - Dwarf stars:
        This relations is independent of metallicity. Notice that 'dwarf' is
        main-sequence stars, subgiants, and giants, where 'giant' is refered
        to supergiants.
        """

        if color=='B-V':
            if star=='dwarf':
                sigma_Teff = 0
                a = [3.979145106714099, -0.654992268598245, 1.740690042385095, -4.608815154057166, \
                     6.792599779944473, -5.396909891322525, 2.192970376522490, -0.359495739295671]
            if star=='giant':
                sigma_Teff = 0
                a = [4.012559732366214, -1.055043117465989, 2.133394538571825, -2.459769794654992, \
                     1.349423943497744, -0.283942579112032]


    if method=='Sousa2008':
        """
        Sousa et al. (2008) - :
        Relation for dwarf stars observed with HARP
        """
        if color=='B-V':
            if 0.510<=C<=1.200 and -0.85<M<+0.40: a = [9114, -6827, 2638, 368]; sigma_Teff = 47



    #-------------------------------------#
    #          MODEL FIT TO DATA          #
    #-------------------------------------#

    # If outside applicable ranges:
    if method=='alonso' or method=='torres' or method=='sausa':
        try: a
        except NameError:
            print('ERROR: COLOR OR METALLICITY IS NOT INSIDE APPLICABLE RANGE!'); exit(1)
    if method=='remirez':
        try: P
        except NameError:
            print('ERROR: COLOR OR METALLICITY IS NOT INSIDE APPLICABLE RANGE!'); exit(1)

    # Calculate Teff if inside applicable range:
    # FIXME: rewrite theta_eff such this is not needed!
    if len(a)==3: theta_eff = a[0] + a[1]*C + a[2]*C**2
    if len(a)==4: theta_eff = a[0] + a[1]*C + a[2]*C**2 + a[3]*M
    if len(a)==6: theta_eff = a[0] + a[1]*C + a[2]*C**2 + a[3]*C*M + a[4]*M + a[5]*M**2
    if len(a)==7: theta_eff = a[0] + a[1]*C + a[2]*C**2 + a[3]*C*M + a[4]*M + a[5]*M**2+a[6]*C**3

    # Using Alonso:
    if method=='alonso':
        Teff = 5040/theta_eff
    # Using Remirez:
    if method=='remirez':
        PP = np.zeros(len(P))
        for i in range(len(P)): PP[i] = P[i]*C**i
        Teff = 5040/theta_eff + sum(PP)
    # Using Torres:
    if method=='torres':
        logT = np.zeros((len(a), len(C)))
        for i in range(len(a)): logT[i] = a[i]*C**i
        Teff = 10**sum(logT)
    # Using Sausa:
    if method=='sausa':
        for i in range(len(P)): PP[i] = P[i]*C**i
        Teff = theta_eff

    # Print results:
    #print('Teff ({}) = {} K +/- {}'.format(color, Teff, sigma_Teff))
    return Teff, sigma_Teff
