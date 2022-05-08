#!/usr/bin/env python3

"""
Written Marts 2018 - Nicholas Jannsen

This small piece of software can be used to quickly use imperical scaling relations for both observational and fundamental stellar parameters. These relations is directly from literature. An opdate for further use will follow, however, for no the software can be used to find bolometric corrections from colors only.
"""

import numpy as np

def BC_V(method, logT):
    """
    The following is a empirical calibration of the bolometric correction in the V-band (BC_V). 
    @param method         : At the moment only Torres results can be used. Typing 'info' print function info.   
    @param logT
    BC_V           : Bolometric corrections
    """

    if method=='torres':
        """
        Guillermo Torres (2010)
        """

        # Sort stars after effective temperature:
        i0 = (logT < 3.70)
        i1 = (logT > 3.70) * (logT <= 3.90)
        i2 = (logT > 3.90)
        # Make list of stars:
        logT0 = logT[i0]
        logT1 = logT[i1]
        logT2 = logT[i2]

        # Coefficients:
        a0 = [-19053.72914964556, 15514.4866764412, -4212.78819301717, 381.476328422343]
        a1 = [-37051.0203809015,  38567.2629965804, -15065.1486316025, 2617.24637119416, -170.623810323864]
        a2 = [-118115.450538963,  137145.973583929, -63623.3812100225, 14741.2923562646, -1705.87278406872, \
              78.8731721804990]
        # If it is calculate:
        BC_V0 = np.zeros((len(a0), len(logT0)))
        BC_V1 = np.zeros((len(a1), len(logT1)))
        BC_V2 = np.zeros((len(a2), len(logT2)))
        for i in range(len(a0)): BC_V0[i] = a0[i]*logT0**i
        for i in range(len(a1)): BC_V1[i] = a1[i]*logT1**i
        for i in range(len(a2)): BC_V2[i] = a2[i]*logT2**i
        # Sum them togther:
        BC_V0 = np.sum(BC_V0, axis=0)
        BC_V1 = np.sum(BC_V1, axis=0)
        BC_V2 = np.sum(BC_V2, axis=0)
        # Insert in the right order:
        BC_V = np.zeros(len(logT))
        BC_V.flat[np.where(i0)[0]] = BC_V0
        BC_V.flat[np.where(i1)[0]] = BC_V1
        BC_V.flat[np.where(i2)[0]] = BC_V2

    # PRINT FINAL RESULTS:
    #print 'BC_V = {}'.format(BC_V)
    return BC_V
