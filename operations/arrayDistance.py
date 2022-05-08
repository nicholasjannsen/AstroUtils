#!/usr/bin/env python3

import numpy as np
from scipy.spatial.distance import cdist


def arrayDistance(array0, array1, threshold, mask=False):
    """
    Function to check the distance between all the points in two arrays of coordinates
    """
    # cdist returns 1 if match or 0 if new star:
    dist = cdist(array0, array1)

    # Perform distance measure using a mask
    if mask is True:
        mdist = np.ma.masked_where(np.tril(dist)==0, dist)
        indices = np.nonzero(mdist < threshold)

    # Perform distance measure without a mask
    elif mask is False: indices = np.nonzero(dist < threshold)

    # Output
    return np.transpose(indices)
