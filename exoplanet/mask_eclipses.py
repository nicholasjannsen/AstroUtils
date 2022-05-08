def mask_eclipses(self, time, flux):
    """
    Function to mask the primary and secondary eclipse for the avoidance of
    excluding data points at eclipse.
    """
    # In case the first transit is removed:
    self.t0 = self.t0_eph + int(val( (time[0] - self.t0_eph) / self.P )) * self.P
    self.orbit_no = int(val( (time[-1] - time[0]) /self.P )) + 1

    # Simplify expressions:
    P = val(self.P); t0 = val(self.t0); T2 = val(self.T2)

    # Start masking transits and out-of-transit:
    time_mask_tra = []; time_mask_occ = []; time_mask_out = []
    flux_mask_tra = []; flux_mask_occ = []; flux_mask_out = []

    for j in range(self.orbit_no+1):

        # Transits:
        index_tra = np.where(np.logical_and( (time > t0 + j*P - T2),\
                                             (time < t0 + j*P + T2) ))
        if (len(time[index_tra]) > 0):
            time_mask_tra = np.concatenate((time_mask_tra, time[index_tra]))
            flux_mask_tra = np.concatenate((flux_mask_tra, flux[index_tra]))

        # Occultations:
        index_occ = np.where(np.logical_and((time > t0 + j*P + P/2 - T2),\
                                            (time < t0 + j*P + P/2 + T2) ))
        if (len(time[index_occ]) > 0):
            time_mask_occ = np.concatenate((time_mask_occ, time[index_occ]))
            flux_mask_occ = np.concatenate((flux_mask_occ, flux[index_occ]))

        # Out-of-eclipses:
        index_tra_out = np.where(np.logical_and( (time > t0 + j*P       + T2 ),\
                                                 (time < t0 + j*P + P/2 - T2 )))
        index_occ_out = np.where(np.logical_and( (time > t0 + j*P + P/2 + T2 ),\
                                                 (time < t0 + (j+1)*P   - T2 )))
        index_out = np.concatenate((index_tra_out + index_occ_out))
        if (len(time[index_out]) > 0):
            time_mask_out = np.concatenate((time_mask_out, time[index_out]))
            flux_mask_out = np.concatenate((flux_mask_out, flux[index_out]))

    # Save out-of-eclipses:
    self.time_mask_out = time_mask_out
    self.flux_mask_out = flux_mask_out
    # # Save transits:
    self.time_mask_tra = time_mask_tra
    self.flux_mask_tra = flux_mask_tra
    # Save occultations:
    self.time_mask_occ = time_mask_occ
    self.flux_mask_occ = flux_mask_occ
    # Save eclipses (transits + occultations):
    self.time_mask_in = np.concatenate((self.time_mask_tra, self.time_mask_occ))
    self.flux_mask_in = np.concatenate((self.flux_mask_tra, self.flux_mask_occ))
