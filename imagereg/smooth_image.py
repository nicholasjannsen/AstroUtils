

def smooth_image(self, data, smooth):
        """ This function takes all loaded flat-frames, dark-frames, bias-frames and combine them to a one
        master-flat-field image.
        ----------INPUT:
        path           : Directory path to data.
        LF_name        : Name of Light Frames (LF) except end number.
        FF_name        : Name of Flat  Frames (FF) except end number.
        DF_name        : Name of Dark  Frames (DF) except end number.
        BF_name        : Name of Bias  Frames (BF) except end number.
        N              : Number of [LF, FF, DF, BF] frames.
        plot           : If plot=1: plots all relevant frames..
        save           : If save=1: saves all corrected light frames in seperate files.
        ---------OUTPUT:
        CF_i           : Cube of Corrected Frames (CF)."""
        print '--------------------------------------------------------------------- smooth_image'

        # Data definition:
        CF_i = data
        
        # Clean data for hot pixels:
        if smooth=='median':
            HF_i = copy(CF_i)
            HF_i = median_filter(HF_i, size=3)

        # Replace hot and dead pixels only:
        if smooth=='hot':
        HF_i = copy(CF_i)
        n, h, w = shape(CF_i) # Number of frames, Height and Width
        for i in range(n):
            # Find good threshold to locate pixels:
            blur_i = median_filter(CF_i[i], size=2)
            diff_i = CF_i[i] - blur_i
            threshold = 2*std(diff_i)
            # Find the hot pixels (ignoring the edges):
            hot_pixels = nonzero((abs(diff_i[1:-1, 1:-1])>threshold))
            hot_pixels = array(hot_pixels) + 1  # +1 because 1. row and column was ignored
            print 'Frame {}: Hot/Dead pixels: {}'.format(i, len(hot_pixels[1]))
            for y,x in zip(hot_pixels[0], hot_pixels[1]):
                HF_i[i,y,x] = blur_i[y,x]
            # Now get the pixels on the edges (but not the corners):
            # Left and right sides:
            for index in range(1, h-1):
                # Left side:
                med  = median(CF_i[i, index-1:index+2, 0:2])
                diff = abs(CF_i[i, index, 0] - med)
                if diff>threshold: 
                    hot_pixels = hstack((hot_pixels, [[index],[0]]))
                    HF_i[i, index, 0] = med
                # Right side:
                med  = median(CF_i[i, index-1:index+2, -2:])
                diff = abs(CF_i[i, index,-1] - med)
                if diff>threshold: 
                    hot_pixels = hstack((hot_pixels, [[index], [w-1]] ))
                    HF_i[i, index, -1] = med
            # Then the top and bottom:
            for index in range(1, w-1):
                # Bottom:
                med  = median(CF_i[i, 0:2, index-1:index+2])
                diff = abs(CF_i[i, 0, index] - med)
                if diff>threshold: 
                    hot_pixels = hstack((hot_pixels, [[0], [index]] ))
                    HF_i[i, 0, index] = med
                # Top:
                med  = median(CF_i[i, -2:, index-1:index+2])
                diff = abs(CF_i[i, -1, index] - med)
                if diff>threshold: 
                    hot_pixels = hstack((hot_pixels, [[h-1], [index]] ))
                    HF_i[i, -1, index] = med
            # Then the corners:
            # Bottom left:
            med  = median(CF_i[i, 0:2, 0:2])
            diff = abs(CF_i[i, 0, 0] - med)
            if diff>threshold: 
                hot_pixels = hstack((hot_pixels, [[0], [0]] ))
                HF_i[i, 0, 0] = med
            # Bottom right:
            med  = median(CF_i[i, 0:2, -2:])
            diff = abs(CF_i[i, 0, -1] - med)
            if diff>threshold: 
                hot_pixels = hstack((hot_pixels, [[0], [w-1]] ))
                HF_i[i, 0, -1] = med
            # Top left:
            med  = median(CF_i[i, -2:, 0:2])
            diff = abs(CF_i[i, -1, 0] - med)
            if diff>threshold: 
                hot_pixels = hstack((hot_pixels, [[h-1], [0]] ))
                HF_i[i, -1, 0] = med
            # Top right:
            med  = median(CF_i[i, -2:, -2:])
            diff = abs(CF_i[i, -1, -1] - med)
            if diff>threshold: 
                hot_pixels = hstack((hot_pixels, [[h-1], [w-1]] ))
                HF_i[i, -1, -1] = med

    # Save frames in new fits files:
    if save==1:
        if smooth!=None:
            for i in range(1, N[0]+1):
                pyfits.writeto(('{}{}_%03d.fits'.format(path, 'CF') %i), HF_i[i-1], clobber=True)
        else:
            for i in range(1, N[0]+1):
                pyfits.writeto(('{}{}_%03d.fits'.format(path, 'CF') %i), CF_i[i-1], clobber=True)

    print ('Filter done in time: %s s' % (time.time() - start_time))
    return CF_i 
