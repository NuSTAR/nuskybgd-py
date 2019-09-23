import numpy as np
import astropy.io.fits as pf


def make_aspecthist_img(coords, dt):

    asphistimg = np.zeros((1000, 1000), dtype=np.float64)
    x_min, y_min, x_max, y_max = 999, 999, 0, 0
    for i in range(len(dt)):
        x, y = int(np.floor(coords[i][0])), int(np.floor(coords[i][1]))
        asphistimg[y, x] += dt[i]
        if y < y_min:
            y_min = y
        elif y > y_max:
            y_max = y
        if x < x_min:
            x_min = x
        elif x > x_max:
            x_max = x

    print('Image subspace: ', x_min, x_max, y_min, y_max)

    return asphistimg, x_min, x_max, y_min, y_max


def write_aspecthist_img(asphistimg, outfilename, aimpointext, offsets,
                         overwrite=False):
    import datetime

    x_min, x_max, y_min, y_max = offsets

    # Write some useful info into the header
    aspecthistfh = pf.HDUList(
        pf.PrimaryHDU(asphistimg[y_min:y_max + 1, x_min:x_max + 1])
    )
    aspecthistfh[0].header['EXTNAME'] = 'ASPECT_HISTOGRAM'
    aspecthistfh[0].header['X_OFF'] = (x_min, 'x offset in pixels')
    aspecthistfh[0].header['Y_OFF'] = (y_min, 'y offset in pixels')
    aspecthistfh[0].header['EXPOSURE'] = (
        np.float32(np.sum(asphistimg[y_min:y_max + 1, x_min:x_max + 1])),
        'seconds, total exposure time'
    )
    aspecthistfh[0].header['COMMENT'] = (
        'Add the x/y offset to image coordinates to recover aimpoint '
        'in detector coordinates.'
    )
    aspecthistfh[0].header['DATE'] = (
        datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%S'),
        'File creation date (YYYY-MM-DDThh:mm:ss UTC)'
    )
    for hdrkey in ('TELESCOP', 'INSTRUME', 'OBS_ID', 'OBJECT', 'TARG_ID',
                   'RA_OBJ', 'DEC_OBJ', 'RA_NOM', 'DEC_NOM',
                   'RA_PNT', 'DEC_PNT', 'TIMESYS', 'MJDREFI', 'MJDREFF',
                   'CLOCKAPP', 'TSTART', 'TSTOP', 'DATE-OBS', 'DATE-END'):
        if hdrkey in aimpointext.header:
            aspecthistfh[0].header[hdrkey] = (
                aimpointext.header[hdrkey],
                aimpointext.header.comments[hdrkey]
            )

    aspecthistfh.writeto(outfilename, overwrite=overwrite)
