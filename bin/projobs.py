#!/usr/bin/env python3
"""
projobs.py

Create an aspect histogram image from pointing info after filtering by GTI.

projobs.py aimpoints.fits gtifile=gti.fits out=aspecthist.fits

Example:

projobs.py nu90201039002A_det1.fits gtifile=nu90201039002A01_gti.fits \\
    out=aspecthistA.fits

projobs.py nu90201039002B_det1.fits gtifile=nu90201039002B01_gti.fits \\
    out=aspecthistB.fits

The output file has an image representing the 2D histogram of pointing
position with time. Gets pointing position from nu%obsid%?_det?.fits and good
time intervals from nu%obsid%?0?_gti.fits. nu%obsid%?_det?.fits (e.g.
nu90201039002A_det1.fits) has a table block named 'DET?_REFPOINT', with 3
columns (TIME, X_DET?, Y_DET?) where ? is the detector number (1-4).
nu%obsid%?0?_gti.fits (e.g. nu90201039002A01_gti.fits) has a table block named
'STDGTI' with two columns (START, STOP) listing intervals of good times.

The image is in an extension with the name ASPECT_HISTOGRAM. Any zero padding
has been cropped, and the x and y offsets are recorded in the header keywords
X_OFF and Y_OFF.
"""

import astropy.io.fits as pf
import os
import sys
import numpy as np


def check_header_aimpoint(hdr):
    """
    Check a FITS header to see if it does not have the expected data table
    columns for telescope aspect.

    Require the following keywords:

    TELESCOP = 'NuSTAR'
    INSTRUME = 'FPM?'           # ? is 'A' or 'B'
    EXTNAME  = 'DET?_REFPOINT'  # ? is a digit 1 to 4
    HDUCLAS1 = 'TEMPORALDATA'
    NAXIS    = 2

    Require the following fields:

    TTYPE1   = 'TIME'
    TUNIT1   = 's'
    TTYPE2   = 'X_DET1'
    TUNIT2   = 'pixel'
    TTYPE3   = 'Y_DET1'
    TUNIT3   = 'pixel'
    """
    try:
        if not (
            hdr['TELESCOP'] == 'NuSTAR' and
            (hdr['INSTRUME'] == 'FPMA' or hdr['INSTRUME'] == 'FPMB') and
            (hdr['EXTNAME'] in
                ['DET%d_REFPOINT' % i for i in (1, 2, 3, 4)]) and
            hdr['HDUCLAS1'] == 'TEMPORALDATA' and
            hdr['NAXIS'] == 2 and
            hdr['TTYPE1'] == 'TIME' and hdr['TUNIT1'] == 's' and
            (hdr['TTYPE2'] in
                ['X_DET%d' % i for i in (1, 2, 3, 4)]) and
            hdr['TUNIT2'] == 'pixel' and
            (hdr['TTYPE3'] in
                ['Y_DET%d' % i for i in (1, 2, 3, 4)]) and
            hdr['TUNIT3'] == 'pixel'
        ):
            return False
    except KeyError:
        return False

    return True


def check_header_gti(hdr):
    """
    Check FITS header to see if it has the expected data table columns for GTI
    extension.

    Require the following keywords:

    TELESCOP = 'NuSTAR'
    EXTNAME  = 'STDGTI'  # ? is a digit 1 to 4
    NAXIS    = 2

    Require the following fields:

    TTYPE1   = 'START'
    TUNIT1   = 's' or 'sec'
    TTYPE2   = 'STOP'
    TUNIT2   = 's' or 'sec'
    """
    try:
        if not (
            hdr['TELESCOP'] == 'NuSTAR' and
            hdr['EXTNAME'] == 'STDGTI' and
            hdr['NAXIS'] == 2 and
            hdr['TTYPE1'] == 'START' and
            hdr['TUNIT1'] in ('s', 'sec') and
            hdr['TTYPE2'] == 'STOP' and
            hdr['TUNIT2'] in ('s', 'sec')
        ):
            return False
    except KeyError:
        return False

    return True


def filter_gti(aimpointext, gtiext):

    gtiarr = np.sort(gtiext.data, order='START', kind='mergesort')
    aimpointarr = np.sort(aimpointext.data, order='TIME', kind='mergesort')

    n_aimpoint = len(aimpointarr)
    i_gti = 0
    flagarr = np.zeros(n_aimpoint, dtype=np.int8)

    print('Processing %d aimpoints...' % n_aimpoint)

    for i_aimpoint in range(n_aimpoint - 1):
        # Ref pointings have a single time while GTI intervals have two.
        # Original logic in projobs.pro considers pointing time against
        # [start, stop), i.e. >= for start time, < for stop.

        aim_time = aimpointarr[i_aimpoint][0]

        # Ensure current GTI interval ends after current pointing time
        try:
            while aim_time >= gtiarr[i_gti][1]:
                i_gti += 1
        except IndexError:
            # Reached the end of the GTI list
            break

        # Flag this pointing if it is inside current GTI interval
        flagarr[i_aimpoint] = int(aim_time >= gtiarr[i_gti][0])

    flagged = np.where(flagarr == 1)[0]
    view1 = aimpointarr[flagged]
    view2 = aimpointarr[flagged + 1]
    coords = np.array([view1['X_DET1'], view1['Y_DET1']]).T
    dt = view2['TIME'] - view1['TIME']

    print('%d aimpoints (%f / %f s)' % (
        len(coords), np.sum(dt),
        np.sum(gtiarr.field('STOP') - gtiarr.field('START'))))

    return coords, dt


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
    aspecthistfh[0].header['HISTORY'] = (
        'Aspect histogram image created by filtering %s using %s.' % (
            aimpointfile, gtifile
        )
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


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(0)

    aimpointfile = sys.argv[1]

    keywords = {
        'gtifile': None,
        'out': None
    }

    for _ in sys.argv[1:]:
        arg = _.split('=')
        if arg[0] in keywords:
            keywords[arg[0]] = arg[1]

    gtifile = keywords['gtifile']
    outfilename = keywords['out']

    # Check the arguments
    scriptname = os.path.basename(__file__)
    print('Running %s. Performing checks...' % scriptname)

    halt = False

    if not os.path.exists(aimpointfile):
        print('aimpointfile not found: %s' % aimpointfile)
        halt = True

    if gtifile is None:
        print('gtifile= missing')
        halt = True
    elif not os.path.exists(gtifile):
        print('GTI file not found: %s' % gtifile)
        halt = True

    if outfilename is None:
        print('out= missing')
        halt = True
    elif os.path.exists(outfilename) and keywords['out'][0] != '!':
        print('Output file exists: %s' % outfilename)
        halt = True

    if halt:
        print('%s did not complete. See error messages.' % scriptname)
        sys.exit(1)

    # Check contents of aimpoint and gti files
    aimpointext = None
    detfh = pf.open(aimpointfile)
    for ext in detfh:
        if check_header_aimpoint(ext.header):
            aimpointext = ext
            break
    if aimpointext is None:
        print('No aspect info in the specified file %s.' % aimpointfile)
    else:
        print('Found extension %s.' % aimpointext.header['EXTNAME'])

    gtiext = None
    gtifh = pf.open(gtifile)
    for ext in gtifh:
        if check_header_gti(ext.header):
            gtiext = ext
            break
    if gtiext is None:
        print('No GTI info in the specified file %s.' % gtifile)
    else:
        print('Found extension %s.' % gtiext.header['EXTNAME'])

    # Check output file
    if os.path.exists(outfilename):
        if outfilename[0] == '!':
            outfilename = outfilename[1:]
            print('Overwriting file %s...')
        else:
            print('Output file %s exists (prefix with ! to overwrite)')
            sys.exit(1)

    coords, dt = filter_gti(aimpointext, gtiext)

    asphistimg, x_min, x_max, y_min, y_max = make_aspecthist_img(coords, dt)

    write_aspecthist_img(
        asphistimg, outfilename, aimpointext,
        (x_min, x_max, y_min, y_max), overwrite=True)

    sys.exit(0)
