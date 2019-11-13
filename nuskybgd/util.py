"""
Utility functions for sharing.
"""
import datetime
import pyregion
import numpy as np


def format_timestamp():
    return datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%S')


def hdu_timestamp_write(hdu):
    """
    Update the DATE keyword of the input HDU header to current time. This
    modifies the input object 'hdu'.
    """
    timestamp = format_timestamp()
    hdu.header['DATE'] = (
        timestamp,
        'File creation date (YYYY-MM-DDThh:mm:ss UTC)'
    )


def hdu_history_write(hdu, message):
    """
    Append HISTORY keyword to the input HDU header. This modifies the input
    object 'hdu'.
    """
    hdu.header['HISTORY'] = message


def fits_write(hdulist, outfile, comment=None, history=None, overwrite=False):
    """
    Write the HDUList to file.

    Control some actions to be taken for all files written by nuskybgd.
    """
    hdulist.writeto(outfile, checksum=True, overwrite=overwrite)


def fits_checkkeyword(fitsfile, keyword, ext=0, silent=False):
    """
    Check the keyword value of a FITS extension.

    Inputs:

    fitsfile -- Path to the FITS file.
    keyword -- The keyword to check.
    ext -- Extension index (int) or key (str)

    Return:

    Value of the keyword.

    Both the specified extension and keyword exist. If either condition is not
    met, a KeyError exception will be raised.

    If silent=True, return None without raising KeyError, unless the
    specified file cannot be found, in which case astropy.io.fits will raise
    an OSError.
    """
    import astropy.io.fits as pf

    fh = pf.open(fitsfile)
    try:
        return fh[ext].header[keyword]
    except KeyError as e:
        if silent:
            return None
        else:
            print('The specified extension or keyword is not found.')
            raise e


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

    # Store references to the columns (numpy arrays), order of magnitude
    # faster for loop operation compared to querying gtiarr and aimpointarr
    # rows (FITS record objects).
    gtistart = gtiarr['START']
    gtistop = gtiarr['STOP']
    aimtimes = aimpointarr['TIME']
    flagarr = np.zeros(n_aimpoint, dtype=np.int8)

    filter_gti_loop(aimtimes, gtistart, gtistop, flagarr)

    flagged = np.where(flagarr == 1)[0]
    view1 = aimpointarr[flagged]
    view2 = aimpointarr[flagged + 1]
    coords = np.array([view1['X_DET1'], view1['Y_DET1']]).T
    dt = view2['TIME'] - view1['TIME']

    print('%d aimpoints (%f / %f s)' % (
        len(coords), np.sum(dt),
        np.sum(gtiarr.field('STOP') - gtiarr.field('START'))))

    return coords, dt


def filter_gti_loop(times, gtistart, gtistop, flags):
    """
    Logic for filtering aimpoints by GTI intervals. Modify flags array
    in-place.
    """
    n_aimpoint = len(times)
    i_gti = 0

    print('Processing %d aimpoints...' % n_aimpoint)

    for i_aimpoint in range(n_aimpoint - 1):
        # Ref pointings have a single time while GTI intervals have two.
        # Original logic in projobs.pro considers pointing time against
        # [start, stop), i.e. >= for start time, < for stop.

        aim_time = times[i_aimpoint]

        # Ensure current GTI interval ends after current pointing time
        try:
            while aim_time >= gtistop[i_gti]:
                i_gti += 1
        except IndexError:
            # Reached the end of the GTI list
            break

        # Flag this pointing if it is inside current GTI interval
        flags[i_aimpoint] = int(aim_time >= gtistart[i_gti])


def fpm_parse(keyword):
    if keyword not in ('FPMA', 'FPMB'):
        return False
    else:
        return keyword[-1]


def mask_from_region(regfile, refimg):
    """
    Create a pixel mask for regfile based on the image WCS info in refimg.

    Uses the pyregion module. Please keep to using circle, box, and ellipse
    shapes in DS9, fk5 format to avoid unexpected behavior.
    """
    import astropy.io.fits as pf

    fh = pf.open(refimg)

    if isinstance(regfile, str):
        reg = pyregion.open(regfile)
        return reg.get_mask(hdu=fh[0])

    elif isinstance(regfile, list):
        masks = []
        for f in regfile:
            reg = pyregion.open(f)
            masks.append(reg.get_mask(hdu=fh[0]))
        return masks

    else:
        raise ValueError('regfile should be a filename or list of filenames.')
