# (c) 2019 Qian Wang
# This file is part of nuskybgd released under MIT License.
# See LICENSE file for full details.

import datetime
import pyregion
import numpy as np


def format_timestamp():
    return datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%S')


def hdu_timestamp_write(hdu):
    """
    Update the DATE keyword of the input HDU header to current time. This
    modifies the input object 'hdu'.

    Parameters
    ----------
    hdu : astropy.io.fits.Header
        FITS header object to update.
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

    Parameters
    ----------
    hdu : astropy.io.fits.Header
        FITS header object to update.
    message : str
        Text to add to the HISTORY key.
    """
    hdu.header['HISTORY'] = message


def fits_write(hdulist, outfile, comment=None, history=None, overwrite=False):
    """
    Write the HDUList to file.

    This controls the actions to be taken for all files written by nuskybgd.
    Currently, this includes the option ``checksum=True``.

    Parameters
    ----------
    hdulist : astropy.io.fits.HDUList
        The FITS file object to write.

    outfile : str
        The output file path.

    comment
        (Optional) Currently unused.

    history
        (Optional) Currently unused.

    overwrite : bool
        (Optional, default: False) Whether to overwrite existing output file,
        this is passed on to HDUList.writeto().
    """
    hdulist.writeto(outfile, checksum=True, overwrite=overwrite)


def fits_checkkeyword(fitsfile, keyword, ext=0, silent=False):
    """
    Check the keyword value of a FITS extension.

    Parameters
    ----------
    fitsfile : str
        Path to the FITS file.
    keyword : str
        The keyword to check.
    ext : int or str
        Extension index (int) or key (str).

    Returns
    -------
    Header key value
        If both the specified extension and keyword exist.

    ``None``
        If a ``KeyError`` exception would have been raised and ``silent=True``
        is set.

    Raises
    ------
    KeyError
        If either the specified extension or the keyword cannot be found, and
        ``silent=False``, a KeyError exception will be raised.
    OSError
        If the specified file cannot be found, astropy.io.fits will raise
        OSError.
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
    Check a FITS header to see if it has the expected data table columns for
    telescope aspect.

    Require the following keywords::

        TELESCOP = 'NuSTAR'
        INSTRUME = 'FPM?'           # ? is 'A' or 'B'
        EXTNAME  = 'DET?_REFPOINT'  # ? is a digit 1 to 4
        HDUCLAS1 = 'TEMPORALDATA'
        NAXIS    = 2

    Require the following fields::

        TTYPE1   = 'TIME'
        TUNIT1   = 's'
        TTYPE2   = 'X_DET1'
        TUNIT2   = 'pixel'
        TTYPE3   = 'Y_DET1'
        TUNIT3   = 'pixel'

    Parameters
    ----------
    hdr : astropy.io.fits.Header
        FITS header to check.

    Returns
    -------
    bool
        Whether the header passes the check.
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

    Require the following keywords::

        TELESCOP = 'NuSTAR'        / Telescope (mission) name
        HDUCLAS1  = 'GTI'          / File contains Good Time Intervals
        NAXIS    = 2               / 2-dimensional binary table

    Require the following fields::

        TTYPE1   = 'START'         / label for field   1
        TUNIT1   = 's' or 'sec'    / physical unit of field
        TTYPE2   = 'STOP'          / label for field   2
        TUNIT2   = 's' or 'sec'    / physical unit of field

    Parameters
    ----------
    hdr : astropy.io.fits.Header
        FITS header to check.

    Returns
    -------
    bool
        Whether the header passes the check.
    """
    try:
        if not (
            hdr['TELESCOP'] == 'NuSTAR' and
            hdr['HDUCLAS1'] == 'GTI' and
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
    """
    Filter an aimpoints data table using good time intervals from a GTI file.

    Parameters
    ----------
    aimpointext : astropy.io.fits.BinTableHDU
        The aimpoints extension.
    gtiext : astropy.io.fits.BinTableHDU
        The GTI extension.

    Returns
    -------
    (numpy.ndarray, np.ndarray)
        A numpy array of [detector x, detector y] coordinates of filtered
        aimpoints, and a corresponding numpy array listing the duration at
        that aimpoint.
    """
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
    """
    Parse the instrument name keyword to determine whether it is FPMA or FPMB.

    Parameters
    ----------
    keyword : str
        The string to check

    Returns
    -------
    str, 'A' or 'B'
    False
        If it is neither.
    """
    if keyword not in ('FPMA', 'FPMB'):
        return False
    else:
        return keyword[-1]


def mask_from_region(regfile, refimg):
    """
    Create a pixel mask for regfile based on the image WCS info in refimg.

    Uses the pyregion module. Please keep to using circle, box, and ellipse
    shapes in DS9, fk5 format to avoid unexpected behavior.

    Parameters
    ----------
    regfile : str, or list of str
        Path to DS9 region file, or a list of paths.
    refimg : str
        Path to the reference image file.

    Returns
    -------
    numpy.ndarray
        The mask image, or an array of mask images if the input ``regfile``
        was a list.

    Raises
    ------
    ValueError
        If the input ``regfile`` is neither a string or a list.
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


def print_hr(width=78):
    """
    Print a horizontal line made from hyphens.

    Parameters
    ----------
    width : int
        Width of the line in characters, must be > 0.
    """
    if width > 1:
        print('-' * width)


def docformat(txt):
    """Format some docstring for display in terminal."""
    return txt.format(b='\033[1;31m', o='\033[0;0m')
