"""
Utility functions for sharing.
"""
import datetime
import astropy.io.fits as pf
import pyregion


def hdu_timestamp_write(hdu):
    """
    Update the DATE keyword of the input HDU header to current time. This
    modifies the input object 'hdu'.
    """
    import datetime

    timestamp = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%S')
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
    fh = pf.open(fitsfile)
    try:
        return fh[ext].header[keyword]
    except KeyError as e:
        if silent:
            return None
        else:
            print('The specified extension or keyword is not found.')
            raise e


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
    fh = pf.open(refimg)
    reg = pyregion.open(regfile)
    mask = reg.get_mask(hdu=fh[0])
    return mask
