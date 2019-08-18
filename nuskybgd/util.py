"""
Utility functions for sharing.
"""
import datetime

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
