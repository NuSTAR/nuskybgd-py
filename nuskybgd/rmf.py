import astropy.io.fits as pf
import os
from . import caldb
from . import util

def add_abs_hdu(detabshdu, rmfhdu):
    """
    Add detector absorption to response matrix. Doing so may obviate the need
    for a separate ARF file. This function works with FITS HDUs.

    **The input 'rmfhdu' is modified.**

    Example uses:

    - Add DETABS from CALDB to the RMF file from CALDB.
    - Add weighed DETABS to a user supplied RMF file.

    Inputs:

    detabshdu -- HDU with BinTable 'DETABS' containing absorption data
    rmfhdu -- HDU with BinTable 'MATRIX' containing RMF data
    """
    if not (isinstance(detabshdu, pf.hdu.table.BinTableHDU) and
            isinstance(rmfhdu, pf.hdu.table.BinTableHDU)):
        raise TypeError('Inputs must be BinTable HDUs.')

    try:
        absmatrix = (rmfhdu.data['MATRIX'] *
                         detabshdu.data['DETABS'])
    except KeyError:
        raise KeyError('MATRIX and DETABS columns are required.')

    rmfhdu.data['MATRIX'] = absmatrix
    util.hdu_timestamp_write(rmfhdu)


def add_abs_fits(detabsfile, rmffile):
    """
    Add detector absorption to response matrix. Doing so may obviate the need
    for a separate ARF file. This function works with FITS files.
    """
    pass



def make_absrmf(evtfile, outfile,
                rmffile='CALDB', detabsfile='CALDB',
                overwrite=False):

    evtfh = pf.open(evtfile)
    evthdr = evtfh['EVENTS'].header
    ab = evthdr['INSTRUME'][-1]
    cal = caldb.CalDB(os.environ['CALDB'])

    outfilename = outfile + '%d%s.rmf'
    # Check if output files exist
    halt = False
    for idet in range(4):
        if os.path.exists(outfilename % (idet, ab)):
            print('File exists: %s' % outfilename % (idet, ab))
            halt = True
    if halt:
        if not overwrite:
            return False
        else:
            print('Overwriting existing files...')

    # Looks up rmf files and detabs files
    rmflist = []
    if rmffile == 'CALDB':
        for idet in range(4):
            detnam = 'DET%d' % idet
            rmflist.append(cal.getRMF(evthdr['INSTRUME'],
                                      detnam, evthdr['DATE-OBS']))
    else:
        rmflist = [rmffile] * 4

    detabslist = []
    if detabsfile == 'CALDB':
        for idet in range(4):
            detnam = 'DET%d' % idet
            detabslist.append('%s/%s' % (
                os.environ['CALDB'],
                cal.getDETABS(evthdr['INSTRUME'], detnam, evthdr['DATE-OBS'])
            ))
    else:
        detabslist = [detabsfile] * 4

    # Print the files used
    print('Calibration files used:')
    for _ in zip([0, 1, 2, 3], rmflist, detabslist):
        print('DET%d %s %s' % (_[0], _[1], _[2]))

    for idet in range(4):
        rmffh = pf.open(rmflist[idet])
        detabsfh = pf.open(detabslist[idet])

        # Modifies rmffh['MATRIX'] in place
        add_abs_hdu(detabsfh[idet + 1], rmffh['MATRIX'])

        # Create the output file
        history_msg = (
            'Modified RMF which includes DETABS, using files %s and %s' % (
                rmflist[idet], detabslist[idet]
            )
        )
        util.hdu_history_write(rmffh['MATRIX'], history_msg)
        util.fits_write(rmffh, outfilename % (idet, ab),
                                 overwrite=overwrite)

    return True
