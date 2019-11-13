# (c) 2019 Qian Wang
# This file is part of nuskybgd released under MIT License.
# See LICENSE file for full details.

import astropy.io.fits as pf
import os
import numpy as np
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


def detabs_hdu(abshdus, weights):
    """
    Coadd absorption spectrum of the detectors using the given weights. The
    weights are normalized to unity before proceeding.

    Example uses:

    - Get detector absorption weighed by area inside a source region.
    - Get absorption of a single detector by giving it 100% weight.

    Inputs:

    abshdus -- HDUList (e.g. CALDB DETABS FITS file)
    weights -- A list of 4 values

    Outputs:

    BinTableHDU 'DETABS' containing absorption data
    """
    if not ((isinstance(weights, list) or
             isinstance(weights, np.ndarray)) and
            len(weights) == 4):
        raise ValueError('weights must be a length-4 array')
    if not isinstance(abshdus, pf.hdu.hdulist.HDUList):
        raise ValueError('detabshdus must be a HDUList object')

    weight_sum = np.sum(weights)

    if not (weight_sum > 0 and np.min(weights) >= 0):
        raise ValueError('weights must be >= 0 and sum to > 0')

    weights = np.float64(weights) / weight_sum

    # Create new HDU for the results
    results = abshdus[1].copy()
    results.data['DETABS'] *= 0.0

    for idet in range(4):
        if weights[idet] > 0:
            results.data['DETABS'] += (
                weights[idet] * abshdus[idet + 1].data['DETABS'])

    util.hdu_timestamp_write(results)
    msg = 'Weighted DETABS: %.8e %.8e %.8e %.8e' % tuple(weights)
    util.hdu_history_write(results, msg)

    return results


def add_weighted_abs(abshdus, weights, rmfhdu):
    """
    Add absorption to an RMF, using given weights for the detectors.

    **The input 'rmfhdu' is modified in-place.**

    Inputs:

    abshdus -- HDUList (e.g. CALDB DETABS FITS file)
    weights -- A list of 4 values
    rmfhdu -- BinTableHDU with RMF data
    """
    weighted_abs = detabs_hdu(abshdus, weights)

    # Modifies rmffh['MATRIX'] in place
    add_abs_hdu(weighted_abs, rmfhdu)

    msg = (
        'RMF includes DETABS weights: '
        '%.8e %.8e %.8e %.8e' % tuple(weights)
    )
    util.hdu_history_write(rmfhdu, msg)


def make_absrmf(evtfile, outfile,
                rmffile='CALDB', detabsfile='CALDB',
                overwrite=False):

    evtfh = pf.open(evtfile)
    evthdr = evtfh['EVENTS'].header
    ab = evthdr['INSTRUME'][-1]
    cal = caldb.CalDB(os.environ['CALDB'])

    outfilename = outfile + '%d%s.rmf'

    # Check if output files exist
    dirname = os.path.dirname(outfile)
    if dirname != '':
        os.makedirs(dirname, exist_ok=True)
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
        weights = [0.0] * 4
        weights[idet] = 1.0

        # Modifies rmffh['MATRIX'] in place
        add_weighted_abs(detabsfh, weights, rmffh['MATRIX'])

        util.fits_write(rmffh, outfilename % (idet, ab),
                                 overwrite=overwrite)

    return True
