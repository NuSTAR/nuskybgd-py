#!/usr/bin/env python3
"""
Create RMF files that includes detector absorption (DETABS).

Usage:

absrmf.py evtfile outfile [rmffile=CALDB] [detabsfile=CALDB]

evtfile is an event file from which the INSTRUME and DATE-OBS keywords are
used.

outfile will be prefixed to the output file names, and can be a file path.

rmffile is the RMF file to multiply by absorption, set it to CALDB (default)
to use the latest CALDB file(s).

detabsfile is the detector absorption file to multiply the RMF with, set it to
CALDB (default) to use the latest CALDB file(s).
"""

import astropy.io.fits as pf
import os
import datetime
import sys

import nuskybgd
import nuskybgd.caldb


def make_absrmf(evtfile, outfile,
                rmffile='CALDB', detabsfile='CALDB',
                overwrite=False):
    evtfh = pf.open(evtfile)
    evthdr = evtfh['EVENTS'].header
    ab = evthdr['INSTRUME'][-1]
    cal = nuskybgd.caldb.CalDB(os.environ['CALDB'])

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
            grprmffile = '%s/%s' % (
                os.environ['CALDB'],
                cal.getGRPRMF(evthdr['INSTRUME'], detnam, evthdr['DATE-OBS'])
            )
            grprmffh = pf.open(grprmffile)
            if len(grprmffh[idet + 1].data) > 1:
                print('Warning: more than one file for GRPRMF %s found. '
                      'Using first entry.' % detnam)
            rmflist.append('%s/%s' % (
                os.path.dirname(grprmffile),
                grprmffh[idet + 1].data[0].field('RMFFILE')
            ))
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

        absmatrix = (rmffh['MATRIX'].data.field('MATRIX') *
                     detabsfh[idet + 1].data.field('DETABS'))

        # Create the output file
        absrmffh = pf.open(rmflist[0])
        timestamp = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%S')
        absrmffh[0].header['DATE'] = (
            timestamp,
            'File creation date (YYYY-MM-DDThh:mm:ss UTC)'
        )
        absrmffh[1].header['DATE'] = (
            timestamp,
            'File creation date (YYYY-MM-DDThh:mm:ss UTC)'
        )
        absrmffh[1].header['HISTORY'] = (
            'RMF multiplied by DETABS using files %s and %s' % (
                rmffile, detabsfile
            )
        )
        absrmffh['MATRIX'].data['MATRIX'] = absmatrix
        absrmffh.writeto(outfilename % (idet, ab),
                         checksum=True, overwrite=overwrite)

    return True


if __name__ == '__main__':

    if len(sys.argv) not in (3, 4, 5):
        print(__doc__)
        sys.exit(0)

    evtfile = sys.argv[1]
    outfile = sys.argv[2]

    keywords = {
        'rmffile': 'CALDB',
        'detabsfile': 'CALDB'
    }

    for _ in sys.argv[3:]:
        arg = _.split('=')
        if arg[0] in keywords:
            keywords[arg[0]] = arg[1]

    # File exist checks
    halt = False
    if not os.path.exists(evtfile):
        print('%s not found.' % evtfile)
        halt = True
    if (keywords['rmffile'] != 'CALDB' and
            not os.path.exists(keywords['rmffile'])):
        print('%s not found.' % keywords['rmffile'])
        halt = True
    if (keywords['detabsfile'] != 'CALDB' and
            not os.path.exists(keywords['detabsfile'])):
        print('%s not found.' % keywords['detabsfile'])
        halt = True

    if halt:
        sys.exit(1)

    # Overwrite output flag?
    overwrite = False
    if outfile[0] == '!':
        overwrite = True
        outfile = outfile[1:]

    make_absrmf(evtfile, outfile,
                rmffile=keywords['rmffile'],
                detabsfile=keywords['detabsfile'],
                overwrite=overwrite)

    sys.exit(0)
