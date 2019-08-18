#!/usr/bin/env python3
"""
instrmap.py

Usage:

instrmap.py A01_cl.evt [usrbpix=morebadpix.fits] [prefix=prefix] [dryrun=yes]

Create an instrument map for a specific observation. The user can specify
additional bad pixels by giving files that have BADPIX extension(s) using the
usrbpix= argument. Multiple files can be given, separated by commas (with no
space in the argument). Bad pixels from CALDB will always be applied.

The prefix= argument lets user add a file name prefix for the output file.
This allows the output to be written to any specified path. If a directory is
intended, it must end with '/'.

If dryrun=yes the task will stop after the bad pixel files have been checked,
before the instrument map is calculated.
"""
import sys
import os
import datetime
import astropy.io.fits as pf
import numpy as np
import scipy.ndimage
import nuskybgd
import nuskybgd.caldb


class CaldbIndx():
    def __init__(self, caldbindxpath):
        self.indxpath = caldbindxpath
        self.indxfits = pf.open(caldbindxpath)
        self.version = self.indxfits['CIF'].header['CALDBVER']


def fpm_parse(keyword):
    if keyword not in ('FPMA', 'FPMB'):
        return False
    else:
        return keyword[-1]


DETNAM = {'DET0': 0, 'DET1': 1, 'DET2': 2, 'DET3': 3}

# grade weighting from NGC253 002 obs.
GRADE_WT = [1.00000, 0.124902, 0.117130, 0.114720, 0.118038,
            0.0114296, 0.0101738, 0.0113617, 0.0122017, 0.0157910,
            0.0144079, 0.0145691, 0.0149934, 0.00165462, 0.00194312,
            0.00156128, 0.00143400, 0.00210433, 0.00180735, 0.00140006,
            0.00169704, 0.00189220, 0.00160371, 0.00150188, 0.00168007,
            0.000296983, 0.000364864]


def shift_image(image, delta):
    return scipy.ndimage.shift(image, delta,
                               mode='wrap', prefilter=False, order=1)


def nuskybgd_timestamp():
    return ('%s nuskybgd' %
            datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))


def get_caldb_instrmap(evthdr):
    """
    Get instrument map image from CALDB, shifted by (-1, -1).

    Returns (image, FITS header).
    """
    fpm = evthdr['INSTRUME']
    obsutctime = evthdr['DATE-OBS']

    caldb = nuskybgd.caldb.CalDB(os.environ['CALDB'])
    imappath = caldb.getINSTRMAP(fpm, obsutctime)
    imapf = pf.open('%s/%s' % (os.environ['CALDB'], imappath))
    try:
        imap = imapf['INSTRMAP'].data
    except KeyError:
        print('Error: INSTRMAP extension not found.')
        return False
    hdr = imapf['INSTRMAP'].header
    hdr['HISTORY'] = nuskybgd_timestamp()
    hdr['HISTORY'] = 'Created inst. map from %s' % os.path.basename(imappath)
    return imap, hdr


def get_caldb_pixpos(evthdr):
    fpm = evthdr['INSTRUME']
    obsutctime = evthdr['DATE-OBS']
    caldb = nuskybgd.caldb.CalDB(os.environ['CALDB'])
    pixpospath = caldb.getPIXPOS(fpm, 'DET0', obsutctime)

    pixposf = pf.open('%s/%s' % (os.environ['CALDB'], pixpospath))

    pixmap = np.full((360, 360), -1, dtype=np.int32)
    detnum = np.full((360, 360), -1, dtype=np.int32)
    allpdf = np.zeros((360, 360), dtype=np.float64)

    for ext in pixposf:
        if (('EXTNAME' not in ext.header) or
            ('PIXPOS' not in ext.header['EXTNAME']) or
                ('DETNAM' not in ext.header)):
            continue
        idet = int(ext.header['DETNAM'].replace('DET', ''))
        pixpos = ext.data
        for ix in np.arange(32):
            for iy in np.arange(32):
                # Get array indices where all of the following are True
                ii = np.where((pixpos.field('REF_DET1X') != -1) *
                              (pixpos.field('RAWX') == ix) *
                              (pixpos.field('RAWY') == iy) *
                              (pixpos.field('GRADE') <= 26))[0]

                thispdf = np.zeros((360, 360), dtype=np.float64)

                for i in ii:
                    if not np.isnan(pixpos.field('PDF')[i]).any():
                        # No nan value in PDF
                        ref_x = pixpos.field('REF_DET1X')[i]
                        ref_y = pixpos.field('REF_DET1Y')[i]
                        thispdf[ref_y:ref_y + 7, ref_x:ref_x + 7] += (
                            pixpos.field('PDF')[i] *
                            GRADE_WT[pixpos.field('GRADE')[i]])

                ii = np.where(thispdf > allpdf)
                if len(ii) > 0:
                    allpdf[ii] = thispdf[ii]
                    pixmap[ii] = ix + iy * 32
                    detnum[ii] = idet

    pixmap = shift_image(pixmap, [-1, -1])
    detnum = shift_image(detnum, [-1, -1])

    return pixmap, detnum


def get_badpix_exts(bpixfiles):
    # Collect list of bad pixel extensions, group into FPMA and B
    bpixlist = {'A': [], 'B': []}
    for _ in bpixfiles:
        ff = pf.open(_)
        # Check FPM
        phdr = ff[0].header
        try:
            if phdr['TELESCOP'] != 'NuSTAR':
                print('%s: TELESCOP != NuSTAR, skipping.' % _)
                continue
        except KeyError:
            print('%s: no TELESCOP keyword, skipping.' % _)
            continue

        try:
            if phdr['INSTRUME'] == 'FPMA':
                ab = 'A'
            elif phdr['INSTRUME'] == 'FPMB':
                ab = 'B'
            else:
                print('%s: unknown INSTRUME: %s' % (_, phdr['INSTRUME']))
                continue
        except KeyError:
            print('%s: no INSTRUME keyword, skipping.' % _)
            continue

        for ext in ff[1:]:
            try:
                ehdr = ext.header
                if not (ehdr['XTENSION'] == 'BINTABLE' and
                        'BADPIX' in ehdr['EXTNAME']):
                    continue
                if ehdr['DETNAM'] in DETNAM:
                    bpixlist[ab].append(ext)
                    print('Added %s %s %s to list.' % (_, ab, ehdr['DETNAM']))
                else:
                    print('%s: unknown DETNAM: %s' % (_, ehdr['DETNAM']))
                    continue
            except KeyError:
                # Extension does not have XTENSION, EXTNAME or DETNAM
                continue

    return bpixlist


def apply_badpix(instrmap, bpixexts, pixmap, detnum):
    output = 1. * instrmap

    for ext in bpixexts:
        hh = ext.header
        if ('TSTART' in hh and 'TSTOP' in hh):
            tobs = hh['TSTOP'] - hh['TSTART']
        else:
            tobs = None

        if not ('EXTNAME' in ext.header and
                'BADPIX' in ext.header['EXTNAME'] and
                'DETNAM' in ext.header and
                ext.header['DETNAM'] in DETNAM):
            continue

        idet = DETNAM[ext.header['DETNAM']]
        print(idet)
        badpix = ext.data

        x = badpix.field('RAWX')
        y = badpix.field('RAWY')
        flags = badpix.field('BADFLAG')
        if tobs is not None:
            dt = badpix.field('TIME_STOP') - badpix.field('TIME')

        for i in np.arange(len(x)):
            if ((tobs is not None and dt[i] > 0.8 * tobs) or
                    flags[i][-2] is True):
                ii = np.where((pixmap == x[i] + y[i] * 32) *
                              (detnum == idet))
                output[ii] = 0

            ii = np.where(output > 0)
            output[ii] = detnum[ii] + 1

    return output


if __name__ == '__main__':

    ERR_CALDB_NOTFOUND = """
Error: checking for calibration files in $CALDB/data/nustar, but cannot find
caldb.indx in this CALDB location. Is the $CALDB path correct?

$CALDB = {caldbpath}

Cannot find the file

{caldbinxpath}
"""

    ERR_BPIXFILE_NOTFOUND = """
Error: bad pixel file not found, skipping:
{file}
"""

    if len(sys.argv) == 1:
        print(__doc__)
        sys.exit(0)

    obsinfofile = sys.argv[1]

    keywords = {
        'usrbpix': None,
        'prefix': '',
        'dryrun': ''
    }

    for _ in sys.argv[1:]:
        arg = _.split('=')
        if arg[0] in keywords:
            keywords[arg[0]] = arg[1]

    scriptname = os.path.basename(__file__)
    print('Running %s. Performing checks...' % scriptname)

    halt = False

    # Inputs File with observation info must exist. Typically this is the
    # event file itself, but this program needs only the information in its
    # FITS header.
    require_hdr_keywords = ['INSTRUME']
    if not os.path.exists(obsinfofile):
        print('Error: %s does not exist.' % obsinfofile)
        halt = True
        sys.exit(1)  # No further checks if the main file doesn't exist

    evtfile = pf.open(obsinfofile)
    if 'EVENTS' not in evtfile:
        print('Error: EVENTS extension not found in %s.' % obsinfofile)
        halt = True
        sys.exit(1)

    evthdr = evtfile['EVENTS'].header
    for _ in require_hdr_keywords:
        if _ not in evthdr:
            halt = True
            print('Error: required keyword %s not in %s' % (_, obsinfofile))

    ab = fpm_parse(evthdr['INSTRUME'])
    if ab is False:
        halt = True
        print('Invalid INSTRUME keyword value: %s' % evthdr['INSTRUME'])

    if halt:
        sys.exit(1)

    # FITS file(s) with BADPIX extensions
    bpixfiles = [obsinfofile]

    if keywords['usrbpix'] is not None:
        for _ in keywords['usrbpix'].split(','):
            if os.path.exists(_):
                bpixfiles.append(_)
            else:
                print(ERR_BPIXFILE_NOTFOUND.format(file=_))

    # Output will be named outprefix_A.fits and/or outprefix_B.fits
    outprefix = keywords['prefix']

    outfilename = outprefix + 'newinstrmap' + ab + '.fits'
    if os.path.exists(outfilename):
        print('Error: output file %s exists!' % outfilename)
        halt = True

    # Check CALDB environment
    try:
        caldb = os.environ['CALDB']
    except KeyError:
        halt = True
        print('Error: environment variable CALDB is undefined.')
    else:
        # Check if it is valid by looking for a file in nustar's fpm data.
        # Doesn't guarantee the other files are OK, but will warn if CALDB
        # path is incorrect altogether.
        caldbinxpath = '$CALDB/data/nustar/fpm/caldb.indx'
        if not os.path.exists(caldbinxpath.replace('$CALDB', caldb)):
            halt = True
            print(ERR_CALDB_NOTFOUND.format(
                caldbpath=caldb, caldbinxpath=caldbinxpath))
        else:
            print('$CALDB=%s' % caldb)
            caldbinx = CaldbIndx(caldbinxpath.replace('$CALDB', caldb))
            print('CALDB version: %s' % caldbinx.version)

    if halt:
        print('%s did not complete. See error messages.' % scriptname)
        sys.exit(1)

    print('Preliminary checks complete, proceeding...')

    print('Creating instrument maps for FPM%s for observation %s' % (
        ab, obsinfofile))

    caldb = nuskybgd.caldb.CalDB(os.environ['CALDB'])
    caldbbpixpath = caldb.getBADPIX(
        evthdr['INSTRUME'], 'DET0', evthdr['DATE-OBS'])
    bpixfiles.append('%s/%s' % (os.environ['CALDB'], caldbbpixpath))

    print('Collecting bad pixel lists...')
    bpixexts = get_badpix_exts(bpixfiles)

    # Stop here if doing a dry run
    if keywords['dryrun'] == 'yes':
        sys.exit(0)

    print('Loading CALDB instrmap...')
    instrmap, header = get_caldb_instrmap(evthdr)
    print('Loading CALDB pixpos...')
    pixmap, detnum = get_caldb_pixpos(evthdr)

    print('Applying bad pixels lists for FPM%s...' % ab)
    masked_instrmap = apply_badpix(instrmap, bpixexts[ab], pixmap, detnum)

    pf.HDUList(
        pf.PrimaryHDU(masked_instrmap, header=header)
    ).writeto(outfilename)


"""
Bad pixels.
https://heasarc.gsfc.nasa.gov/ftools/caldb/help/nuflagbad.html

The list of bad pixels is stored in four distinct new extensions, named
'BADPIX', one for each of the four detectors. Optionally, if requested by the
user through the parameter 'outbpfile', the bad pixels list is also written in
a separate output file.

The bad pixel extensions contain the positional and temporal information and a
16-bit binary number column, named 'BADFLAG', indicating the class of the bad
pixel with the following meaning:

Bad pixels BADFLAG flags
b0000000000000001 Bad pixel from on-ground CALDB Bad Pixel File
b0000000000000010 Disabled pixel from on-board software
b0000000000000100 Bad pixel in the file provided by the user


CALDB file content notes.

Filename: nuApixpos20100101v007.fits
No.    Name      Ver    Type      Cards   Dimensions   Format
  0  PRIMARY       1 PrimaryHDU      12   ()
  1  PIXPOS        1 BinTableHDU     66   32768R x 6C   [1B, 1B, 1I, 1I, 1I, 49E]
  2  PIXPOS        2 BinTableHDU     66   32768R x 6C   [1B, 1B, 1I, 1I, 1I, 49E]
  3  PIXPOS        3 BinTableHDU     66   32768R x 6C   [1B, 1B, 1I, 1I, 1I, 49E]
  4  PIXPOS        4 BinTableHDU     66   32768R x 6C   [1B, 1B, 1I, 1I, 1I, 49E]

PRIMARY header:

SIMPLE  =                    T / file does conform to FITS standard
BITPIX  =                   16 / number of bits per data pixel
NAXIS   =                    0 / number of data axes
EXTEND  =                    T / FITS dataset may contain extensions
COMMENT   FITS (Flexible Image Transport System) format is defined in 'Astronomy
COMMENT   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H
TELESCOP= 'NuSTAR   '          / Telescope (mission) name
INSTRUME= 'FPMA    '           / Instrument name (FPMA or FPMB)
ORIGIN  = 'Caltech   '         / Source of FITS file
CREATOR = 'FTOOLS 6.13 '       / Creator
CHECKSUM= 'Bdh3EZe2Bde2BZe2'   / HDU checksum updated 2013-07-25T17:08:37
DATASUM = '         0'         / data unit checksum updated 2013-07-25T17:08:37

PIXPOS header:

XTENSION= 'BINTABLE'           / binary table extension
BITPIX  =                    8 / 8-bit bytes
NAXIS   =                    2 / 2-dimensional binary table
NAXIS1  =                  204 / width of table in bytes
NAXIS2  =                32768 / number of rows in table
PCOUNT  =                    0 / size of special data area
GCOUNT  =                    1 / one data group (required keyword)
TFIELDS =                    6 / number of fields in each row
TTYPE1  = 'RAWX    '           / X-position in Raw Detector coordinates
TFORM1  = '1B      '           / data format of field: BYTE
TUNIT1  = 'pixel   '           / physical unit of field
TTYPE2  = 'RAWY    '           / Y-position in Raw Detector coordinates
TFORM2  = '1B      '           / data format of field: BYTE
TUNIT2  = 'pixel   '           / physical unit of field
TTYPE3  = 'GRADE   '           / Event Grade
TFORM3  = '1I      '           / data format of field: 2-byte INTEGER
TTYPE4  = 'REF_DET1X'          / Focal Plane X value of lower-left corner of PDF
TFORM4  = '1I      '           / data format of field: 2-byte INTEGER
TUNIT4  = 'pixel   '           / physical unit of field
TTYPE5  = 'REF_DET1Y'          / Focal Plane Y value of lower-left corner of PDF
TFORM5  = '1I      '           / data format of field: 2-byte INTEGER
TUNIT5  = 'pixel   '           / physical unit of field
TTYPE6  = 'PDF     '           / Probability Distribution Function Map
TFORM6  = '49E     '           / data format of field: 4-byte REAL
TDIM6   = '(7, 7)  '           / Array dimensions for column 6
EXTNAME = 'PIXPOS  '           / Name of the binary table extension
EXTVER  =                    1 / There shall be one instance of this extension f
DETNAM  = 'DET0    '           / CZT Detector ID (0,1,2 or 3)
TELESCOP= 'NuSTAR  '           / Telescope (mission) name
INSTRUME= 'FPMA    '           / Instrument name (FPMA or FPMB)
ORIGIN  = 'Caltech   '         / Source of FITS file
CREATOR = 'FTOOLS 6.13 '       / Creator
VERSION =                    7 / Extension version number
FILENAME= 'nuApixpos20100101v007.fits' / File name
CONTENT = 'NuSTAR Pixel Positions' / File content
TIMESYS = 'TT      '           / Terrestrial Time: synchronous with, but 32.184
MJDREFI =                55197 / MJD reference day 01 Jan 2010 00:00:00 UTC
MJDREFF =        0.00076601852 / MJD reference (fraction part: 32.184 secs + 34
COMMENT MJDREFI + MJDREFF is the epoch January 1.0, 2010, in the TT time system.
CLOCKAPP=                    F / Set to TRUE if correction has been applied to t
CCLS0001= 'BCF     '           / Daset is a Basic Calibration File
CDTP0001= 'DATA    '           / Calibration file contains data
CCNM0001= 'PIXPOS  '           / Type of calibration data
CVSD0001= '2010-01-01'         / Date when this file should first be used
CVST0001= '00:00:00'           / Time of day when this file should first be used
CDES0001= 'NuSTAR Pixel Position Coefficients' / Description
COMMENT
COMMENT   This extension provides, for each pixel of the detector, the probabili
COMMENT   distribution for events with Grades 0 through 13 obtained from X-ray
COMMENT   generator scans during on-ground calibratiions.
COMMENT   For each pixel and grade, a 7x7 grid is stored in the CALDB. These bin
COMMENT   represent the Focal Plane Bench Frame (DET1X) bins. The REF_DET1X and
COMMENT   REF_DET1Y columns give the address of the lower-left corner of the 7x7
COMMENT   grid in DET1X/Y integer units (DET1X/Y runs from 1 to 360).
COMMENT   The first 7 values represent the probability that a pixel occurred in
COMMENT   pixels (REF_DET1X, REF_DET1Y), (REF_DET1X+1, REF_DET1Y) ...
COMMENT   (REF_DET1X+6, REF_DET1Y). The second seven are (REF_DET1X, REF_DET1Y+1
COMMENT   through (REF_DET1X+6, REF_DET1Y+1) etc.
COMMENT
COMMENT  Nominal pixel size (microns): PIXSIZE = 604.80
COMMENT  Nominal Software Pixel Size (microns): CELLSIZE = 120.96
COMMENT  CELLSIZE is defined to be PIXSIZE / 5
COMMENT
DATE    = '2013-07-25T17:08:08' / file creation date (YYYY-MM-DDThh:mm:ss UT)
CHECKSUM= 'gq9Rjo7Pgo7Pgo7P'   / HDU checksum updated 2013-07-25T17:08:37
DATASUM = '1600104595'         / data unit checksum updated 2013-07-25T17:08:37


PIXPOS columns:

ColDefs(
    name = 'RAWX'; format = '1B'; unit = 'pixel'
    name = 'RAWY'; format = '1B'; unit = 'pixel'
    name = 'GRADE'; format = '1I'
    name = 'REF_DET1X'; format = '1I'; unit = 'pixel'
    name = 'REF_DET1Y'; format = '1I'; unit = 'pixel'
    name = 'PDF'; format = '49E'; dim = '(7, 7)'
)

dtype=(numpy.record, [('RAWX', 'u1'), ('RAWY', 'u1'), ('GRADE', '>i2'),
('REF_DET1X', '>i2'), ('REF_DET1Y', '>i2'), ('PDF', '>f4', (7, 7))])


*_bp.fits

--------------------------------------------------------------------------------
Dataset: nu90201039002B_bp.fits
--------------------------------------------------------------------------------

     Block Name                          Type         Dimensions
--------------------------------------------------------------------------------
Block    1: PRIMARY                        Null
Block    2: BADPIX                         Table         5 cols x 8        rows
Block    3: BADPIX2                        Table         5 cols x 11       rows
Block    4: BADPIX3                        Table         5 cols x 4        rows
Block    5: BADPIX4                        Table         5 cols x 5        rows


*_cl.evt

--------------------------------------------------------------------------------
Dataset: nu90201039002B06_cl.evt
--------------------------------------------------------------------------------

     Block Name                          Type         Dimensions
--------------------------------------------------------------------------------
Block    1: PRIMARY                        Null
Block    2: EVENTS                         Table        38 cols x 18296    rows
Block    3: GTI                            Table         1 cols x 53       rows
Block    4: BADPIX                         Table         5 cols x 45       rows
Block    5: BADPIX2                        Table         5 cols x 263      rows
Block    6: BADPIX3                        Table         5 cols x 114      rows
Block    7: BADPIX4                        Table         5 cols x 323      rows


BADPIX columns:

--------------------------------------------------------------------------------
Columns for Table Block BADPIX
--------------------------------------------------------------------------------

ColNo  Name                 Unit        Type             Range
   1   RAWX                 pixel        Byte                                X-position in raw detector coordinates
   2   RAWY                 pixel        Byte                                Y-position in raw detector coordinates
   3   TIME                 s            Real8          214279884.0:214367318.2499879897 Start Time of Bad Pixel Interval
   4   TIME_STOP            s            Real8          -Inf:+Inf            Stop Time of Bad Pixel Interval
   5   BADFLAG[2]                        Bit(2)                              Bad Pixel flag
"""
