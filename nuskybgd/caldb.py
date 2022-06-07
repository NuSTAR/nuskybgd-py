# (c) 2019 Qian Wang
# This file is part of nuskybgd released under MIT License.
# See LICENSE file for full details.

import astropy.io.fits as pf
import os
import numpy as np
import datetime
import sys


class CalDB:
    """
    CalDB class to look up appropriate CALDB file. The logic for selecting the
    file is in ``getcaldbfile()``. Various functions such as ``getDETABS()``
    exist, calling ``getcaldbfile()`` with certain preset arguments.

    Examples
    --------

    >>> caldb_path = '/usr/local/heasarc/CALDB'
    >>> telescope = 'NUSTAR'
    >>> instrument = 'FPMA'  # FPMA or FPMB
    >>> cnam = 'PIXPOS'
    >>> obsutctime = '2014-01-01T00:00:00'
    >>> detnam = 'DET0'  # DET0 - DET3
    >>> cal = CalDB(caldb_path)
    >>> caldbfile = cal.getcaldbfile(caldb_path, telescope, instrument,
    ...                              cnam, detnam, obsutctime)
    CALDB files for PIXPOS FPMA DET0:
    /usr/local/heasarc/CALDB/
    data/nustar/fpm/bcf/pixpos/nuApixpos20100101v001.fits
    data/nustar/fpm/bcf/pixpos/nuApixpos20100101v002.fits
    data/nustar/fpm/bcf/pixpos/nuApixpos20100101v003.fits
    data/nustar/fpm/bcf/pixpos/nuApixpos20100101v004.fits
    data/nustar/fpm/bcf/pixpos/nuApixpos20100101v005.fits
    data/nustar/fpm/bcf/pixpos/nuApixpos20100101v006.fits
    data/nustar/fpm/bcf/pixpos/nuApixpos20100101v007.fits ***SELECTED***
    >>> caldbfile
    'data/nustar/fpm/bcf/pixpos/nuApixpos20100101v007.fits'

    >>> caldbfile = cal.getcaldbfile(caldb_path, telescope, instrument,
    ...                              'BADPIX', detnam, obsutctime)
    CALDB files for BADPIX FPMA DET0:
    /usr/local/heasarc/CALDB/
    data/nustar/fpm/bcf/badpix/nuAbadpix20100101v001.fits
    data/nustar/fpm/bcf/badpix/nuAbadpix20100101v002.fits ***SELECTED***
    >>> caldbfile
    'data/nustar/fpm/bcf/badpix/nuAbadpix20100101v002.fits'
    """

    def __init__(self, path):
        self._CalDBPath = path
        self._index = self._loadCalDBIndex()
        self.telescope = 'NUSTAR'

    def _loadCalDBIndex(self):
        indexfile = '%s/data/nustar/fpm/caldb.indx' % self._CalDBPath
        if not os.path.exists(indexfile):
            print('Error: CALDB index file not found.')
            return None
        return pf.open(indexfile)

    # ------------------------------------

    def getGRPPSF(self, instrument, obsutctime):
        """Call ``getcaldbfile()`` with ``cnam='GRPPSF'``."""
        return self.getcaldbfile(
            self.telescope, instrument, 'GRPPSF', None, obsutctime)

    def getDETABS(self, instrument, detnam, obsutctime):
        """Call ``getcaldbfile()`` with ``cnam='DETABS'``."""
        return self.getcaldbfile(
            self.telescope, instrument, 'DETABS', detnam, obsutctime)

    def getBEABS(self, instrument, obsutctime):
        """Call ``getcaldbfile()`` with ``cnam='BEABSPAR'``."""
        return self.getcaldbfile(
            self.telescope, instrument, 'BEABSPAR', None, obsutctime)
    
    def getGRPRMF(self, instrument, detnam, obsutctime):
        """Call ``getcaldbfile()`` with ``cnam='GRPRMF'``."""
        return self.getcaldbfile(
            self.telescope, instrument, 'GRPRMF', detnam, obsutctime)

    def getSPECRESP(self, instrument, obsutctime):
        """Call ``getcaldbfile()`` with ``cnam='SPECRESP'``."""
        return self.getcaldbfile(
            self.telescope, instrument, 'SPECRESP', None, obsutctime)

    def getTVIGNET(self, instrument, obsutctime):
        """Call ``getcaldbfile()`` with ``cnam='TVIGNET'``."""
        return self.getcaldbfile(
            self.telescope, instrument, 'TVIGNET', None, obsutctime)

    def getINSTRMAP(self, instrument, obsutctime):
        """Call ``getcaldbfile()`` with ``cnam='INSTRMAP'``."""
        return self.getcaldbfile(
            self.telescope, instrument, 'INSTRMAP', None, obsutctime)

    def getBADPIX(self, instrument, detnam, obsutctime):
        """Call ``getcaldbfile()`` with ``cnam='BADPIX'``."""
        return self.getcaldbfile(
            self.telescope, instrument, 'BADPIX', detnam, obsutctime)

    def getPIXPOS(self, instrument, detnam, obsutctime):
        """Call ``getcaldbfile()`` with ``cnam='PIXPOS'``."""
        return self.getcaldbfile(
            self.telescope, instrument, 'PIXPOS', detnam, obsutctime)

    def getAPERTURE(self, instrument, obsutctime):
        """Call ``getcaldbfile()`` with ``cnam='APERTURE'``."""
        return self.getcaldbfile(
            self.telescope, instrument, 'APERTURE', None, obsutctime)

    def getRMF(self, instrument, detnam, obsutctime):
        """
        First look for the GRPRMF file using ``getGRPRMF()``, then look for
        the RMF file for the specified DETNAM.
        """
        idet = int(detnam[-1])
        grprmffile = self.getGRPRMF(instrument, detnam, obsutctime)
        grprmfpath = '%s/%s' % (self._CalDBPath, grprmffile)
        grprmffh = pf.open(grprmfpath)
        if len(grprmffh[idet + 1].data) > 1:
            print('Warning: more than one RMF in GRPRMF %s for %s. '
                  'Using first entry.' % (grprmffile, detnam))
        rmfdir = os.path.dirname(grprmfpath)
        return ('%s/%s' % (
            rmfdir,
            grprmffh[idet + 1].data[0].field('RMFFILE')
        ))

    # ------------------------------------

    def getcaldbfile(
            self, telescope, instrument, cnam,
            detnam, obsutctime, verbose=False, diag=False):
        """
        The CALDB file is selected in the following manner:

        - Match TELESCOP, INSTRUME, CAL_CNAM keywords in header
        - Select files whose CAL_VSD and CAL_VST pre-date the requested time
        - Filter by DETNAM, if specified
        - (GRPRMF only) select files with CAL_CBD containing DEPTHCUT(NOMINAL)
        - Sort file names and pick the last one

        The CALDB files are typically named with dates and versions in such a
        way that the last entry after sorting is the latest file.
        """
        table = self._index[1].data

        halt = False
        for kw in ('TELESCOP', 'INSTRUME', 'CAL_CNAM',
                   'CAL_DIR', 'CAL_FILE', 'DETNAM', 'CAL_VSD', 'CAL_VST'):
            if kw not in table.names:
                print('Error: column %s not found in CALDB table.')
                halt = True
        if halt:
            print('The following columns exist in the specified file:')
            print(' '.join(table.names))
            return False

        inx = np.where(
            (table.field('CAL_CNAM') == cnam) *
            (table.field('TELESCOP') == telescope) *
            (table.field('INSTRUME') == instrument)
        )[0]

        if len(inx) == 0:
            print('No suitable file for %s in CALDB for %s.' % (
                cnam, instrument))
            return False

        obsdt = datetime.datetime.strptime(obsutctime, '%Y-%m-%dT%H:%M:%S')
        startdt = np.array([
            datetime.datetime.strptime('%sT%s' % (
                vs[0], vs[1]), '%Y-%m-%dT%H:%M:%S')
            for vs in zip(
                table[inx].field('CAL_VSD'),
                table[inx].field('CAL_VST'))
        ])

        # Filter start date later than obsdate
        inx2 = np.where(startdt < obsdt)[0]
        startdt = startdt[inx2]
        inx = inx[inx2]
        if len(inx) == 0:
            print('No suitable file for this obsdate in CALDB.')
            return False

        # GRPRMF filter --- CAL_CBD starts with DEPTHCUT(NOMINAL)
        if cnam == 'GRPRMF':
            inx22 = []
            for i in range(len(inx)):
                if (table[inx[i]].field('CAL_CBD')[:17] ==
                        'DEPTHCUT(NOMINAL)'):
                    inx22.append(i)
            inx = inx[np.array(inx22)]

        # Filter by DETNAM
        if detnam is not None:
            inx3 = np.where(table[inx].field('DETNAM') == detnam)[0]
            inx = inx[inx3]

        filelist = sorted(['%s/%s' % (_[0], _[1]) for _ in
                           zip(table[inx].field('CAL_DIR'),
                               table[inx].field('CAL_FILE'))])

        if len(filelist) == 0:
            print('No file matched DETNAM filter; try using detnam=None.')
            return False

        if verbose:
            print('CALDB files for %s %s %s %s:' % (
                cnam, instrument, detnam, obsutctime))
            print('%s/' % self._CalDBPath)
            print('\n'.join(filelist[:-1]))
            print('%s ***SELECTED***' % filelist[-1])

        if not diag:
            return filelist[-1]
        else:
            return table[inx]
