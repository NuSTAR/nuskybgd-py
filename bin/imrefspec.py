#!/usr/bin/env python3
"""
Creates reference spectra of the aperture background and fcxb background
components, used by nuskybgd_image, to produce counts images for arbitrary
energy bands.  Only needs to be run once unless the nuabs parameters change in
the auxil directory.

Usage:

imrefspec.py AB 0123

Creates fake spectra for all CCDs of both modules.

imrefspec.py A 2

Creates fake spectrum for det 2 of module A.

Required files:

RMF files that include detector absorption (see absrmf.py).

Output:

Simulated spectra using Xspec's fakeit.
"""

import xspec
import os
import sys

_NUSKYBGD_AUX_ENV = 'NUSKYBGD_AUXIL'


def imrefspec(fpmlist, idetlist, auxil,
              cxbp={'pl': 1.29, 'ecut': 41.13, 'norm': 0.002353},
              texp=1.0e9):

    queue = []

    # Check for RMF files
    for idet in idetlist:
        idet = int(idet)
        if idet not in (0, 1, 2, 3):
            continue
        for ab in fpmlist:
            if ab not in ('A', 'B'):
                continue
            rmffile = 'det%d%s.rmf' % (idet, ab)
            if os.path.exists(rmffile):
                queue.append((idet, ab))
                print('Queued FPM%s det%d: %s' %
                      (ab, idet, rmffile))
            else:
                print('Skipped FPM%s det%d: %s not found.' %
                      (ab, idet, rmffile))

    if len(queue) == 0:
        print('Nothing to be done.')
        return True

    m = xspec.Model('cutoffpl')
    m.cutoffpl.PhoIndex = cxbp['pl']
    m.cutoffpl.HighECut = cxbp['ecut']
    m.cutoffpl.norm = cxbp['norm']

    for idet, ab in queue:
        # If there is existing spectrum loaded, RMF and ARF in FakeSettings
        # will be ignored (the 'from existing spectra' flow). We need to do
        # the 'from scratch' flow.
        xspec.AllData.clear()

        fakesettings_ap = xspec.FakeitSettings(
            response='det%d%s.rmf' % (idet, ab),
            arf='%s/be.arf' % auxil,
            exposure=texp,
            correction=0.0,
            fileName='aperbgd%s_det%d.pha' % (ab, idet)
        )
        fakesettings_fcxb = xspec.FakeitSettings(
            response='det%d%s.rmf' % (idet, ab),
            arf='%s/fcxb%s.arf' % (auxil, ab),
            exposure=texp,
            correction=0.0,
            fileName='fcxbbgd%s_det%d.pha' % (ab, idet)
        )
        xspec.AllData.fakeit(
            nSpectra=2,
            settings=[fakesettings_ap, fakesettings_fcxb],
            applyStats=False,
            filePrefix='')

    return True


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(0)

    # Check auxil dir setting
    if _NUSKYBGD_AUX_ENV not in os.environ:
        print('Please set the NUSKYBGD_AUXIL environment variable first.')
        sys.exit(1)

    # Check auxil dir is OK...
    if not os.path.exists('%s/be.arf' % os.environ[_NUSKYBGD_AUX_ENV]):
        print('Error: be.arf not in %s' % os.environ[_NUSKYBGD_AUX_ENV])

    imrefspec(
        list(sys.argv[1]),
        list(sys.argv[2]),
        os.environ[_NUSKYBGD_AUX_ENV])

    sys.exit(0)
