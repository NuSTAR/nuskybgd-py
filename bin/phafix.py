#!/usr/bin/env python3

# (c) 2019 Qian Wang
# This file is part of nuskybgd released under MIT License.
# See LICENSE file for full details.

"""
phafix.py

Fix relative paths of response and auxiliary file keywords in spectrum files.

In all use cases for nuskybgd, any response (RMF) and auxiliary (ARF) files
are co-located with the PI spectrum file. When loading the spectrum file in
the current working directory, the file paths referencing those other files
should not point somewhere else. This quick fix tool strips the 'dirname' part
of file paths leaving only the 'basename' part, in the BACKFILE, CORRFILE,
RESPFILE, and ANCRFILE keywords in the header of any extensions that have
HDUCLAS1='SPECTRUM' or EXTNAME='SPECTRUM'.
"""

import sys
import os
import astropy.io.fits as pf


def phafix(specfile):
    t = pf.open(specfile)

    for ext in t:
        if (
            ('HDUCLAS1' in ext.header and
                ext.header['HDUCLAS1'] == 'SPECTRUM') or
            ('EXTNAME' in ext.header and
                ext.header['EXTNAME'] == 'SPECTRUM')):
            for kw in ('BACKFILE', 'CORRFILE', 'RESPFILE', 'ANCRFILE'):
                if kw in ext.header:
                    old = ext.header[kw]
                    ext.header[kw] = os.path.basename(ext.header[kw])
                    # Provide some feedback
                    if ext.header[kw] != old:
                        print('%s -> %s' % (old, ext.header[kw]))
    t.writeto(specfile, overwrite=True)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(0)

    specfile = sys.argv[1]
    if not os.path.exists(specfile):
        print('File %s does not exist.' % specfile)
        sys.exit(1)
    else:
        phafix(specfile)

    sys.exit(0)
