#!/usr/bin/env python3
import os
import sys
import numpy as np
import re
import json


# Input params

MSG_WARN_1 = """
WARNING: NUSKYBGD_FITAB is assuming you have supplied corresponding
   !     spectra for A and B in corresponding order, e.g.,
   !       bgdreg=["bgd1A.reg","bgd2A.reg","bgd1B.reg","bgd2B.reg"]
   !     where bgd1A.reg and bgd1B.reg, etc., cover nearly the same
   !     area on the sky and all the A spectra/regions are listed first.
"""


def nuskybgd_fitab(indir, obsid, bgdreg, specdir, specname, ab,
                   bgddir=None, header,
                   clobber=clobber, pa=pa, paramfile=paramfile,
                   nocheck=nocheck, fixfcxb=fixfcxb, grxe=None, iisrc=None,
                   fixap=fixap, nofit=nofit, tieap=tieap,
                   srcarr=None, srcdir=None, runxcm=runxcm,
                   fix_line_ratios=fix_line_ratios, no_solar=no_solar):

    absave = ab

    if srcarr is None:
        srcarr = [0] * len(bgdreg)

    if srcdir is None:
        srcdir = ''

    auxildir = os.environ['NUSKYBGD_AUXIL'] + '/'
    caldbdir = os.environ['CALDB'] + '/'

    dir_ = None
    cldir = []
    clspecdir = []
    clsrcdir = None

    for i in range(len(indir)):
        tdir = indir[i]
        if tdir[-1] != '/':
            tdir += '/'

        cldir.append(tdir + obsid[i] + '/event_cl/')
        clspecdir.append(tdir + obsid[i] + '/event_cl/' + specdir + '/')
        clsrcdir.append(tdir + obsid[i] + '/event_cl/' + srcdir + '/')

        if bgddir is None:
            dir_.append(cldir[i])
        else:
            dir_.append(cldir[i] + bgddir + '/')
            if not os.path.exists(dir[i]):
                os.makedirs(dir[i])  # Makes parent dirs if needed

    if grxe is None:
        grxe = 0

    if iisrc is None:
        iisrc = [0] * len(bgdreg)

    if not isinstance(ab, str):
        raise TypeError('Input ab must be a string.')

    iidir = None

    if ab.lower() == 'ab':
        if len(bgdreg) % (2 * len(dir)) != 0:
            print('Error: bgdreg is not multiple of 2 * dir, stopping')
            sys.exit(1)

        abarr = (['A'] * len(bgdreg) / 2 / len(dir) +
                 ['B'] * len(bgdreg) / 2 / len(dir)) * len(dir)

        for i in range(len(dir)):
            iidir.append([i] * (len(bgdreg) / len(dir)))

        print(MSG_WARN_1)

        if len(bgdreg) % 2 != 0 or len(specname) % 2 != 0:
            print('Error: Odd number of inputs in bgdreg or specname')
            sys.exit(1)

    elif ab.lower() in ('a', 'b'):
        abarr = [ab] * len(bgdreg)

        for i in range(len(dir)):
            iidir.append([i] * (len(bgdreg) / len(dir)))

    else:
        print('Error: FPM input must be A, B, or AB')
        sys.exit(1)

    xcmfile = clspecdir[0] + specdir + '.xcm'

    output = []









    tb_A = np.loadtxt(auxildir + 'ratiosA.dat')


    # 7 columns
    aeline, awidth, af0, af1, af2, af3 = (
        tb_A[:, 0], tb_A[:, 1], tb_A[:, 2], tb_A[:, 3], tb_A[:, 4], tb_A[:, 5])

    nlines = len(aeline) - 2

    eline = np.zeros((nlines, 2), dtype=np.float32)
    width = np.zeros((nlines, 2), dtype=np.float32)
    ifactors = np.zeros((len(aeline), 4, 2), dtype=np.float32)
    index1 = np.zeros((2, 2), dtype=np.float32)
    index2 = np.zeros((2, 2), dtype=np.float32)
    ebreak = np.zeros((2, 2), dtype=np.float32)


    for iab in ['A', 'B']:
        tb = np.loadtxt('%sratios%s.dat' % (auxildir, iab))




def read_model_table(fpm):
    if fpm not in ('A', 'B'):
        print('Error: FPM must be A or B')
        return False
    model = []
    fh = open('%sratios%s.dat' % (auxildir, fpm))
    lines = fh.read().splitlines()
    for line in lines:
        pars = re.sub(r'\s+', ',', line.strip()).split(',')
        model.append(pars)
    return model



def read_model_template(fpm):
    if fpm not in ('A', 'B'):
        print('Error: FPM must be A or B')
        return False

    auxildir = os.environ['NUSKYBGD_AUXIL']
    fh = open('%s/ratios%s.json' % (auxildir, fpm))

    models = json.loads(fh.read())

    print('Loaded XSPEC models template: %s' % models['name'])
    print('It has the following components:')

    for m in models['models']:
        print(m['name'])
        print(m['formula'])



