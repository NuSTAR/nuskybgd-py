#!/usr/bin/env python3
import os
import re


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


import numpy as np
import json
import astropy.io.fits as pf
import xspec


def fpm_parse(keyword):
    if keyword not in ('FPMA', 'FPMB'):
        return False
    else:
        return keyword[-1]


def mask_from_region(regfile, refimg):
    """
    Create a pixel mask for regfile based on the image WCS info in refimg.
    """
    pass


# os.environ['NUSKYBGD_AUXIL']
auxildir = '/Users/qw/astro/nustar/nuskybgd-idl/auxil'


ratios = json.loads(open('%s/ratios.json' % auxildir).read())

# In [45]: ratiosA.keys()
# Out[45]: dict_keys(['name', 'comment', 'models'])


# bgfiles and regfiles must have the same ordering
bgfiles = ['bgd1A_sr_g30.pha', 'bgd1B_sr_g30.pha',
           'bgd2A_sr_g30.pha', 'bgd2B_sr_g30.pha',
           'bgd3A_sr_g30.pha', 'bgd3B_sr_g30.pha']

regfiles = ['bgd1A.reg', 'bgd1B.reg',
            'bgd2A.reg', 'bgd2B.reg',
            'bgd3A.reg', 'bgd3B.reg']


bgdapim = {}
bgdapim['A'] = pf.open('bgdapA.fits')[0].data
bgdapim['B'] = pf.open('bgdapB.fits')[0].data

bgddetim = {}
bgddetim['A'] = [
    pf.open('det0Aim.fits')[0].data,
    pf.open('det1Aim.fits')[0].data,
    pf.open('det2Aim.fits')[0].data,
    pf.open('det3Aim.fits')[0].data
]
bgddetim['B'] = [
    pf.open('det0Bim.fits')[0].data,
    pf.open('det1Bim.fits')[0].data,
    pf.open('det2Bim.fits')[0].data,
    pf.open('det3Bim.fits')[0].data
]


spectra = []

# Load each spectrum as a new data group
for i in range(len(bgfiles)):
    spectra.append(xspec.AllData('{num}:{num} {file}'.format(
        num=i + 1,
        file=bgfiles[i])))

# Check for valid INSTRUME keyword in spectrum header: need this to determine
# which FPM is used.

# Reference spectrum's index for the linked model parameters
# 0-based; note that xspec.AllData() is 1-based
refspec = {'A': None, 'B': None}

for i in range(xspec.AllData.nSpectra):
    spec = xspec.AllData(i + 1)

    try:
        fpm = fpm_parse(spec.fileinfo('INSTRUME'))
        if fpm is False:
            print('Spectrum %s does not have valid INSTRUME keyword.' %
                  spec.fileName)
        else:
            if fpm == 'A' and refspec['A'] is None:
                refspec['A'] = i
            elif fpm == 'B' and refspec['B'] is None:
                refspec['B'] = i

    except Exception:
        print(Exception)
        print('Spectrum %s does not have the INSTRUME keyword.'
              % spec.fileName)


# Specify the RMF and ARF files for the different model sources
for i in range(xspec.AllData.nSpectra):
    spec = xspec.AllData(i + 1)
    spec.multiresponse[1] = spec.response.rmf  # 2:apbgd
    spec.multiresponse[2] = spec.response.rmf  # 3:intbgd
    spec.multiresponse[3] = spec.response.rmf  # 4:fcxb
    spec.multiresponse[4] = '%s/diag.rmf' % auxildir  # 5:intn
    spec.multiresponse[5] = spec.response.rmf  # 6:grxe
    spec.multiresponse[1].arf = '%s/be.arf' % auxildir  # 2:apbgd
    spec.multiresponse[3].arf = '%s/fcxb%s.arf' % (
        auxildir, fpm_parse(spec.fileinfo('INSTRUME')))  # 4:fcxb
    spec.multiresponse[5].arf = '%s/be.arf' % auxildir  # 6:grxe


# Compute aperture image and detector mask based weights using each background
# region's mask

# Each background spectrum has a list of 4 values associated with each CCD:
# number of pixels in the region mask.
bgddetimsum = []

# Each background spectrum has a single value that is the sum of the aperture
# image in the region mask.
bgdapimwt = []

for i in range(xspec.AllData.nSpectra):
    spec = xspec.AllData(i + 1)
    fpm = fpm_parse(spec.fileinfo('INSTRUME'))

    regmask = mask_from_region(regfiles[i])

    detnpix = [np.sum(regmask * detim) for detim in bgddetim[fpm]]
    bgddetimsum.append(detnpix)

    bgdapimwt.append(np.sum(regmask * bgdapim[fpm]))


# Load models

# Show model 'apbgd' for data group 1
# xspec.AllModels(1, 'apbgd').show()
# ========================================================================
# Model apbgd:cutoffpl<1> Source No.: 2   Active/On
# Model Model Component  Parameter  Unit     Value
#  par  comp
#                            Data group: 1
#    1    1   cutoffpl   PhoIndex            1.00000      +/-  0.0
#    2    1   cutoffpl   HighECut   keV      15.0000      +/-  0.0
#    3    1   cutoffpl   norm                1.00000      +/-  0.0
# ________________________________________________________________________

# Get values array for parameter (3) norm of that model
# xspec.AllModels(1, 'apbgd')(3).values
# [1.0, 0.01, 0.0, 0.0, 1e+20, 1e+24]

# Toggle freezing parameters
# xspec.AllModels(1, 'apbgd')(3).frozen = True

# Link parameters (view)
# xspec.AllModels(3, 'apbgd')(3).link
# '= apbgd:p3'
# This is linked to the norm of the first data group model

# Set (omit the '=')
# xspec.AllModels(3, 'apbgd')(3).link = 'apbgd:p6'
# This is now linked to the norm of the second data group model


"""
Aperture image (apbgd)

cutoffpl

p1 (PhoIndex) and p2 (HighCut) are frozen.

p3 (norm) is defined as

0.002353
--------  x  bgdapimwt
   32

bgdapimwt = np.sum(bgdapim * regmask), for each FPM and region.

This is calculated for the reference spectra for each FPM. The other spectra
are scaled using the second part, sum(bgdapim * regmask).
"""

mod_apbgd = ratios['models'][0]['components']

m_apbgd_1 = xspec.Model('cutoffpl', 'apbgd', 2)

for i in range(xspec.AllData.nSpectra):
    spec = xspec.AllData(i + 1)
    fpm = fpm_parse(spec.fileinfo('INSTRUME'))
    m = xspec.AllModels(i + 1, 'apbgd')
    if i == refspec['A'] or i == refspec['B']:
        m.cutoffpl.PhoIndex.values = mod_apbgd['cutoffpl'][fpm]['phoindex']
        m.cutoffpl.HighECut.values = mod_apbgd['cutoffpl'][fpm]['highecut']
        m.cutoffpl.norm.values = 0.002353 / 32 * bgdapimwt[i]
    else:
        m.cutoffpl.PhoIndex.link = 'apbgd:p%d' % (3 * refspec[fpm] + 1)
        m.cutoffpl.HighECut.link = 'apbgd:p%d' % (3 * refspec[fpm] + 2)
        m.cutoffpl.norm.link = '%f * apbgd:p%d' % (
            bgdapimwt[i] / bgdapimwt[refspec[fpm]],
            3 * refspec[fpm] + 3)


"""
Instrument background (intbgd)

An apec, then many lorentz lines,

Each lorentz has 3 parameters (LineE, Width, norm); apec has 4 parameters (kT,
Abundanc, Redshift, norm).

For the reference spectra:

apec params are loaded from preset.

Lines 1-3 are loaded from preset.

Line 4 (19 keV) is loaded from preset:

norm = sum (ifactor * bgddetimsum)
       ---------------------------
            sum (bgddetimsum)


The other lines are scaled to this one using:

sum(ifactor * bgddetimsum)

Other spectra:

Lines 1-3 scale to refspec lines 1-3.

Line 4 scale to refspec line 4.

Other lines scale to line 4.

apec scale to refspec apec norm.
"""

mod_intbgd = ratios['models'][1]['components']

m_intbgd_1 = xspec.Model('apec' + '+lorentz' * (len(mod_intbgd) - 1),
                         'intbgd', 3)


# There are in total (len(mod_intbgd) * 3 + 1) parameters per spectrum. 4
# parameters for apec and 3 parameters for lorentz.
for i in range(xspec.AllData.nSpectra):
    spec = xspec.AllData(i + 1)
    fpm = fpm_parse(spec.fileinfo('INSTRUME'))
    m = xspec.AllModels(i + 1, 'intbgd')

    # Reference spectrum
    if i == refspec['A'] or i == refspec['B']:

        m.apec.kT.values = mod_intbgd['apec'][fpm]['kt']
        m.apec.Abundanc.values = mod_intbgd['apec'][fpm]['abundanc']
        m.apec.Redshift.values = mod_intbgd['apec'][fpm]['redshift']
        m.apec.norm.values = np.sum(
            bgddetimsum[i] * np.array(
                mod_intbgd['apec'][fpm]['ifactors'])
        ) / np.sum(bgddetimsum[i])

        # 3 solar lines and 19 keV line -- load them from preset
        for attr in ['lorentz', 'lorentz_2', 'lorentz_3', 'lorentz_4']:
            m_line = getattr(m, attr)
            m_line.LineE.values = mod_intbgd[attr][fpm]['linee']
            m_line.Width.values = mod_intbgd[attr][fpm]['width']
            m_line.norm.values = np.sum(
                bgddetimsum[i] * np.array(
                    mod_intbgd[attr][fpm]['ifactors'])
            ) / np.sum(bgddetimsum[i])

        # All the other lines -- lorentz_5 through lorentz_(components-1) --
        # load them from preset and link to 19 keV line (lorentz_4). There are
        # 4+3*3=13 parameters before lorentz_4.
        for attr_n in range(5, len(mod_intbgd)):
            attr = 'lorentz_%d' % attr_n
            m_line = getattr(m, attr)
            m_line.LineE.values = mod_intbgd[attr][fpm]['linee']
            m_line.Width.values = mod_intbgd[attr][fpm]['width']

            preset = np.sum(
                bgddetimsum[i] * np.array(
                    mod_intbgd[attr][fpm]['ifactors'])
            ) / np.sum(bgddetimsum[i])

            lorentz_4_norm_npar = i * m.nParameters + 16

            m_line.norm.link = '%f * intbgd:p%d' % (
                preset / m.lorentz_4.norm.values,
                refspec_norm_npar)

    else:
        m_ref = xspec.AllModels(refspec[fpm] + 1, 'intbgd')
        m_ref_npar_offset = m.nParameters * refspec[fpm]

        # Link apec to refspec
        m.apec.kT.link = 'intbgd:p%d' % (m_ref_npar_offset + 1)
        m.apec.Abundanclink = 'intbgd:p%d' % (m_ref_npar_offset + 2)
        m.apec.Redshift.link = 'intbgd:p%d' % (m_ref_npar_offset + 3)

        preset = np.sum(
            bgddetimsum[i] * np.array(
                mod_intbgd['apec'][fpm]['ifactors'])
        ) / np.sum(bgddetimsum[i])

        m.apec.norm.link = '%f * intbgd:p%d' % (
            preset / m_ref.apec.norm.values,
            m_ref_npar_offset + 4
        )

        # All lines --- load values from preset and link to refspec
        for attr_n in range(5, len(mod_intbgd)):
            if attr_n == 1:
                attr = 'lorentz'  # special case
            else:
                attr = 'lorentz_%d' % attr_n

            m_line = getattr(m, attr)
            m_line.LineE.link = 'intbgd:p%d' % (
                m_ref_npar_offset + 3 * (attr_n - 1) + 5)
            m_line.Width.link = 'intbgd:p%d' % (
                m_ref_npar_offset + 3 * (attr_n - 1) + 6)

            preset = np.sum(
                bgddetimsum[i] * np.array(
                    mod_intbgd[attr][fpm]['ifactors'])
            ) / np.sum(bgddetimsum[i])

            m_line.norm.link = '%f * intbgd:p%d' % (
                preset / getattr(m_ref, attr).norm.values,
                m_ref_npar_offset + 3 * (attr_n - 1) + 7
            )


"""
Focused CXB (fcxb)

cutoffpl, with first two parameters (PhoIndex, HighCut) tied to apbgd model
source.

The norm of cutoffpl is to be computed from

0.002353 * 1.5 * (2.45810736/3500/1000)^2 * backscal

backscal = np.sum(regmask * bgddetim) / 1000^2
"""

mod_fxcb = ratios['models'][2]['components']

m_fcxb_1 = xspec.Model('cutoffpl', 'fcxb', 4)

mod_fcxb_factor = 0.002353 * 1.5 * (2.45810736 / 3500 / 1000)**2
for i in range(xspec.AllData.nSpectra):
    spec = xspec.AllData(i + 1)
    fpm = fpm_parse(spec.fileinfo('INSTRUME'))
    m = xspec.AllModels(i + 1, 'fcxb')

    m.cutoffpl.PhoIndex.link = 'apbgd:p%d' % (3 * refspec[fpm] + 1)
    m.cutoffpl.HighECut.link = 'apbgd:p%d' % (3 * refspec[fpm] + 2)

    preset = mod_fcxb_factor / 1000**2 * np.sum(bgddetimsum[i])

    if i == refspec['A'] or i == refspec['B']:
        m.cutoffpl.norm.values = preset
    else:
        m_ref = xspec.AllModels(refspec[fpm] + 1, 'fcxb')
        m.cutoffpl.norm.link = '%f * fcxb:p%d' % (
            preset / m_ref.cutoffpl.norm.values,
            3 * refspec[fpm] + 3)


"""
intn

bknpower

First 3 parameters read from preset (PhoIndx1, BreakE, PhoIndx2).

norm = np.sum(bgddetimsum * ifactor)
"""

mod_intn = ratios['models'][3]['components']

m_intn_1 = xspec.Model('bknpower', 'intn', 5)

for i in range(xspec.AllData.nSpectra):
    spec = xspec.AllData(i + 1)
    fpm = fpm_parse(spec.fileinfo('INSTRUME'))
    m = xspec.AllModels(i + 1, 'fcxb')

    m.bknpower.PhoIndx1.values = mod_intn['bknpower'][fpm]['phoindx1']
    m.bknpower.BreakE.values = mod_intn['bknpower'][fpm]['breake']
    m.bknpower.PhoIndx2.values = mod_intn['bknpower'][fpm]['phoindx2']

    preset = np.sum(bgddetimsum * mod_intn['bknpower'][fpm]['ifactors'])

    if i == refspec['A'] or i == refspec['B']:
        m.bknpower.norm.values = preset
    else:
        m_ref = xspec.AllModels(refspec[fpm] + 1, 'intn')
        m.bknpower.norm.link = '%f * intn:p%d' % (
            preset / m_ref.bknpower.norm.values,
            4 * refspec[fpm] + 4)





