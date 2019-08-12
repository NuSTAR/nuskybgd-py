#!/usr/bin/env python3
import os
import json
import re
import numpy as np
import xspec
# !!Load xspec before astropy!! When testing using latest xspec (6.26.1) and
# astropy (3.2.1), xspec has cfitsio with SONAME 8, and loading astropy before
# it throws an error when loading spectrum file, stating that the loaded
# cfitsio library has SONAME 7. Is it astropy that's been compiled with an
# older version of cfitsio, causing a problem if the older library loads
# first? Because the only library in the rpath is heasoft's latest cfitsio
# library.
import astropy.io.fits as pf
import pyregion


# Input params
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

refimgf = 'bgdapA.fits'

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


MSG_WARN_1 = """
WARNING: NUSKYBGD_FITAB is assuming you have supplied corresponding
   !     spectra for A and B in corresponding order, e.g.,
   !       bgdreg=["bgd1A.reg","bgd2A.reg","bgd1B.reg","bgd2B.reg"]
   !     where bgd1A.reg and bgd1B.reg, etc., cover nearly the same
   !     area on the sky and all the A spectra/regions are listed first.
"""


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


def fpm_parse(keyword):
    if keyword not in ('FPMA', 'FPMB'):
        return False
    else:
        return keyword[-1]


def mask_from_region(regfile, refimg):
    """
    Create a pixel mask for regfile based on the image WCS info in refimg.

    Uses the pyregion module. Please keep to using circle, box, and ellipse
    shapes in DS9, fk5 format to avoid unexpected behavior.
    """
    fh = pf.open(refimg)
    reg = pyregion.open(regfile)
    mask = reg.get_mask(hdu=fh[0])
    return mask


def nuskybgd_fitab(indir, obsid, bgdreg, specdir, specname, ab,
                   bgddir=None, header,
                   clobber=clobber, pa=pa, paramfile=paramfile,
                   nocheck=nocheck, fixfcxb=fixfcxb, grxe=None, iisrc=None,
                   fixap=fixap, nofit=nofit, tieap=tieap,
                   srcarr=None, srcdir=None, runxcm=runxcm,
                   fix_line_ratios=fix_line_ratios, no_solar=no_solar):
    pass


spectra = []

xspec.DataManager.clear(0)  # Clear any existing loaded data

# Load each spectrum as a new data group
for i in range(len(bgfiles)):
    spectra.append(xspec.AllData('{num}:{num} {file}'.format(
        prefix=prefix,
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

    regmask = mask_from_region(regfiles[i], refimgf)

    detnpix = [np.sum(regmask * detim) for detim in bgddetim[fpm]]
    bgddetimsum.append(detnpix)

    bgdapimwt.append(np.sum(regmask * bgdapim[fpm]))


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
        m.cutoffpl.PhoIndex.frozen = True
        m.cutoffpl.HighECut.frozen = True
        m.cutoffpl.norm.values = 0.002353 / 32 * bgdapimwt[i]
    else:
        m.cutoffpl.PhoIndex.link = 'apbgd:p%d' % (3 * refspec[fpm] + 1)
        m.cutoffpl.HighECut.link = 'apbgd:p%d' % (3 * refspec[fpm] + 2)
        m.cutoffpl.norm.link = '%e * apbgd:p%d' % (
            bgdapimwt[i] / bgdapimwt[refspec[fpm]],
            3 * refspec[fpm] + 3)


"""
Instrument background (intbgd)

An apec, then many lorentz lines. Each lorentz has 3 parameters (LineE, Width,
norm); apec has 4 parameters (kT, Abundanc, Redshift, norm).

For the reference spectra:

apec params are loaded from preset.

Lines 1-3 (lorentz, lorentz_3, lorentz_4) are loaded from preset.

Line 4 (19 keV, lorentz_5) is loaded from preset:

norm = sum (ifactor * bgddetimsum)
       ---------------------------
            sum (bgddetimsum)


The other lines are scaled to this one using:

sum(ifactor * bgddetimsum)

Other spectra:

Lines 1-3 scale to refspec lines 1-3 (lorentz, lorentz_3, lorentz_4).

Line 4 scale to refspec line 4 (lorentz_5).

Other lines scale to line 4 (lorentz_5).

apec scale to refspec apec norm.

Note: due to perculiar model component numbering by XSPEC,
as in the following,

Model intbgd:apec<1> + lorentz<2> + lorentz<3> + lorentz<4> + ...

In this case, the lorentz components have references lorentz, lorentz_3,
lorentz_4, ... etc. (skipping over lorentz_2). This numbering appears to be
for all components, not just repeated ones, except the first occurrence omits
the label.
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
        m.apec.kT.frozen = True
        m.apec.Abundanc.frozen = True
        m.apec.Redshift.frozen = True
        m.apec.norm.values = np.sum(
            bgddetimsum[i] * np.array(
                mod_intbgd['apec'][fpm]['ifactors'])
        ) / np.sum(bgddetimsum[i])

        # 3 solar lines and 19 keV line -- load them from preset
        for attr in ['lorentz', 'lorentz_3', 'lorentz_4', 'lorentz_5']:
            m_line = getattr(m, attr)
            m_line.LineE.values = mod_intbgd[attr][fpm]['linee']
            m_line.Width.values = mod_intbgd[attr][fpm]['width']
            m_line.LineE.frozen = True
            m_line.Width.frozen = True
            m_line.norm.values = np.sum(
                bgddetimsum[i] * np.array(
                    mod_intbgd[attr][fpm]['ifactors'])
            ) / np.sum(bgddetimsum[i])

        # All the other lines -- lorentz_6 through lorentz_(components) --
        # load them from preset and link to 19 keV line (lorentz_5). There are
        # 4+3*3=13 parameters before lorentz_5.
        for attr_n in range(5, len(mod_intbgd)):
            attr = 'lorentz_%d' % (attr_n + 1)
            m_line = getattr(m, attr)
            m_line.LineE.values = mod_intbgd[attr][fpm]['linee']
            m_line.Width.values = mod_intbgd[attr][fpm]['width']
            m_line.LineE.frozen = True
            m_line.Width.frozen = True

            preset = np.sum(
                bgddetimsum[i] * np.array(
                    mod_intbgd[attr][fpm]['ifactors'])
            ) / np.sum(bgddetimsum[i])

            lorentz_5_norm_npar = i * m.nParameters + 16

            m_line.norm.link = '%e * intbgd:p%d' % (
                preset / m.lorentz_5.norm.values[0],
                lorentz_5_norm_npar)

    else:
        m_ref = xspec.AllModels(refspec[fpm] + 1, 'intbgd')
        m_ref_npar_offset = m.nParameters * refspec[fpm]

        # Link apec to refspec
        m.apec.kT.link = 'intbgd:p%d' % (m_ref_npar_offset + 1)
        m.apec.Abundanc.link = 'intbgd:p%d' % (m_ref_npar_offset + 2)
        m.apec.Redshift.link = 'intbgd:p%d' % (m_ref_npar_offset + 3)

        preset = np.sum(
            bgddetimsum[i] * np.array(
                mod_intbgd['apec'][fpm]['ifactors'])
        ) / np.sum(bgddetimsum[i])

        m.apec.norm.link = '%e * intbgd:p%d' % (
            preset / m_ref.apec.norm.values[0],
            m_ref_npar_offset + 4
        )

        # All lines --- load values from preset and link to refspec
        for attr_n in range(1, len(mod_intbgd)):
            if attr_n == 1:
                attr = 'lorentz'  # special case
            else:
                attr = 'lorentz_%d' % (attr_n + 1)

            m_line = getattr(m, attr)
            m_line.LineE.link = 'intbgd:p%d' % (
                m_ref_npar_offset + 3 * (attr_n - 1) + 5)
            m_line.Width.link = 'intbgd:p%d' % (
                m_ref_npar_offset + 3 * (attr_n - 1) + 6)

            preset = np.sum(
                bgddetimsum[i] * np.array(
                    mod_intbgd[attr][fpm]['ifactors'])
            ) / np.sum(bgddetimsum[i])

            m_line.norm.link = '%e * intbgd:p%d' % (
                preset / getattr(m_ref, attr).norm.values[0],
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
        m.cutoffpl.norm.link = '%e * fcxb:p%d' % (
            preset / m_ref.cutoffpl.norm.values[0],
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
    m = xspec.AllModels(i + 1, 'intn')

    m.bknpower.PhoIndx1.values = mod_intn['bknpower'][fpm]['phoindx1']
    m.bknpower.BreakE.values = mod_intn['bknpower'][fpm]['breake']
    m.bknpower.PhoIndx2.values = mod_intn['bknpower'][fpm]['phoindx2']
    m.bknpower.PhoIndx1.frozen = True
    m.bknpower.BreakE.frozen = True
    m.bknpower.PhoIndx2.frozen = True

    preset = np.sum(
        bgddetimsum[i] * np.float64(mod_intn['bknpower'][fpm]['ifactors']))

    if i == refspec['A'] or i == refspec['B']:
        m.bknpower.norm.values = preset
    else:
        m_ref = xspec.AllModels(refspec[fpm] + 1, 'intn')
        m.bknpower.norm.link = '%e * intn:p%d' % (
            preset / m_ref.bknpower.norm.values[0],
            4 * refspec[fpm] + 4)


xspec.AllData.ignore('**-3. 150.-**')
xspec.Fit.method = 'leven 30000 1e-4'
xspec.Fit.statMethod = 'chi'
xspec.Fit.perform()


i = 1
while i < 101:  # Retry 100 times...
    if i > 1:
        savefile = 'bgdparams%d.xcm' % i
    else:
        savefile = 'bgdparams.xcm'
    if not os.path.exists(savefile):
        xspec.Xset.save(savefile, info='a')
        print('Saved results to %s.' % savefile)
        break
    else:
        i += 1
