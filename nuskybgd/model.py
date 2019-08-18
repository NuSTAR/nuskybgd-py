"""
Functions for nuskybgd modelling, many involving interacting with xspec.
"""
import xspec
import os
import json
import numpy as np
from . import util


def check_bgdinfofile(bgdinfofile):
    """
    Check that the background info file has the required items.
    """
    if not os.path.exists(bgdinfofile):
        print('Error: file %s not found.' % bgdinfofile)
        return False

    bgdinfo = json.loads(open(bgdinfofile).read())

    problem = False
    for key in (
        'bgfiles', 'regfiles', 'refimgf', 'bgdapfiles', 'bgddetfiles'
    ):
        if key not in bgdinfo:
            problem = True
            print('%s not found in background info file.' % key)

    # Same number of items in bgfiles and regfiles:
    if len(bgdinfo['bgfiles']) != len(bgdinfo['regfiles']):
        problem = True
        print('bgfiles and regfiles must have the same number of entries.')

    # A and B keys in bgdapfiles and bgddetfiles
    if ('A' not in bgdinfo['bgdapfiles'] or
            'B' not in bgdinfo['bgdapfiles'] or
            'A' not in bgdinfo['bgddetfiles'] or
            'B' not in bgdinfo['bgddetfiles']):
        problem = True
        print('bgdapfiles and bgddetfiles must have A and B keys.')

    # Check files exist
    queue = []
    for _ in bgdinfo['bgfiles']:
        if isinstance(_, str):
            queue.append(_)

    for _ in bgdinfo['regfiles']:
        if isinstance(_, str):
            queue.append(_)

    for _ in bgdinfo['bgddetfiles']['A']:
        if isinstance(_, str):
            queue.append(_)

    for _ in bgdinfo['bgddetfiles']['B']:
        if isinstance(_, str):
            queue.append(_)

    queue.append(bgdinfo['refimgf'])
    queue.append(bgdinfo['bgdapfiles']['A'])
    queue.append(bgdinfo['bgdapfiles']['B'])

    for _ in queue:
        if not os.path.exists(_):
            problem = True
            print('Error: file %s not found.' % _)

    if problem:
        return False
    else:
        return bgdinfo


def addmodel_apbgd(presets, refspec, bgdapimwt, model_num, model_name='apbgd'):
    """
    Xspec model component 2: apbgd (aperture image background)

    cutoffpl

    p1 (PhoIndex) and p2 (HighCut) are frozen.

    p3 (norm) is defined as

    0.002353
    --------  x  bgdapimwt
       32

    bgdapimwt = np.sum(bgdapim * regmask), for each FPM and region.

    This is calculated for the reference spectra for each FPM. The other
    spectra are scaled using the second part, sum(bgdapim * regmask).

    Inputs:

    presets - Preset model parameter values read from file.

    refspec - Dictionary with 'A' and 'B' keys whose values are the 0-based
    index of the reference spectrum for FPMA and FPMB.

    bgdapimwt - List of sum of aperture image inside background region.

    model_num - Model component number in Xspec, cannot be the same as another
    model or it will replace.

    model_name - Model component name in Xspec, default 'apbgd', cannot be the
    same as another model.
    """
    mod_apbgd = presets['models'][0]['components']

    xspec.Model('cutoffpl', model_name, model_num)

    for i in range(xspec.AllData.nSpectra):
        spec = xspec.AllData(i + 1)
        fpm = util.fpm_parse(spec.fileinfo('INSTRUME'))
        m = xspec.AllModels(i + 1, model_name)
        if i == refspec['A'] or i == refspec['B']:
            m.cutoffpl.PhoIndex.values = mod_apbgd['cutoffpl'][fpm]['phoindex']
            m.cutoffpl.HighECut.values = mod_apbgd['cutoffpl'][fpm]['highecut']
            m.cutoffpl.PhoIndex.frozen = True
            m.cutoffpl.HighECut.frozen = True
            m.cutoffpl.norm.values = 0.002353 / 32 * bgdapimwt[i]
        else:
            m.cutoffpl.PhoIndex.link = '%s:p%d' % (
                model_name, 3 * refspec[fpm] + 1)
            m.cutoffpl.HighECut.link = '%s:p%d' % (
                model_name, 3 * refspec[fpm] + 2)
            m.cutoffpl.norm.link = '%e * %s:p%d' % (
                bgdapimwt[i] / bgdapimwt[refspec[fpm]],
                model_name,
                3 * refspec[fpm] + 3)


def addmodel_intbgd(presets, refspec, bgddetimsum, model_num,
                    model_name='intbgd'):
    """
    Xspec model component 3: intbgd (instrument background)

    An apec, then many lorentz lines. Each lorentz has 3 parameters (LineE,
    Width, norm); apec has 4 parameters (kT, Abundanc, Redshift, norm).

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
    lorentz_4, ... etc. (skipping over lorentz_2). This numbering appears to
    be for all components, not just repeated ones, except the first occurrence
    omits the label.

    Inputs:

    presets - Preset model parameter values read from file.

    refspec - Dictionary with 'A' and 'B' keys whose values are the 0-based
    index of the reference spectrum for FPMA and FPMB.

    bgddetimsum - List of area inside background region for each CCD.

    model_num - Model component number in Xspec, cannot be the same as another
    model or it will replace.

    model_name - Model component name in Xspec, default 'intbgd', cannot be
    the same as another model.
    """
    mod_intbgd = presets['models'][1]['components']

    xspec.Model('apec' + '+lorentz' * (len(mod_intbgd) - 1),
                model_name, model_num)

    # There are in total (len(mod_intbgd) * 3 + 1) parameters per spectrum. 4
    # parameters for apec and 3 parameters for lorentz.
    for i in range(xspec.AllData.nSpectra):
        spec = xspec.AllData(i + 1)
        fpm = util.fpm_parse(spec.fileinfo('INSTRUME'))
        m = xspec.AllModels(i + 1, model_name)

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
            # load them from preset and link to 19 keV line (lorentz_5). There
            # are 4+3*3=13 parameters before lorentz_5.
            for attr_n in range(5, len(mod_intbgd)):
                attr = 'lorentz_%d' % (attr_n + 1)
                m_line = getattr(m, attr)
                m_line.LineE.values = mod_intbgd[attr][fpm]['linee']
                m_line.Width.values = mod_intbgd[attr][fpm]['width']
                m_line.LineE.frozen = True
                m_line.Width.frozen = True

                norm_preset = np.sum(
                    bgddetimsum[i] * np.array(
                        mod_intbgd[attr][fpm]['ifactors'])
                ) / np.sum(bgddetimsum[i])

                lorentz_5_norm_npar = i * m.nParameters + 16

                m_line.norm.link = '%e * %s:p%d' % (
                    norm_preset / m.lorentz_5.norm.values[0],
                    model_name,
                    lorentz_5_norm_npar)

        else:
            # Link to ref spectrum
            m_ref = xspec.AllModels(refspec[fpm] + 1, model_name)
            m_ref_npar_offset = m.nParameters * refspec[fpm]

            # Link apec to refspec
            m.apec.kT.link = '%s:p%d' % (
                model_name, m_ref_npar_offset + 1)
            m.apec.Abundanc.link = '%s:p%d' % (
                model_name, m_ref_npar_offset + 2)
            m.apec.Redshift.link = '%s:p%d' % (
                model_name, m_ref_npar_offset + 3)

            norm_preset = np.sum(
                bgddetimsum[i] * np.array(
                    mod_intbgd['apec'][fpm]['ifactors'])
            ) / np.sum(bgddetimsum[i])

            m.apec.norm.link = '%e * %s:p%d' % (
                norm_preset / m_ref.apec.norm.values[0],
                model_name,
                m_ref_npar_offset + 4
            )

            # All lines --- load values from preset and link to refspec
            for attr_n in range(1, len(mod_intbgd)):
                if attr_n == 1:
                    attr = 'lorentz'  # special case
                else:
                    attr = 'lorentz_%d' % (attr_n + 1)

                m_line = getattr(m, attr)
                m_line.LineE.link = '%s:p%d' % (
                    model_name,
                    m_ref_npar_offset + 3 * (attr_n - 1) + 5)
                m_line.Width.link = '%s:p%d' % (
                    model_name,
                    m_ref_npar_offset + 3 * (attr_n - 1) + 6)

                norm_preset = np.sum(
                    bgddetimsum[i] * np.array(
                        mod_intbgd[attr][fpm]['ifactors'])
                ) / np.sum(bgddetimsum[i])

                m_line.norm.link = '%e * %s:p%d' % (
                    norm_preset / getattr(m_ref, attr).norm.values[0],
                    model_name,
                    m_ref_npar_offset + 3 * (attr_n - 1) + 7
                )


def addmodel_fcxb(refspec, bgddetimsum, model_num,
                  model_name='fcxb', apbgd_name='apbgd'):
    """
    Xspec model component 4: fcxb (focused CXB)

    Dependency: the aperture image background model must have already been
    added.

    cutoffpl, with first two parameters (PhoIndex, HighCut) tied to apbgd model
    source.

    The norm of cutoffpl is to be computed from

    0.002353 * 1.5 * (2.45810736/3500/1000)^2 * backscal

    backscal = np.sum(regmask * bgddetim) / 1000^2

    Inputs:

    refspec - Dictionary with 'A' and 'B' keys whose values are the 0-based
    index of the reference spectrum for FPMA and FPMB.

    bgddetimsum - List of area inside background region for each CCD.

    model_num - Model component number in Xspec, cannot be the same as another
    model or it will replace.

    model_name - Model component name in Xspec, default 'fcxb', cannot be the
    same as another model.

    apbgd_name - Model component name of the aperture background model, must
    have been already added.
    """
    xspec.Model('cutoffpl', model_name, model_num)

    mod_fcxb_factor = 0.002353 * 1.5 * (2.45810736 / 3500 / 1000)**2

    for i in range(xspec.AllData.nSpectra):
        spec = xspec.AllData(i + 1)
        fpm = util.fpm_parse(spec.fileinfo('INSTRUME'))
        m = xspec.AllModels(i + 1, model_name)

        m.cutoffpl.PhoIndex.link = '%s:p%d' % (
            apbgd_name, 3 * refspec[fpm] + 1)
        m.cutoffpl.HighECut.link = '%s:p%d' % (
            apbgd_name, 3 * refspec[fpm] + 2)

        norm_preset = mod_fcxb_factor / 1000**2 * np.sum(bgddetimsum[i])

        if i == refspec['A'] or i == refspec['B']:
            m.cutoffpl.norm.values = norm_preset
        else:
            m_ref = xspec.AllModels(refspec[fpm] + 1, model_name)
            m.cutoffpl.norm.link = '%e * %s:p%d' % (
                norm_preset / m_ref.cutoffpl.norm.values[0],
                model_name,
                3 * refspec[fpm] + 3)


def addmodel_intn(presets, refspec, bgddetimsum, model_num, model_name='intn'):
    """
    Xspec model component 5: intn (bknpower)

    The 3 parameters are read from preset and frozen (PhoIndx1, BreakE,
    PhoIndx2).

    norm is calculated as the detector mask area (bgddetimsum) times the
    multiplier for that detector, summed over all detectors.

    Inputs:

    presets - Preset model parameter values read from file.

    refspec - Dictionary with 'A' and 'B' keys whose values are the 0-based
    index of the reference spectrum for FPMA and FPMB.

    bgddetimsum - List of area inside background region for each CCD.

    model_num - Model component number in Xspec, cannot be the same as another
    model or it will replace.

    model_name - Model component name in Xspec, default 'intn', cannot be the
    same as another model.
    """
    mod_intn = presets['models'][3]['components']

    xspec.Model('bknpower', model_name, model_num)

    for i in range(xspec.AllData.nSpectra):
        spec = xspec.AllData(i + 1)
        fpm = util.fpm_parse(spec.fileinfo('INSTRUME'))
        m = xspec.AllModels(i + 1, model_name)

        m.bknpower.PhoIndx1.values = mod_intn['bknpower'][fpm]['phoindx1']
        m.bknpower.BreakE.values = mod_intn['bknpower'][fpm]['breake']
        m.bknpower.PhoIndx2.values = mod_intn['bknpower'][fpm]['phoindx2']
        m.bknpower.PhoIndx1.frozen = True
        m.bknpower.BreakE.frozen = True
        m.bknpower.PhoIndx2.frozen = True

        norm_preset = np.sum(
            bgddetimsum[i] * np.float64(mod_intn['bknpower'][fpm]['ifactors'])
        )

        if i == refspec['A'] or i == refspec['B']:
            m.bknpower.norm.values = norm_preset
        else:
            m_ref = xspec.AllModels(refspec[fpm] + 1, model_name)
            m.bknpower.norm.link = '%e * %s:p%d' % (
                norm_preset / m_ref.bknpower.norm.values[0],
                model_name,
                4 * refspec[fpm] + 4
            )


def run_fit():
    """
    Tell Xspec to perform the fitting using these settings.
    """
    xspec.AllData.ignore('**-3. 150.-**')
    xspec.Fit.method = 'leven 30000 1e-4'
    xspec.Fit.statMethod = 'chi'
    xspec.Fit.perform()


def save_xcm(prefix='bgdparams'):
    """
    Save current Xspec state to an xcm file.

    Input:

    prefix - Specify file name prefix for the saved file. By default this is
    'bgdparams', and if a file with this name exists, will attempt
    bgdparams2.xcm through bgdparams100.xcm.
    """
    i = 1
    while i < 101:  # Retry 100 times...
        if i > 1:
            savefile = '%s%d.xcm' % (prefix, i)
        else:
            savefile = '%s.xcm' % prefix
        if not os.path.exists(savefile):
            xspec.Xset.save(savefile, info='a')
            print('\n*** Saved results to %s. ***' % savefile)
            return True
            break
        else:
            i += 1
    return False
