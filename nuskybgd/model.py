# (c) 2019 Qian Wang
# This file is part of nuskybgd released under MIT License.
# See LICENSE file for full details.

import xspec
import os
import json
import numpy as np
from . import util
from . import conf


class ModelSource():
    def __init__(self, xspec_model):
        self.xsmod = xspec_model

    def __str__(self):
        mod = {}

        mod['name'] = self.xsmod.name
        mod['expression'] = self.xsmod.expression
        mod['componentNames'] = self.xsmod.componentNames
        mod['nParameters'] = self.xsmod.nParameters
        mod['startParIndex'] = self.xsmod.startParIndex
        mod['components'] = {}

        for component in self.xsmod.componentNames:
            mod['components'][component] = {}
            comp = self.xsmod.__getattribute__(component)

            mod['components'][component]['parameterNames'] = comp.parameterNames

            for param in comp.parameterNames:
                mod['components'][component][param] = (
                    comp.__getattribute__(param).values[0])

        return json.dumps(mod, indent=4)

    def zero_norm_pars(self):
        """
        Return a dictionary for xspec.AllModels.setPars() to update all
        normalization parameters (those that match the name "norm") in this
        model to zero.

        Note: In Xspec, parameters for different data groups for a single model
        source are numbered sequentially, such that if there were 12 parameters
        in the model, then 1-12 apply to data group 1, 13-24 apply to data group
        2, and so on. The PyXspec model object is for one data group; the total
        number of parameters is given by model.nParameters and the starting
        parameter number of the data group is given by model.startParIndex. When
        performing parameter updates using model.setPars(), the numbering always
        starts at 1 (do not account for preceding data groups). Similarly when
        using xspec.AllModels.setPars() parameters for each data group starts at
        1.

        In practice, this means that the same dictionary that lists the
        parameter numbers to set to zero can be used for all data groups of a
        model source, since they have the same set of parameters.

        Example:

        ```
        # Get information for the intbgd source using data group 2
        test = ModelSource(xspec.AllModels(2, 'intbgd'))
        update_pars = test.zero_norm_pars()
        # Use this to zero normalization for all data groups

        allpars = []
        for i_grp in range(1, xspec.AllData.nGroups + 1):
            m = xspec.AllModels(i_grp, 'intbgd')
            allpars.extend([m, update_pars])
        xspec.AllModels.setPars(*allpars)
        ```
        """
        pars = {}

        pnum = 1

        for comp_name in self.xsmod.componentNames:
            param_names = self.xsmod.__getattribute__(comp_name).parameterNames
            try:
                norm_inx = param_names.index('norm')
                pars.update({
                    pnum + norm_inx: '0'
                })
            except ValueError:
                print('Param "norm" missing for component %s' % comp_name)

            pnum += len(param_names)

        return pars


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


def check_spec_order(bgdinfo):
    """
    Check if the order of the loaded spectra files is the same as that in
    bgdinfo. Run this when loading an Xspec save file for the background model
    as a consistency check, because if the number or order of spectra files
    don't match, subsequent steps are invalid. Returns True if the loaded
    spectra and entries in bgdinfo are the same.

    Input:

    bgdinfo - Object loaded from the background info JSON file.
    """
    print('Checking for consistency between bgdinfo and loaded spectra...')
    if len(bgdinfo['bgfiles']) != xspec.AllData.nSpectra:
        print('Number of loaded spectra differs from bgdinfo entries!')
        return False

    problem = False
    for i in range(xspec.AllData.nSpectra):
        specfile = xspec.AllData(i + 1).fileName
        if specfile != bgdinfo['bgfiles'][i]:
            problem = True
            print('%s\tError: non-matching %s' % (specfile, bgdinfo['bgfiles'][i]))
        else:
            print('%s\t%s\tOK' % (specfile, bgdinfo['regfiles'][i]))

    if problem:
        return False
    else:
        return True


def load_bgdimgs(bgdinfo):
    """
    Retrieve and return the aspect-projected background maps.

    Input:

    bgdinfo - This should come from check_bgdinfofile().

    Output:

    (bgdapim, bgddetim) - Both dictionaries.

    bgdapim = {'A': numpy.ndarray, 'B': numpy.ndarray}
    bgddetim = {'A': [numpy.ndarray], 'B': [numpy.ndarray]}
    """
    import astropy.io.fits as pf

    bgdapfiles = bgdinfo['bgdapfiles']
    bgddetfiles = bgdinfo['bgddetfiles']

    bgdapim = {}
    bgdapim['A'] = pf.open(bgdapfiles['A'])[0].data
    bgdapim['B'] = pf.open(bgdapfiles['B'])[0].data

    bgddetim = {}
    bgddetim['A'] = [
        pf.open(bgddetfiles['A'][0])[0].data,
        pf.open(bgddetfiles['A'][1])[0].data,
        pf.open(bgddetfiles['A'][2])[0].data,
        pf.open(bgddetfiles['A'][3])[0].data
    ]
    bgddetim['B'] = [
        pf.open(bgddetfiles['B'][0])[0].data,
        pf.open(bgddetfiles['B'][1])[0].data,
        pf.open(bgddetfiles['B'][2])[0].data,
        pf.open(bgddetfiles['B'][3])[0].data
    ]

    return bgdapim, bgddetim


def get_keyword_specfiles(specfiles, keyword, ext=0):
    """
    Get values of specified keyword from a list of spectra files.

    Input:

    specfiles -- List of spectra file names.
    keyword -- Header keyword to check.
    ext -- Extension index or name.

    Output:

    List of keyword values.

    If either ext or keyword could not be found, that entry will be None. If
    an input spectrum file does not exist, an OSError bubbles up.
    """
    values = []

    for i in range(len(specfiles)):

        values.append(
            util.fits_checkkeyword(
                specfiles[i], keyword, ext=ext, silent=True))

    return values


def get_keyword_xspecdata(keyword):
    """
    Get values of specified keyword from currently loaded Xspec data files.

    Input:

    keyword -- Header keyword to check.

    Output:

    List of keyword values.

    If the specified keyword is not found, that entry will be None.
    """
    values = []

    for i in range(xspec.AllData.nSpectra):
        spec = xspec.AllData(i + 1)

        try:
            values.append(spec.fileinfo(keyword))

        except KeyError:
            print('Spectrum %s does not have the %s keyword.'
                  % (spec.fileName, keyword))
            values.append(None)

    return values


def get_refspec(instlist):
    """
    Determine the list indices of reference spectra files. The first entry of
    each module (inferred from the INSTRUME keyword) is designated the
    reference spectrum. If no spectrum from a module is encountered, it will
    have an entry None. Note that in xspec.AllData(i), the index i is
    1-based.

    Input:

    instlist -- List of INSTRUME keyword values. This can be retrieved from a
        list of spectra files using get_keyword_specfiles() or from the loaded
        xspec.AllData using get_keyword_xspecdata().

    Output:

    {'A': None, 'B': None} where None is replaced with the list index of the
    {'first spectrum encountered for each module.

    An Exception is raised if the focal plane module could not be determined
    for any of the entries in specfiles.
    """
    refspec = {'A': None, 'B': None}

    for i in range(len(instlist)):

        fpm = util.fpm_parse(instlist[i])

        if fpm is False:
            raise Exception(
                'Could not determine focal plane module for spectrum #%d.' %
                i)
        else:
            if fpm == 'A' and refspec['A'] is None:
                refspec['A'] = i
            elif fpm == 'B' and refspec['B'] is None:
                refspec['B'] = i

    return refspec


def calc_det_weights(detmasks, regmasks, fpmlist):
    """
    Calculate weights by detector area in the masks.

    Input:

    detmasks -- Detector masks {'A': [], 'B': []} with lists of images for
        each detector.
    regmasks -- List of region mask images.
    fpmlist -- List of modules corresponding to each region.

    Output:

    {'sum': [], 'fraction': [[]]} List of detector area in each region per
    detector; list of area fraction on each detector for every region.
    """
    res = {'sum': [], 'fraction': []}

    for reg, instr in zip(regmasks, fpmlist):
        fpm = util.fpm_parse(instr)
        det_areas = [np.sum(reg * detim) for detim in detmasks[fpm]]
        tot = np.sum(det_areas)
        det_frac = det_areas / tot
        res['sum'].append(det_areas)
        res['fraction'].append(det_frac)

    return res


def calc_ap_weights(apimgs, regmasks, fpmlist):
    """
    Calculate sum of aperture image in the masks.

    Input:

    apimgs -- Detector masks {'A': np.ndarray, 'B': np.ndarray} with the
        aperture image of each module.
    regmasks -- List of region mask images.
    fpmlist -- List of modules corresponding to each region.

    Output:

    {'sum': []} List of sum of aperture image in each region.
    """
    res = {'sum': []}

    for reg, instr in zip(regmasks, fpmlist):
        fpm = util.fpm_parse(instr)
        res['sum'].append(np.sum(reg * apimgs[fpm]))

    return res


def addspec(specfiles, fresh=True):
    """
    Add spectral data files in Xspec, each in their own data group. Afterward
    check the number of loaded spectra against length of input list. If not
    all spectra loaded, an Exception is raised.

    Inputs:

    specfiles - List of files to load.

    fresh - (Optional, default True) All existing spectra are cleared before
        adding. If set to False, new spectra files are appended.
    """
    if fresh:
        xspec.DataManager.clear(0)  # Clear any existing loaded data
        count_before = 0
    else:
        count_before = xspec.AllData.nSpectra

    # Load each spectrum as a new data group
    for i in range(len(specfiles)):
        xspec.AllData('{num}:{num} {file}'.format(
            num=i + 1 + count_before,  # Xspec spectrum numbering starts at 1
            file=specfiles[i]))

    if xspec.AllData.nSpectra != count_before + len(specfiles):
        raise Exception('Not all requested spectra loaded, cannot proceed!')


def applymodel_apbgd(presets, refspec, bgdapimwt, model_num, src_number=None,
                     model_name='apbgd'):
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

    src_number - (Optional) Position of first spectrum to process, if not
        starting at 1. This is useful for rescaling background model parameters
        for a newly added spectrum.

    model_name - (Optional) Model component name in Xspec, default 'apbgd', cannot be the
        same as another model.
    """
    mod_apbgd = presets['models'][0]['components']

    if isinstance(src_number, int):
        spec_start = src_number
    else:
        spec_start = 1

    if not (xspec.AllData.nSpectra >= spec_start > 0):  # Must be within 1...nspec
        raise Exception('Cannot apply model for spectrum: spectrum number out of range.')

    # Process all spectra from spec_start to the end
    spec_count = xspec.AllData.nSpectra - (spec_start - 1)

    # Add the response for this source
    for i in range(spec_count):
        s = xspec.AllData(i + spec_start)
        s.multiresponse[model_num - 1] = s.response.rmf
        s.multiresponse[model_num - 1].arf = '%s/be.arf' % conf._AUX_DIR

    if model_num not in xspec.AllModels.sources:  # Adding new, or updating?
        xspec.Model('cutoffpl', model_name, model_num)
    elif model_name != xspec.AllModels.sources[model_num]:
        print('Error: the requested model number exists and does not match'
              'specified model name. Cannot proceed with updating parameters!')
        return False

    # Parameters for cutoffpl: 3
    # PhoIndex HighECut norm

    allpars = []

    for i in range(spec_count):
        spec_number = i + spec_start  # Spectrum number in Xspec
        spec_arrinx = spec_number - 1  # Index in bgdapimwt array etc.
        spec = xspec.AllData(spec_number)
        fpm = util.fpm_parse(spec.fileinfo('INSTRUME'))
        m = xspec.AllModels(spec_number, model_name)

        newpar = {}

        _pars = mod_apbgd['cutoffpl'][fpm]

        if spec_arrinx == refspec['A'] or spec_arrinx == refspec['B']:

            norm = 0.002353 / 32 * bgdapimwt[spec_arrinx]

            newpar.update({
                1: f"{_pars['phoindex']}, -0.01",
                2: f"{_pars['highecut']}, -0.01",
                3: f"{norm:e}, 0.01"
            })

        else:

            refparoffset = 3 * refspec[fpm]

            ratio = bgdapimwt[spec_arrinx] / bgdapimwt[refspec[fpm]]

            newpar.update({
                1: f"={model_name}:{refparoffset+1}",
                2: f"={model_name}:{refparoffset+2}",
                3: f"={ratio:e} * {model_name}:{refparoffset+3}"
            })

        allpars.extend([m, newpar])

    xspec.AllModels.setPars(*allpars)


def applymodel_intbgd(presets, refspec, bgddetimsum, model_num,
                      src_number=None, model_name='intbgd',
                      fix_line_ratios=False):
    """
    Xspec model component 3: intbgd (instrument background)

    This model consists of an APEC component followed by many lorentz lines.
    APEC has 4 parameters (kT, Abundanc, Redshift, norm). Each lorentz
    component has 3 parameters (LineE, Width, norm);

    For the reference spectra (one from each telescope):

    *   APEC params are loaded from preset.

    *   Lines 1-3 (lorentz, lorentz_3, lorentz_4) are loaded from preset. The
        first two (3 and 6 keV) are probably solar related, and the third (10
        keV) is poorly fit.

    *   From Line 4 (19 keV, lorentz_5) onward, they are loaded from preset
        (detector area weighted sum of pre-determined values for each
        detector):

            norm = sum (ifactor * bgddetimsum)
                   ---------------------------
                        sum (bgddetimsum)

    *   If the option fix_line_ratios=True, Line 5 onward are linked to Line 4
        using:

            sum(ifactor * bgddetimsum)

    For the other spectra:

    *   All of the normalizations are linked to their reference spec
        counterparts.


    Note: due to perculiar model component numbering in XSPEC, as in the
    following,

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

        model_num - Model component number in Xspec, cannot be the same as
            another model or it will replace.

        src_number - (Optional) Position of first spectrum to process, if not
            starting at 1. This is useful for rescaling background model
            parameters for a newly added spectrum. (Default: None)

        model_name - (Optional) Model component name in Xspec. Must not be the
            same as another model. (Default: 'intbgd')

        fix_line_ratios - (Optional) Scale Gaussian line norms after the 19
            keV line to the 19 keV norm. (Default: False)
    """
    mod_intbgd = presets['models'][1]['components']

    # Starting spectrum number
    if isinstance(src_number, int):
        spec_start = src_number
    else:
        spec_start = 1

    # Validate the requested startspec number
    if not (xspec.AllData.nSpectra >= spec_start > 0):  # Range is 1...nspec
        raise Exception(
            'Cannot apply model for spectrum: spectrum number out of range.')

    # Number of spectra to process
    spec_count = xspec.AllData.nSpectra - (spec_start - 1)

    # Add the response for this source.
    for i in range(spec_count):
        s = xspec.AllData(i + spec_start)
        s.multiresponse[model_num - 1] = s.response.rmf

    # Are we adding the model, or updating it?
    if model_num not in xspec.AllModels.sources:
        xspec.Model('apec' + '+lorentz' * (len(mod_intbgd) - 1),
                    model_name, model_num)
    elif model_name != xspec.AllModels.sources[model_num]:
        # Model already exists when source spectrum is added. Proceed to adjust
        # the params.
        print('Error: the requested model number (%d, %s) exists and does not match '
              'specified model name (%s). Cannot proceed with updating parameters!' % (
                  model_num, xspec.AllModels.sources[model_num], model_name))
        return False

    # Parameters for apec:
    # kT Abundanc Redshift norm
    # Parameters for lorentz
    # LineE Width norm

    allpars = []

    # There are in total (len(mod_intbgd) * 3 + 1) parameters per spectrum. 4
    # parameters for apec and 3 parameters for lorentz.
    for i in range(spec_count):
        spec_number = i + spec_start  # Spectrum number in Xspec
        spec_arrinx = spec_number - 1  # Index in bgdapimwt array etc.
        spec = xspec.AllData(spec_number)
        fpm = util.fpm_parse(spec.fileinfo('INSTRUME'))
        m = xspec.AllModels(spec_number, model_name)

        newpar = {}

        # Reference spectrum
        if spec_arrinx == refspec['A'] or spec_arrinx == refspec['B']:
            m_ref_npar_offset = m.nParameters * refspec[fpm]

            _pars = mod_intbgd['apec'][fpm]

            apec_norm = np.sum(
                bgddetimsum[spec_arrinx] * np.array(_pars['ifactors'])
            ) / np.sum(bgddetimsum[spec_arrinx])

            newpar.update({
                1: f"{_pars['kt']}, -0.01",
                2: f"{_pars['abundanc']}, -0.001",
                3: f"{_pars['redshift']}, -0.01",
                4: f"{apec_norm:e}"
            })

            parnum = 5

            # 3 solar lines and 19 keV line -- load them from preset
            # Leave these free, zero them out if you don't want them.
            for attr in ['lorentz', 'lorentz_3', 'lorentz_4', 'lorentz_5']:

                _pars = mod_intbgd[attr][fpm]

                line_norm = np.sum(
                    bgddetimsum[spec_arrinx] * np.array(_pars['ifactors'])
                ) / np.sum(bgddetimsum[spec_arrinx])

                newpar.update({
                    parnum: f"{_pars['linee']}, -0.05",
                    parnum + 1: f"{_pars['width']}, -0.05",
                    parnum + 2: f"{line_norm}, 0.01"
                })

                parnum += 3

            # All the other lines -- lorentz_6 through lorentz_(components) --
            # if fix_line_ratios=True, scale their initial norm to 19 keV line
            # (lorentz_5) using preset ratios. There are 4+3*3=13 parameters
            # before lorentz_5.
            for attr_n in range(5, len(mod_intbgd)):
                attr = 'lorentz_%d' % (attr_n + 1)

                _pars = mod_intbgd[attr][fpm]

                norm_preset = np.sum(
                    bgddetimsum[spec_arrinx] * np.array(_pars['ifactors'])
                ) / np.sum(bgddetimsum[spec_arrinx])

                newpar.update({
                    parnum: f"{_pars['linee']}, -0.05",
                    parnum + 1: f"{_pars['width']}, -0.05"
                })

                if fix_line_ratios:
                    lorentz_5_norm_npar = spec_arrinx * m.nParameters + 16

                    line_ratio = norm_preset / m.lorentz_5.norm.values[0]

                    newpar.update({
                        parnum + 2: f"={line_ratio:e} * {model_name}:{lorentz_5_norm_npar}"
                    })

                else:
                    newpar.update({
                        parnum + 2: f"{norm_preset:e}, 0.01"
                    })

                parnum += 3

        else:
            # Link to ref spectrum
            # m_ref = xspec.AllModels(refspec[fpm] + 1, model_name)
            m_ref_npar_offset = m.nParameters * refspec[fpm]

            # Link apec to refspec

            _pars = mod_intbgd['apec'][fpm]

            norm_preset = np.sum(
                bgddetimsum[spec_arrinx] * np.array(_pars['ifactors'])
            ) / np.sum(bgddetimsum[spec_arrinx])

            ##################
            # Special consideration if not startig with spectrum 1.
            # The reference preset may have changed after fitting, so we must
            # calculate its original value for scaling.
            ref_preset = np.sum(
                bgddetimsum[refspec[fpm]] * np.array(_pars['ifactors'])
            ) / np.sum(bgddetimsum[refspec[fpm]])
            ##################

            line_ratio = norm_preset / ref_preset

            newpar.update({
                1: f"={model_name}:{m_ref_npar_offset+1}",
                2: f"={model_name}:{m_ref_npar_offset+2}",
                3: f"={model_name}:{m_ref_npar_offset+3}",
                4: f"={line_ratio:e} * {model_name}:{m_ref_npar_offset+4}"
            })

            parnum = 5

            # All lines --- load values from preset and link to refspec
            for attr_n in range(1, len(mod_intbgd)):
                if attr_n == 1:
                    attr = 'lorentz'  # special case
                else:
                    attr = 'lorentz_%d' % (attr_n + 1)

                _pars = mod_intbgd[attr][fpm]

                ######################
                # Special consideration if not startig with spectrum 1.
                # Calculate the original value for scaling
                ########################
                ref_preset = np.sum(
                    bgddetimsum[refspec[fpm]] * np.array(
                        _pars['ifactors'])
                ) / np.sum(bgddetimsum[refspec[fpm]])
                ######################

                norm_preset = np.sum(
                    bgddetimsum[i] * np.array(
                        _pars['ifactors'])
                ) / np.sum(bgddetimsum[i])

                line_ratio = norm_preset / ref_preset

                newpar.update({
                    parnum: f"={model_name}:{parnum+m_ref_npar_offset}",
                    parnum + 1: f"={model_name}:{parnum+1+m_ref_npar_offset}",
                    parnum + 2: f"={line_ratio:e} * {model_name}:{parnum+2+m_ref_npar_offset}"
                })

                parnum += 3

        allpars.extend([m, newpar])

    xspec.AllModels.setPars(*allpars)


def applymodel_fcxb(refspec, bgddetimsum, model_num, src_number=None,
                  model_name='fcxb', apbgd_name='apbgd'):
    """
    Xspec model component 4: fcxb (focused CXB)

    Dependency: the aperture image background model must have already been
    added.

    cutoffpl, with first two parameters (PhoIndex, HighCut) tied to apbgd model
    source.

    The norm of cutoffpl is to be computed from

    0.002353 * 1.5 * (2.45810736/3600*1000)^2 * backscal

    backscal = np.sum(regmask * bgddetim) / 1000^2

    Inputs:

    refspec - Dictionary with 'A' and 'B' keys whose values are the 0-based
        index of the reference spectrum for FPMA and FPMB.

    bgddetimsum - List of area inside background region for each CCD.

    model_num - Model component number in Xspec, cannot be the same as another
        model or it will replace.

    model_name - Model component name in Xspec, default 'fcxb', cannot be the
        same as another model.

    src_number - (Optional) Position of first spectrum to process, if not
        starting at 1. This is useful for rescaling background model parameters
        for a newly added spectrum.

    apbgd_name - Model component name of the aperture background model, must
        have been already added.
    """
    # Starting spectrum number
    if isinstance(src_number, int):
        spec_start = src_number
    else:
        spec_start = 1

    # Validate the requested startspec number
    if not (xspec.AllData.nSpectra >= spec_start > 0):  # Must be within 1...nspec
        raise Exception('Cannot apply model for spectrum: spectrum number out of range.')

    # Number of spectra to process
    spec_count = xspec.AllData.nSpectra - (spec_start - 1)

    # Add the response for this source.
    for i in range(spec_count):
        s = xspec.AllData(i + spec_start)
        s.multiresponse[model_num - 1] = s.response.rmf
        s.multiresponse[model_num - 1].arf = '%s/fcxb%s.arf' % (
            conf._AUX_DIR, util.fpm_parse(s.fileinfo('INSTRUME')))

    # Are we adding the model, or updating it?
    if model_num not in xspec.AllModels.sources:
        xspec.Model('cutoffpl', model_name, model_num)
    elif model_name != xspec.AllModels.sources[model_num]:
        # Model already exists when source spectrum is added. Proceed to adjust the params.
        print('Error: the requested model number exists and does not match'
              'specified model name. Cannot proceed with updating parameters!')
        return False

    # Parameters for cutoffpl: 3
    # PhoIndex HighECut norm

    mod_fcxb_factor = 0.002353 * 1.5 * (2.45810736 / 3600 * 1000)**2

    allpars = []

    for i in range(spec_count):
        spec_number = i + spec_start  # Spectrum number in Xspec
        spec_arrinx = spec_number - 1  # Index in bgdapimwt array etc.
        spec = xspec.AllData(spec_number)
        fpm = util.fpm_parse(spec.fileinfo('INSTRUME'))
        m = xspec.AllModels(spec_number, model_name)

        refparoffset = 3 * refspec[fpm]

        norm_preset = mod_fcxb_factor / 1000**2 * np.sum(bgddetimsum[spec_arrinx])

        newpar = {
            1: f"={apbgd_name}:{refparoffset+1}",
            2: f"={apbgd_name}:{refparoffset+2}",
            3: f"{norm_preset:e}, 0.01"
        }

        allpars.extend([m, newpar])

        # if spec_arrinx == refspec['A'] or spec_arrinx == refspec['B']:
        #     m.cutoffpl.norm.values = norm_preset
        # else:
        #     m_ref = xspec.AllModels(refspec[fpm] + 1, model_name)

        #     ######################
        #     # Special consideration if not startig with spectrum 1.
        #     # Calculate the original value for scaling
        #     ########################
        #     ref_preset = mod_fcxb_factor / 1000**2 * np.sum(bgddetimsum[refspec[fpm]])
        #     #######################

        #     m.cutoffpl.norm.link = '%e * %s:p%d' % (
        #         norm_preset / ref_preset,
        #         model_name,
        #         3 * refspec[fpm] + 3)

    xspec.AllModels.setPars(*allpars)


def fcxb_linkab(links, model_name='fcxb'):
    """
    Xspec model component 4: fcxb (focused CXB)

    Tie normalizations between A and B regions that cover the same region of
    sky.

    Inputs:

    links - List of arrays, [ispec1, ispec2], to link Xspec spectrum number 2
    to spectrum number 1.

    model_name - Model component name in Xspec, default 'fcxb', cannot be the
        same as another model.
    """
    npars = xspec.AllModels(1, model_name).nParameters

    for i in range(len(links)):
        # These are expected to be Xspec spectrum numbers (first one is 1)
        link_ref = links[i][0]
        link_tied = links[i][1]
        norm_offset = npars * (link_ref - 1)
        m = xspec.AllModels(link_tied, model_name)
        m.cutoffpl.norm.link = '%s:p%d' % (
            model_name,
            norm_offset + 3
            )


def applymodel_intn(presets, refspec, bgddetimsum, model_num, src_number=None,
                  model_name='intn'):
    """
    Xspec model component 5: intn (bknpower)

    The 3 parameters are read from preset and frozen (PhoIndx1, BreakE,
    PhoIndx2).

    Normalization is calculated by multiplying the detector mask area
    (bgddetimsum) in arcmin^2 (1 px is 2.45 arcsec) and the preset for that
    detector, summed over all detectors.

    Inputs:

    presets - Preset model parameter values read from file.

    refspec - Dictionary with 'A' and 'B' keys whose values are the 0-based
        index of the reference spectrum for FPMA and FPMB.

    bgddetimsum - List of area inside background region for each CCD.

    model_num - Model component number in Xspec, cannot be the same as another
        model or it will replace.

    src_number - (Optional) Position of first spectrum to process, if not
        starting at 1. This is useful for rescaling background model parameters
        for a newly added spectrum.

    model_name - Model component name in Xspec, default 'intn', cannot be the
        same as another model.
    """
    mod_intn = presets['models'][3]['components']

    # Starting spectrum number
    if isinstance(src_number, int):
        spec_start = src_number
    else:
        spec_start = 1

    # Validate the requested startspec number
    if not (xspec.AllData.nSpectra >= spec_start > 0):  # Must be within 1...nspec
        raise Exception('Cannot apply model for spectrum: spectrum number out of range.')

    # Number of spectra to process
    spec_count = xspec.AllData.nSpectra - (spec_start - 1)

    # Add the response for this source.
    for i in range(spec_count):
        s = xspec.AllData(i + spec_start)
        s.multiresponse[model_num - 1] = '%s/diag.rmf' % conf._AUX_DIR

    # Are we adding the model, or updating it?
    if model_num not in xspec.AllModels.sources:
        xspec.Model('bknpower', model_name, model_num)
    elif model_name != xspec.AllModels.sources[model_num]:
        # Model already exists when source spectrum is added. Proceed to adjust the params.
        print('Error: the requested model number exists and does not match'
              'specified model name. Cannot proceed with updating parameters!')
        return False

    # Parameters for bknpower:
    # PhoIndx1 BreakE PhoIndx2 norm

    allpars = []

    for i in range(spec_count):
        spec_number = i + spec_start  # Spectrum number in Xspec
        spec_arrinx = spec_number - 1  # Index in bgdapimwt array etc.
        spec = xspec.AllData(spec_number)
        fpm = util.fpm_parse(spec.fileinfo('INSTRUME'))
        m = xspec.AllModels(spec_number, model_name)

        _pars = mod_intn['bknpower'][fpm]

        norm_preset = np.sum(
            bgddetimsum[spec_arrinx] * np.float64(_pars['ifactors'])
        ) / 599.75  # 1 arcsec^2 = 599.75 px

        newpar = {}

        newpar.update({
            1: f"{_pars['phoindx1']}, -0.01",
            2: f"{_pars['breake']}, -0.01",
            3: f"{_pars['phoindx2']}, -0.01"
        })

        if spec_arrinx == refspec['A'] or spec_arrinx == refspec['B']:
            newpar.update({
                4: f"{norm_preset:e}, 0.01"
            })
        else:
            ######################
            # Special consideration if not startig with spectrum 1.
            # Calculate the original value for scaling
            ########################
            ref_preset = np.sum(
                bgddetimsum[refspec[fpm]] * np.float64(_pars['ifactors'])
            ) / 599.75  # 1 arcsec^2 = 599.75 px
            #######################

            ratio = norm_preset / ref_preset

            newpar.update({
                4: f"={ratio:e} * {model_name}:{4 * refspec[fpm] + 4}"
            })

        allpars.extend([m, newpar])

    xspec.AllModels.setPars(*allpars)


def addmodel_grxe(presets, refspec, model_num, model_name='grxe'):
    """
    Galactic ridge X-ray emission model, placeholder for now.
    """
    # Add the response for this source.
    for i in range(xspec.AllData.nSpectra):
        s = xspec.AllData(i + 1)
        s.multiresponse[model_num - 1] = s.response.rmf
        s.multiresponse[model_num - 1].arf = '%s/be.arf' % conf._AUX_DIR


def read_xspec_model_norms(spec_nums=None):
    """
    Make a record of the normalizations of the current XSPEC model. Only the
    values are saved, not parameter link information. (Note that if
    normalization value is set directly, the parameter link is overwritten.)

    Input:

    spec_nums -- (Optional) An array of spectrum numbers (starts at 1). If not
        specified, all spectra are included.

    Output:

    A dictionary with the format
    {spec_num: {
        'source_name':{
            'component1': norm1,
            'component2': norm2,
            ...},...},...}
    """
    model_norms = {}

    print('Recording current XSPEC model normalizations...')
    print('=' * 80)

    if spec_nums is None:
        # Include all spectra
        spec_nums = range(1, 1 + xspec.AllData.nSpectra)

    for spec_num in spec_nums:

        model_norms[spec_num] = {}

        for model_source in xspec.AllModels.sources.items():
            print('-' * 80)
            print('Spectrum number, model source: ', spec_num, model_source)

            model_norms[spec_num][model_source[1]] = {}
            current_model = xspec.AllModels(spec_num, model_source[1])
            source_components = current_model.componentNames

            print('Source components: ', source_components)
            print('Listing component_name, norm: ')

            for component_name in source_components:
                current_norm = current_model.__getattribute__(component_name).norm.values[0]
                model_norms[spec_num][model_source[1]][component_name] = current_norm
                print(component_name, current_norm)

            print('-' * 80)
        print('End of spectrum number', spec_num)
        print('=' * 80)
    return model_norms


def zero_xspec_model_norms(model_norms):
    """
    Iterate through sources and model components in model_norms (obtained by
    read_xspec_model_norms) and set all of the normalizations to zero.

    Note: setting the values to zero breaks any link it has to another
    parameter. The model thus modified should not be used for fitting.
    """
    print('Setting XSPEC model normalizations to zero...')
    print('=' * 80)

    allpars = []

    for _, source_name in xspec.AllModels.sources.items():
        print('Source: %s' % source_name)
        # Use first data group to get the list of parameter numbers to zero
        source_pars = ModelSource(xspec.AllModels(1, source_name)).zero_norm_pars()

        # Apply to every data group by adding them to allpars.
        # Repeated for all sources.
        for grp_num in range(1, xspec.AllData.nGroups + 1):
            print('Data group %d' % grp_num)
            m = xspec.AllModels(grp_num, source_name)
            allpars.extend([m, source_pars])
        print('-' * 80)

    xspec.AllModels.setPars(*allpars)

    # for spec_num, sources in model_norms.items():
    #     for source_name, source_components in sources.items():
    #         print('-' * 80)
    #         print('Current spectrum number, model source: ', spec_num, source_name)
    #         current_model = xspec.AllModels(spec_num, source_name)
    #         for component_name, _ in source_components.items():
    #             current_model.__getattribute__(component_name).norm.values = 0.0
    #             print('Set to zero: ', component_name)
    #         print('-' * 80)
    #     print('End of spectrum ', spec_num)
    #     print('=' * 80)


def bgimg_apbgd(bgdinfo, model_norms, bgdapweights, refspec,
                model_name='apbgd', ignore='**-3. 20.-**',
                outprefix=''):
    """
    Create two images, one for each module, based on the aperture images.

    The spectral models were constructed by summing the aperture images in the
    background regions, and using that as the weight to scale different
    background regions to the reference spectra (the first spectrum for each
    of A and B).

    The image model of the background is created through the following steps:

    - 1. For each FPM reference region, sum the aperture images in the
      background regions;
    - 2. Ask XSPEC for the model predicted counts of the source `apbgd` (in
      user-specified energy band);
    - 3. Multiply the aperture image by `model_count_rate / bgdap_region_sum`.

    Inputs:

    bgdinfo - Object returned by model.check_bgdinfofile().

    model_norms - Object returned by read_xspec_model_norms().

    bgdapweights - Object returned by model.calc_ap_weights() containing
        aperture image sums on each detector, for each spectrum.

    refspec - A dictionary indicating the reference spectra index,
        e.g. {'A': 0, 'B': 1}.

    model_name - (Optional) Source name used in XSPEC. Default: 'apbgd'.

    ignore - (Optional) The XSPEC 'ignore' command to use before asking for
        the model predicted count rate. Default: '**-3. 20.-**'.

    outprefix - (Optional) Prefix for the output file name.

    Output:

    True, upon completion.
    """
    import astropy.io.fits as pf

    # Set the energy interval
    xspec.AllData.notice('**')
    xspec.AllData.ignore(ignore)

    for fpm in ('A', 'B'):
        output = pf.open(bgdinfo['bgdapfiles'][fpm])

        spec_num = 1 + refspec[fpm]
        current_model = xspec.AllModels(spec_num, model_name)
        current_spec = xspec.AllData(spec_num)
        current_norms = model_norms[spec_num][model_name]

        # Set model norms for the reference spectra
        component_name = 'cutoffpl'
        current_model.__getattribute__(component_name).norm.values = current_norms[component_name]
        print('Set ', component_name, current_norms[component_name])

        total_counts = current_spec.rate[3] * current_spec.exposure
        output[0].data *= (total_counts / bgdapweights['sum'][refspec[fpm]])
        print(total_counts)

        output.writeto('%sbgd_%s_%s.fits' % (outprefix, model_name, fpm),
            overwrite=True)

        # Unset the model norms
        xspec.AllModels(1 + refspec[fpm], 'apbgd').cutoffpl.norm.values = 0.0

    return True


def bgimg_intbgd(presets, refspec, bgdinfo, model_norms, bgddetweights, bgddetim,
                 model_name='intbgd', ignore='**-3. 20.-**',
                 outprefix=''):
    """
    The model assigns constant background values for each detector (4
    detectors x 2 modules).

    In the spectral model, there are optional parameter links among the
    spectral components, but the detector weighting (the 'ifactors' in the
    auxil file) are fixed. However, the different spectral components have
    different sets of 'ifactors' weights, so for any given region, the
    contributions from the detectors are different for each component.

    The image model of the background is created through the following steps:

    - 1. For each module, go through each spectral model component separately
      (cache the model norms, then isolate the current component and zero the
      others);
    - 2. Ask XSPEC for the model count rate of `intbgd` (in user-specified
      energy band);
    - 3. Distribute the model count rate based on weights `ifactor *
      det_area`:
        - Because the region may not cover all detectors (i.e. some `det_area`
          is zero), we shall calculate the flux on the detector that has the
          most counts and then scale it for the other detectors using their
          `ifactor` values.
        - `model_count_rate * ifactor_m / sum_i(ifactor_i * det_area_i)`,
          where `ifactor_m * det_area_m` is the largest weight, will be the
          pixel value of `det_m`.
        - Calculate the pixel values of the other detectors using their
          relative `ifactor` values.
    - 4. Add the background image of the current component to the result.

    Inputs:

    presets - Preset information from auxil/ratios.json.

    refspec - A dictionary indicating the reference spectra index,
        e.g. {'A': 0, 'B': 1}.

    bgdinfo - Object returned by model.check_bgdinfofile().

    model_norms - Object returned by read_xspec_model_norms().

    bgddetweights - Object returned by model.calc_det_weights() containing
        detector image sums on each detector, for each spectrum.

    bgddetim - Object containing seperate detector images returned by
        model.load_bgdimgs().

    model_name - (Optional) Source name used in XSPEC. Default: 'intbgd'.

    ignore - (Optional) The XSPEC 'ignore' command to use before asking for
        the model predicted count rate. Default: '**-3. 20.-**'.

    outprefix - (Optional) Prefix for the output file name.

    Output:

    True, upon completion.
    """
    import astropy.io.fits as pf

    mod_intbgd = presets['models'][1]['components']

    # Set the energy interval
    xspec.AllData.notice('**')
    xspec.AllData.ignore(ignore)

    for fpm in ('A', 'B'):

        intbgd_values = [0.] * 4  # Tally intbgd counts / px on each detector

        # Use the detector with greatest area as reference for calculations
        det_areas = bgddetweights['sum'][refspec[fpm]]
        det_fracs = bgddetweights['fraction'][refspec[fpm]]
        refdet = det_areas.index(max(det_areas))

        spec_num = 1 + refspec[fpm]
        current_model = xspec.AllModels(spec_num, model_name)
        current_spec = xspec.AllData(spec_num)
        current_norms = model_norms[spec_num][model_name]

        # Check and warn about model components mismatch with presets
        if not (current_norms.keys() == mod_intbgd.keys()):
            print('Warning: XSPEC intbgd components do not match preset!')
            print('Preset:')
            print(', '.join(mod_intbgd.keys()))
            print('Loaded model:')
            print(', '.join(current_norms.keys()))

        for component_name in mod_intbgd.keys():

            # Set current component to fitted norm value
            current_model.__getattribute__(component_name).norm.values = current_norms[component_name]
            print('Set ', component_name, current_norms[component_name])

            total_counts = current_spec.rate[3] * current_spec.exposure

            # Zero current component norm
            current_model.__getattribute__(component_name).norm.values = 0.0

            # Determine intbgd counts / px on the reference detector
            ifactors = np.array(mod_intbgd[component_name][fpm]['ifactors'])
            refdet_bgd = total_counts * ifactors[refdet] / np.sum(ifactors * det_areas)

            # Use ifactors to scale to all detectors
            for i_det in range(4):
                intbgd_values[i_det] += refdet_bgd * ifactors[i_det] / ifactors[refdet]

        # Multiply detector images by their intbgd_values[i] and create composite image
        output = pf.open(bgdinfo['refimgf'])
        output[0].data *= 0.0

        for i_det in range(4):
            output[0].data += bgddetim[fpm][i_det] * intbgd_values[i_det]

        output.writeto('%sbgd_%s_%s.fits' % (outprefix, model_name, fpm),
            overwrite=True)

    return True


def bgimg_fcxb(bgdinfo, model_norms, bgddetweights, bgddetim, regmask,
               model_name='fcxb', ignore='**-3. 20.-**',
               outprefix=''):
    """
    All of the regions are allowed to have their own flux values.

    The image model of the background is created through the following steps:

    - 1. For each spectrum, ask XSPEC for the model count rate of `fcxb` (in
      user-specified energy band);
    - 2. Calculate the flux as `model_count_rate / n_pixel`.
    - 3. Fill the detector mask pixels in this spectral region with this flux
      value, output an image for this region;
    - 4. After all spectral regions have been processed, create an image for
      each module based on the fluxes in all of the background regions.

    The last step above currently produces an image assuming a constant
    background over the field of view. Users may make use of the images of
    individual regions to create their own models of this.

    Inputs:

    bgdinfo - Object returned by model.check_bgdinfofile().

    model_norms - Object returned by read_xspec_model_norms().

    bgddetweights - Object returned by model.calc_det_weights() containing
        detector image sums on each detector, for each spectrum.

    bgddetim - Object containing seperate detector images returned by
        model.load_bgdimgs().

    regmask - Image masks for the spectral regions returned by
        util.mask_from_region().

    model_name - (Optional) Source name used in XSPEC. Default: 'fcxb'.

    ignore - (Optional) The XSPEC 'ignore' command to use before asking for
        the model predicted count rate. Default: '**-3. 20.-**'.

    outprefix - (Optional) Prefix for the output file name.

    Output:

    True, upon completion.
    """
    import astropy.io.fits as pf

    fpmlist = [util.fpm_parse(_) for _ in get_keyword_xspecdata('INSTRUME')]

    # Set the energy interval
    xspec.AllData.notice('**')
    xspec.AllData.ignore(ignore)

    # Track simple mean flux
    overall_counts = 0.
    overall_area = 0.

    for i_spec in range(xspec.AllData.nSpectra):

        fpm = fpmlist[i_spec]

        # Use the detector with greatest area as reference for calculations
        det_areas = bgddetweights['sum'][i_spec]

        spec_num = 1 + i_spec
        current_model = xspec.AllModels(spec_num, model_name)
        current_spec = xspec.AllData(spec_num)
        current_norms = model_norms[spec_num][model_name]


        # Set current component to fitted norm value
        component_name = 'cutoffpl'
        current_model.__getattribute__(component_name).norm.values = current_norms[component_name]
        print('Set ', component_name, current_norms[component_name])

        total_counts = current_spec.rate[3] * current_spec.exposure
        overall_counts += total_counts
        print(total_counts)

        # Zero current component norm
        current_model.__getattribute__(component_name).norm.values = 0.0

        # Determine fcxb counts / px
        total_area = np.sum(det_areas)
        overall_area += total_area
        fcxb_flux = total_counts / total_area

        output = pf.open(bgdinfo['refimgf'])
        output[0].data *= 0.0

        for i_det in range(4):
            output[0].data += bgddetim[fpm][i_det] * fcxb_flux

        output[0].data *= regmask[i_spec]

        output.writeto('bgd_fcxb_%d.fits' % (1 + i_spec), overwrite=True)


    # Create detector images filled by mean flux across all spectral regions
    # (overlapping areas are treated as independent)
    overall_flux = overall_counts / overall_area
    for fpm in ('A', 'B'):
        output = pf.open(bgdinfo['refimgf'])
        output[0].data *= 0.0
        for i_det in range(4):
            output[0].data += bgddetim[fpm][i_det] * overall_flux
        output.writeto('%sbgd_%s_%s.fits' % (outprefix, model_name, fpm),
            overwrite=True)

    return True


def bgimg_intn(presets, refspec, bgdinfo, model_norms, bgddetweights, bgddetim,
               model_name='intn', ignore='**-3. 20.-**',
               outprefix=''):
    """
    The model uses a single component, a broken power law.

    There are only two free normalizations, one for each module. Relative
    weights of detectors are fixed by preset values of `ifactors`. The model
    image is created through an identical process to `intbgd` except we need
    only look at one model component.

    Inputs:

    presets - Preset information from auxil/ratios.json.

    refspec - A dictionary indicating the reference spectra index,
        e.g. {'A': 0, 'B': 1}.

    bgdinfo - Object returned by model.check_bgdinfofile().

    model_norms - Object returned by read_xspec_model_norms().

    bgddetweights - Object returned by model.calc_det_weights() containing
        detector image sums on each detector, for each spectrum.

    bgddetim - Object containing seperate detector images returned by
        model.load_bgdimgs().

    model_name - (Optional) Source name used in XSPEC. Default: 'intn'.

    ignore - (Optional) The XSPEC 'ignore' command to use before asking for
        the model predicted count rate. Default: '**-3. 20.-**'.

    outprefix - (Optional) Prefix for the output file name.

    Output:

    True, upon completion.
    """
    import astropy.io.fits as pf

    mod_intn = presets['models'][3]['components']

    # Set the energy interval
    xspec.AllData.notice('**')
    xspec.AllData.ignore(ignore)

    for fpm in ('A', 'B'):

        intbgd_values = [0.] * 4  # Tally intbgd counts / px on each detector

        # Use the detector with greatest area as reference for calculations
        det_areas = bgddetweights['sum'][refspec[fpm]]
        det_fracs = bgddetweights['fraction'][refspec[fpm]]
        refdet = det_areas.index(max(det_areas))

        spec_num = 1 + refspec[fpm]
        current_model = xspec.AllModels(spec_num, model_name)
        current_spec = xspec.AllData(spec_num)
        current_norms = model_norms[spec_num][model_name]

        # Check and warn about model components mismatch with presets
        if not (current_norms.keys() == mod_intn.keys()):
            print('Warning: the loaded XSPEC model for intn does not match components in preset!')
            print('Preset:')
            print(', '.join(mod_intn.keys()))
            print('Loaded model:')
            print(', '.join(current_norms.keys()))

        for component_name in mod_intn.keys():

            # Set current component to fitted norm value
            current_model.__getattribute__(component_name).norm.values = current_norms[component_name]
            print('Set ', component_name, current_norms[component_name])

            total_counts = current_spec.rate[3] * current_spec.exposure

            # Zero current component norm
            current_model.__getattribute__(component_name).norm.values = 0.0

            # Determine intbgd counts / px on the reference detector
            ifactors = np.array(mod_intn[component_name][fpm]['ifactors'])
            refdet_bgd = total_counts * ifactors[refdet] / np.sum(ifactors * det_areas)

            # Use ifactors to scale to all detectors
            for i_det in range(4):
                intbgd_values[i_det] += refdet_bgd * ifactors[i_det] / ifactors[refdet]

        # Multiply detector images by their intbgd_values[i] and create composite image
        output = pf.open(bgdinfo['refimgf'])
        output[0].data *= 0.0

        for i_det in range(4):
            output[0].data += bgddetim[fpm][i_det] * intbgd_values[i_det]

        output.writeto('%sbgd_%s_%s.fits' % (outprefix, model_name, fpm),
            overwrite=True)

    return True


def remove_ispec(nspec):
    """
    Remove a number of spectra from the front of the loaded spectra.

    Example:

    remove_ispec(nbgd)

        Removes the first nbgd spectra.

    Input:

    nspec - Number of spectra to remove.
    """
    if nspec > xspec.data.AllData.nSpectra:
        print('Error: cannot remove more spectra than there exists.')
    for i in range(nspec):
        xspec.data.AllData -= 1  # Remove the first spectrum many times


def freeze_pars(specinx):
    """
    Freeze spectrum model parameters for all models. Returns a list of parameters
    that are not linked and were not frozen, which can be used to thaw them.

    Example:

    frozenpars = freeze_pars(range(1, nbgd+1))

        All model parameters for spectra 1 through nbgd are frozen.

    Input:

    specinx - An int or list of int, the spectrum number(s) to freeze model
        parameters for.
    """
    frozen_pars = []

    if isinstance(specinx, int):
        specinx = [specinx]

    for modelnum, modelname in xspec.AllModels.sources.items():
        for j in specinx:
            for i in range(1, 1+xspec.AllModels(j, modelname).nParameters):
                if xspec.AllModels(j, modelname)(i).link == '':
                    if xspec.AllModels(j, modelname)(i).frozen is False:
                        frozen_pars.append((j, modelname, i))
                        xspec.AllModels(j, modelname)(i).frozen = True
    return frozen_pars


def thaw_pars(parlist):
    """
    Thaw a list of parameters. Input should be a list of tuples containing
    (spectrum number, model name, parameters number).

    Example:

    thaw_pars(frozenpars)

        Thaws all the parameters in the frozenpars list.

    Input:

    parlist - A list of tuples(3) containing (spectrum number, model name,
        parameter number).
    """
    for j, modname, i in parlist:
        try:
            xspec.AllModels(j, modname)(i).frozen = False
        except Exception:
            print('Error at ', j, modname, i)


def run_fit(statmethod='chi', method='leven 30000 1e-4', ignore='**-3. 150.-**'):
    """
    Change the fit settings in Xspec, and run fit.

    Inputs:

    statmethod - (Optional) The statistical method to be used by Xspec. This
        is passed directly to Xspec via an analogous command to 'statistic'.
        By default this is 'chi', requiring spectra to be grouped into bins
        with typically 30 counts at least.

    method - (Optional) Optimization method for the fit. By default this is
        'leven 30000 1e-4'. Any valid argument for the Xspec command 'fit' is
        acceptable.

    ignore - (Optional) Channels or energies to ignore. By default, this is
        '**-3. 150.-**'. Any valid argument for the Xspec command 'ignore'
        will work here.
    """
    run_fit_settings(statmethod=statmethod, method=method, ignore=ignore)
    xspec.Fit.perform()


def run_fit_settings(
    statmethod='chi', method='leven 30000 1e-4', ignore='**-3. 150.-**'):
    """
    Change the fit settings in Xspec.

    Inputs:

    statmethod - (Optional) The statistical method to be used by Xspec. This
        is passed directly to Xspec via an analogous command to 'statistic'.
        By default this is 'chi', requiring spectra to be grouped into bins
        with typically 30 counts at least.

    method - (Optional) Optimization method for the fit. By default this is
        'leven 30000 1e-4'. Any valid argument for the Xspec command 'fit' is
        acceptable.

    ignore - (Optional) Channels or energies to ignore. By default, this is
        '**-3. 150.-**'. Any valid argument for the Xspec command 'ignore'
        will work here.
    """
    xspec.AllData.ignore(ignore)
    xspec.Fit.method = method
    xspec.Fit.statMethod = statmethod


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
