# (c) 2019 Qian Wang
# This file is part of nuskybgd released under MIT License.
# See LICENSE file for full details.

import xspec
import os
import json
import numpy as np
from . import util
from . import conf


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

    for i in range(spec_count):
        spec_number = i + spec_start  # Spectrum number in Xspec
        spec_arrinx = spec_number - 1  # Index in bgdapimwt array etc.
        spec = xspec.AllData(spec_number)
        fpm = util.fpm_parse(spec.fileinfo('INSTRUME'))
        m = xspec.AllModels(spec_number, model_name)

        if spec_arrinx == refspec['A'] or spec_arrinx == refspec['B']:
            m.cutoffpl.PhoIndex.values = mod_apbgd['cutoffpl'][fpm]['phoindex']
            m.cutoffpl.HighECut.values = mod_apbgd['cutoffpl'][fpm]['highecut']
            m.cutoffpl.PhoIndex.frozen = True
            m.cutoffpl.HighECut.frozen = True
            m.cutoffpl.norm.values = 0.002353 / 32 * bgdapimwt[spec_arrinx]
        else:
            m.cutoffpl.PhoIndex.link = '%s:p%d' % (
                model_name, 3 * refspec[fpm] + 1)
            m.cutoffpl.HighECut.link = '%s:p%d' % (
                model_name, 3 * refspec[fpm] + 2)
            m.cutoffpl.norm.link = '%e * %s:p%d' % (
                bgdapimwt[spec_arrinx] / bgdapimwt[refspec[fpm]],
                model_name,
                3 * refspec[fpm] + 3)


def applymodel_intbgd(presets, refspec, bgddetimsum, model_num, src_number=None,
                    model_name='intbgd', fix_line_ratios=False):
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
            parameters for a newly added spectrum.

        model_name - (Optional) Model component name in Xspec, default 'intbgd',
            cannot be the same as another model.
    """
    mod_intbgd = presets['models'][1]['components']

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

    # Are we adding the model, or updating it?
    if model_num not in xspec.AllModels.sources:
        xspec.Model('apec' + '+lorentz' * (len(mod_intbgd) - 1),
                model_name, model_num)
    elif model_name != xspec.AllModels.sources[model_num]:
        # Model already exists when source spectrum is added. Proceed to adjust
        # the params.
        print('Error: the requested model number (%d, %s) exists and does not match '
              'specified model name (%s). Cannot proceed with updating parameters!' %(
              model_num, xspec.AllModels.sources[model_num], model_name))
        return False

    # There are in total (len(mod_intbgd) * 3 + 1) parameters per spectrum. 4
    # parameters for apec and 3 parameters for lorentz.
    for i in range(spec_count):
        spec_number = i + spec_start  # Spectrum number in Xspec
        spec_arrinx = spec_number - 1  # Index in bgdapimwt array etc.
        spec = xspec.AllData(spec_number)
        fpm = util.fpm_parse(spec.fileinfo('INSTRUME'))
        m = xspec.AllModels(spec_number, model_name)

        # Parameter tracker apec --> 4 params
        
        # Reference spectrum
        if spec_arrinx == refspec['A'] or spec_arrinx == refspec['B']:
            m_ref_npar_offset = m.nParameters * refspec[fpm]

            apec_norm = np.sum(
                               bgddetimsum[spec_arrinx] * np.array(
                                  mod_intbgd['apec'][fpm]['ifactors'])
                                ) / np.sum(bgddetimsum[spec_arrinx])
            par_dict = {1:f"{mod_intbgd['apec'][fpm]['kt']}, -0.1",
                          2:f"{mod_intbgd['apec'][fpm]['abundanc']}, -0.1",
                          3:f"{mod_intbgd['apec'][fpm]['redshift']}, -0.1",
                          4:f"{apec_norm}"}
            parnum =5   
            
#           The old way...this is slow because it evaluates the model every time
#           you set the parameters this way.
#
#
#             m.apec.kT.values = mod_intbgd['apec'][fpm]['kt']
#             m.apec.Abundanc.values = mod_intbgd['apec'][fpm]['abundanc']
#             m.apec.Redshift.values = mod_intbgd['apec'][fpm]['redshift']
#             m.apec.kT.frozen = True
#             m.apec.Abundanc.frozen = True
#             m.apec.Redshift.frozen = True
#             m.apec.norm.values = np.sum(
#                 bgddetimsum[spec_arrinx] * np.array(
#                     mod_intbgd['apec'][fpm]['ifactors'])
#             ) / np.sum(bgddetimsum[spec_arrinx])


            # 3 solar lines and 19 keV line -- load them from preset
            # Leave these free, zero them out if you don't want them.
            for attr in ['lorentz', 'lorentz_3', 'lorentz_4', 'lorentz_5']:
                line_norm = np.sum(
                    bgddetimsum[spec_arrinx] * np.array(
                        mod_intbgd[attr][fpm]['ifactors'])
                ) / np.sum(bgddetimsum[spec_arrinx])
                
                ref_norm = line_norm
                ref_par = int(parnum+2) + m_ref_npar_offset
                this_par = {parnum:f"{mod_intbgd[attr][fpm]['linee']}, -0.1",
                        parnum+1:f"{mod_intbgd[attr][fpm]['width']}, -0.1",
                        parnum+2:f"{line_norm}, 0.1"}
                par_dict.update(this_par)
                parnum += 3

                continue
                # Below is the old way to do this
#                 m_line = getattr(m, attr)
#                 m_line.LineE.values = mod_intbgd[attr][fpm]['linee']
#                 m_line.Width.values = mod_intbgd[attr][fpm]['width']
#                 m_line.LineE.frozen = True
#                 m_line.Width.frozen = True
#                 m_line.norm.values = np.sum(
#                     bgddetimsum[spec_arrinx] * np.array(
#                         mod_intbgd[attr][fpm]['ifactors'])
#                 ) / np.sum(bgddetimsum[spec_arrinx])
                

            # All the other lines -- lorentz_6 through lorentz_(components) --
            # if fix_line_ratios=True, scale their initial norm to 19 keV line
            # (lorentz_5) using preset ratios. There are 4+3*3=13 parameters
            # before lorentz_5.
            ref_set = False
            for attr_n in range(5, len(mod_intbgd)):
                attr = 'lorentz_%d' % (attr_n + 1)
                line_norm = np.sum(
                        bgddetimsum[spec_arrinx] * np.array(
                        mod_intbgd[attr][fpm]['ifactors'])
                ) / np.sum(bgddetimsum[spec_arrinx])
                
                if not ref_set:
                    # If this is set, then enabled the logic to make all lines relative
                    # to the reference line.
                    if fix_line_ratios:
                        ref_set = True
                    ref_norm = line_norm
                    ref_par = int(parnum+2)+ m_ref_npar_offset
                    this_par = {parnum:f"{mod_intbgd[attr][fpm]['linee']}, -0.1",
                            parnum+1:f"{mod_intbgd[attr][fpm]['width']}, -0.1",
                            parnum+2:f"{line_norm}, 0.1"}
                else:
                    rel_norm = line_norm / ref_norm                
                    this_par = {parnum:f"{mod_intbgd[attr][fpm]['linee']}, -0.1",
                                parnum+1:f"{mod_intbgd[attr][fpm]['width']}, -0.1",
                                parnum+2:f"={model_name}:{ref_par}*{rel_norm}"}
                par_dict.update(this_par)
                parnum += 3

                continue
                
                # Replaced by above
#                 attr = 'lorentz_%d' % (attr_n + 1)
#                 m_line = getattr(m, attr)
#                 m_line.LineE.values = mod_intbgd[attr][fpm]['linee']
#                 m_line.Width.values = mod_intbgd[attr][fpm]['width']
#                 m_line.LineE.frozen = True
#                 m_line.Width.frozen = True
# 
#                 norm_preset = np.sum(
#                     bgddetimsum[spec_arrinx] * np.array(
#                         mod_intbgd[attr][fpm]['ifactors'])
#                 ) / np.sum(bgddetimsum[spec_arrinx])
# 
#                 if fix_line_ratios:
#                     lorentz_5_norm_npar = spec_arrinx * m.nParameters + 16
# 
#                     m_line.norm.link = '%e * %s:p%d' % (
#                         norm_preset / m.lorentz_5.norm.values[0],
#                         model_name,
#                         lorentz_5_norm_npar)
#                 else:
#                     m_line.norm.values = norm_preset

            # Actually apply the changes here
            m.setPars(par_dict)


        else:
            # Link to ref spectrum
            m_ref = xspec.AllModels(refspec[fpm] + 1, model_name)
            m_ref_npar_offset = m.nParameters * refspec[fpm]

            # Link apec to refspec
            


            norm_preset = np.sum(
                bgddetimsum[spec_arrinx] * np.array(
                    mod_intbgd['apec'][fpm]['ifactors'])
            ) / np.sum(bgddetimsum[spec_arrinx])

            ##################
            # Special consideration if not startig with spectrum 1.
            # The reference preset may have changed after fitting, so we must
            # calculate its original value for scaling.
            ref_preset = np.sum(
                bgddetimsum[refspec[fpm]] * np.array(
                    mod_intbgd['apec'][fpm]['ifactors'])
            ) / np.sum(bgddetimsum[refspec[fpm]])
            ##################


#             m.apec.kT.link = '%s:p%d' % (
#                 model_name, m_ref_npar_offset + 1)
#             m.apec.Abundanc.link = '%s:p%d' % (
#                 model_name, m_ref_npar_offset + 2)
#             m.apec.Redshift.link = '%s:p%d' % (
#                 model_name, m_ref_npar_offset + 3)
#             m.apec.norm.link = '%e * %s:p%d' % (
#                 norm_preset / ref_preset,
#                 model_name,
#                 m_ref_npar_offset + 4
#             )
# 
            par_dict = {1:f"={model_name}:{m_ref_npar_offset+1}",
                        2:f"={model_name}:{m_ref_npar_offset+2}",
                        3:f"={model_name}:{m_ref_npar_offset+3}",
                        4:f"={norm_preset/ref_preset}*{model_name}:{m_ref_npar_offset+4}"}


            parnum=5
            # All lines --- load values from preset and link to refspec
            for attr_n in range(1, len(mod_intbgd)):
                if attr_n == 1:
                    attr = 'lorentz'  # special case
                else:
                    attr = 'lorentz_%d' % (attr_n + 1)


                ######################
                # Special consideration if not startig with spectrum 1.
                # Calculate the original value for scaling
                ########################
                ref_preset = np.sum(
                    bgddetimsum[refspec[fpm]] * np.array(
                        mod_intbgd[attr][fpm]['ifactors'])
                ) / np.sum(bgddetimsum[refspec[fpm]])
                ######################

                norm_preset = np.sum(
                    bgddetimsum[i] * np.array(
                        mod_intbgd[attr][fpm]['ifactors'])
                ) / np.sum(bgddetimsum[i])
                
#                m_line = getattr(m, attr)
#                 m_line.LineE.link = '%s:p%d' % (
#                     model_name,
#                     m_ref_npar_offset + 3 * (attr_n - 1) + 5)
#                 m_line.Width.link = '%s:p%d' % (
#                     model_name,
#                     m_ref_npar_offset + 3 * (attr_n - 1) + 6)
#                 m_line.norm.link = '%e * %s:p%d' % (
#                     norm_preset / ref_preset,
#                     model_name,
#                     m_ref_npar_offset + 3 * (attr_n - 1) + 7
#                 )

                this_par = {parnum:f"={model_name}:{parnum+m_ref_npar_offset}",
                            parnum+1:f"={model_name}:{parnum+1+m_ref_npar_offset}",
                            parnum+2:f"={model_name}:{parnum+2+m_ref_npar_offset} * {norm_preset/ref_preset}"}
                par_dict.update(this_par)
                parnum += 3

            m.setPars(par_dict)


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

    mod_fcxb_factor = 0.002353 * 1.5 * (2.45810736 / 3600 * 1000)**2

    for i in range(spec_count):
        spec_number = i + spec_start  # Spectrum number in Xspec
        spec_arrinx = spec_number - 1  # Index in bgdapimwt array etc.
        spec = xspec.AllData(spec_number)
        fpm = util.fpm_parse(spec.fileinfo('INSTRUME'))
        m = xspec.AllModels(spec_number, model_name)

        m.cutoffpl.PhoIndex.link = '%s:p%d' % (
            apbgd_name, 3 * refspec[fpm] + 1)
        m.cutoffpl.HighECut.link = '%s:p%d' % (
            apbgd_name, 3 * refspec[fpm] + 2)

        norm_preset = mod_fcxb_factor / 1000**2 * np.sum(bgddetimsum[spec_arrinx])

        m.cutoffpl.norm.values = norm_preset

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

    norm is calculated as the detector mask area (bgddetimsum) times the
    multiplier for that detector, summed over all detectors.

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

    for i in range(spec_count):
        spec_number = i + spec_start  # Spectrum number in Xspec
        spec_arrinx = spec_number - 1  # Index in bgdapimwt array etc.
        spec = xspec.AllData(spec_number)
        fpm = util.fpm_parse(spec.fileinfo('INSTRUME'))
        m = xspec.AllModels(spec_number, model_name)

        m.bknpower.PhoIndx1.values = mod_intn['bknpower'][fpm]['phoindx1']
        m.bknpower.BreakE.values = mod_intn['bknpower'][fpm]['breake']
        m.bknpower.PhoIndx2.values = mod_intn['bknpower'][fpm]['phoindx2']
        m.bknpower.PhoIndx1.frozen = True
        m.bknpower.BreakE.frozen = True
        m.bknpower.PhoIndx2.frozen = True

        norm_preset = np.sum(
            bgddetimsum[spec_arrinx] * np.float64(mod_intn['bknpower'][fpm]['ifactors'])
        )


        if spec_arrinx == refspec['A'] or spec_arrinx == refspec['B']:
            m.bknpower.norm.values = norm_preset
        else:
            ######################
            # Special consideration if not startig with spectrum 1.
            # Calculate the original value for scaling
            ########################
            ref_preset = np.sum(
                bgddetimsum[refspec[fpm]] * np.float64(mod_intn['bknpower'][fpm]['ifactors'])
            )
            #######################

            m_ref = xspec.AllModels(refspec[fpm] + 1, model_name)
            m.bknpower.norm.link = '%e * %s:p%d' % (
                norm_preset / ref_preset,
                model_name,
                4 * refspec[fpm] + 4
            )


def addmodel_grxe(presets, refspec, model_num, model_name='grxe'):
    """
    Galactic ridge X-ray emission model, placeholder for now.
    """
    # Add the response for this source.
    for i in range(xspec.AllData.nSpectra):
        s = xspec.AllData(i + 1)
        s.multiresponse[model_num - 1] = s.response.rmf
        s.multiresponse[model_num - 1].arf = '%s/be.arf' % conf._AUX_DIR


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
