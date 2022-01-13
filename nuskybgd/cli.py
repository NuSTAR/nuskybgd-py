# (c) 2019 Qian Wang
# This file is part of nuskybgd released under MIT License.
# See LICENSE file for full details.

from . import env
from .util import docformat


def run(args=[]):
    """
{b}NAME{o}
    {b}nuskybgd{o}
    Run nuskybgd tasks from command line.

{b}USAGE{o}
    nuskybgd task [arguments for task]

{b}DESCRIPTION{o}
    Without task arguments, the usage message for that task will be printed.
    The functions for these tasks can also be used in Python, either in an
    interactive session or in a script, by passing an arguments list to them
    in lieu of sys.argv.

    (For Python functions, see their __doc__ attribute.)

    The following tasks are available:

        {b}aspecthist   {o}Create aspect histogram image of pointing
        {b}mkinstrmap   {o}Make instrument map image
        {b}projbgd      {o}Create aperture background and detector mask images
                     rotated and convolved with the aspect histogram
        {b}fit          {o}Fit the background models to spectra from background
                     regions
        {b}spec         {o}Scale the fitted background model for a source
                     region
        {b}image        {o}Create images of the background sources
        {b}simplify     {o}Simplify the source + background Xspec save file by
                     removing the spectra from the background regions.
        {b}absrmf       {o}Add detector absorption to RMF files
    """
    tasks = {
        'aspecthist': aspecthist,
        'mkinstrmap': mkinstrmap,
        'projbgd': projbgd,
        'fit': fit,
        'image': image,
        'spec': spec,
        'simplify': simplify,
        'absrmf': absrmf
    }

    if len(args) == 1:
        print(docformat(run.__doc__))
        return 0

    if args[1] in tasks:
        return tasks[args[1]](args[1:])
    else:
        print('%s is unknown. Run \'nuskybgd\' to see all commands.' % args[1])


def absrmf(args=[]):
    """
{b}NAME{o}
    {b}nuskybgd absrmf{o}
    Create RMF files that includes detector absorption (DETABS).

{b}USAGE{o}
    absrmf evtfile outfile [rmffile=CALDB] [detabsfile=CALDB]

{b}DESCRIPTION{o}
    {b}evtfile{o} - An event file from which the INSTRUME and DATE-OBS
        keywords are taken.

    {b}outfile{o} - Will be prefixed to the output file names, and can be a
        file path.

    {b}rmffile{o} - The RMF file to multiply by absorption. Set it to CALDB
        (default) to use the latest CALDB file(s).

    {b}detabsfile{o} - Detector absorption file to multiply the RMF with. Set
        it to CALDB (default) to use the latest CALDB file(s).
    """
    if env.block() is True:
        return 1

    import os
    from . import rmf

    if len(args) not in (3, 4, 5):
        print(docformat(absrmf.__doc__))
        return 0

    evtfile = args[1]
    outfile = args[2]

    keywords = {
        'rmffile': 'CALDB',
        'detabsfile': 'CALDB'
    }

    for _ in args[3:]:
        arg = _.split('=')
        if arg[0] in keywords:
            keywords[arg[0]] = arg[1]

    # File exist checks
    halt = False
    if not os.path.exists(evtfile):
        print('%s not found.' % evtfile)
        halt = True
    if (keywords['rmffile'] != 'CALDB' and
            not os.path.exists(keywords['rmffile'])):
        print('%s not found.' % keywords['rmffile'])
        halt = True
    if (keywords['detabsfile'] != 'CALDB' and
            not os.path.exists(keywords['detabsfile'])):
        print('%s not found.' % keywords['detabsfile'])
        halt = True

    if halt:
        return 1

    # Overwrite output flag?
    overwrite = False
    if outfile[0] == '!':
        overwrite = True
        outfile = outfile[1:]

    rmf.make_absrmf(evtfile, outfile,
                    rmffile=keywords['rmffile'],
                    detabsfile=keywords['detabsfile'],
                    overwrite=overwrite)

    return 0


def fit(args=[]):
    """
{b}NAME{o}
    {b}nuskybgd fit{o}
    Fit NuSTAR background model

{b}USAGE{o}
    nuskybgd fit bgdinfo.json [savefile=bgdparams.xcm]

    nuskybgd fit --help  # print a sample bgdinfo.json

{b}DESCRIPTION{o}
    Generate a multi-component background model for spectra from several
    background regions and save the model containing preset normalizations to
    an xcm file. If the intended save file exists, will retry 99 times with a
    number (2 to 100) appended to the name.

    Required in bgdinfo.json:

        {b}bgfiles{o} - An array of spectra file names, extracted from
            background regions, grouped so that bins have gaussian statistics.

        {b}regfiles{o} - An array of region files for the background regions,
            in the same order as bgfiles.

        {b}refimgf{o} - Reference image file for creating region mask images,
            can be the background aperture image file.

        {b}bgdapfiles{o} - Dictionary with 'A' and 'B' keys, pointing to
            background aperture image files for each focal plane module.

        {b}bgddetfiles{o} - Dictionary with 'A' and 'B' keys, pointing to
            lists of 4 detector mask image files for each focal plane module.
    """
    if env.block() is True:
        return 1

    import os
    import json
    import xspec
    import numpy as np
    import astropy.io.fits as pf
    from . import util
    from . import model as numodel

    EXAMPLE_BGDINFO = """
Sample bgdinfo.json:

{
    "bgfiles": [
        "bgd1A_sr_g30.pha", "bgd1B_sr_g30.pha",
        "bgd2A_sr_g30.pha", "bgd2B_sr_g30.pha",
        "bgd3A_sr_g30.pha", "bgd3B_sr_g30.pha"
    ],

    "regfiles": [
        "bgd1A.reg", "bgd1B.reg",
        "bgd2A.reg", "bgd2B.reg",
        "bgd3A.reg", "bgd3B.reg"
    ],

    "refimgf": "bgdapA.fits",

    "bgdapfiles": {
        "A": "bgdapA.fits",
        "B": "bgdapB.fits"
    },

    "bgddetfiles": {
        "A": [
            "det0Aim.fits",
            "det1Aim.fits",
            "det2Aim.fits",
            "det3Aim.fits"
        ],
        "B": [
            "det0Bim.fits",
            "det1Bim.fits",
            "det2Bim.fits",
            "det3Bim.fits"
        ]
    },

    "fcxb_config": {
        "links": [
            [1, 2],
            [3, 4],
            [5, 6]
        ]
    },

    "intbgd_fix_line_ratios": true
}
"""

    if len(args) not in (2, 3):
        print(docformat(fit.__doc__))
        return 0

    if args[1] == '--help':
        print(docformat(fit.__doc__))
        print(EXAMPLE_BGDINFO)
        return 0

    keywords = {
        'savefile': None
    }

    for _ in args[2:]:
        arg = _.split('=')
        if arg[0] in keywords:
            keywords[arg[0]] = arg[1]

#    auxildir = env._AUX_DIR
    auxildir = env.auxildir

    # Input params
    bgdinfo = numodel.check_bgdinfofile(args[1])
    if bgdinfo is False:
        print('Error: background info file %s problem.' % args[1])
        return 1

    #####################################################
    presets = json.loads(open('%s/ratios.json' % auxildir).read())

    instrlist = numodel.get_keyword_specfiles(bgdinfo['bgfiles'],
                                              'INSTRUME', ext='SPECTRUM')

    # Compute aperture image and detector mask based weights using each
    # background region's mask
    regmask = util.mask_from_region(bgdinfo['regfiles'],
                                    bgdinfo['refimgf'])
    bgdapim, bgddetim = numodel.load_bgdimgs(bgdinfo)

    # Number of det pixels in the region mask.
    bgddetweights = numodel.calc_det_weights(bgddetim, regmask, instrlist)
    bgddetimsum = bgddetweights['sum']

    # Sum of the aperture image in the region mask.
    bgdapweights = numodel.calc_ap_weights(bgdapim, regmask, instrlist)
    bgdapimwt = bgdapweights['sum']

    refspec = numodel.get_refspec(instrlist)
    #####################################################

    # Interact with Xspec
    numodel.addspec(bgdinfo['bgfiles'])
    ##########
    # The model numbers below (they will become source number in Xspec) are
    # arbitrary but the same ones must be used in subsequent processes that
    # will load this Xspec save file.
    numodel.applymodel_apbgd(presets, refspec, bgdapimwt, 2)
    numodel.applymodel_intbgd(presets, refspec, bgddetimsum, 3,
        fix_line_ratios=bgdinfo['intbgd_fix_line_ratios'])
    numodel.applymodel_fcxb(refspec, bgddetimsum, 4)
    numodel.applymodel_intn(presets, refspec, bgddetimsum, 5)
    if 'fcxb_config' in bgdinfo and 'links' in bgdinfo['fcxb_config']:
        numodel.fcxb_linkab(bgdinfo['fcxb_config']['links'])
    ##########

    numodel.run_fit_settings()
    # numodel.run_fit()  # Disabled, to give user xcm file with prefilled
                         # values. Otherwise impossible to troubleshoot
                         # difficult background regions (e.g. containing
                         # additional sources)

    if keywords['savefile'] is None:
        numodel.save_xcm()
    else:
        numodel.save_xcm(prefix=keywords['savefile'].strip())

    return 0


def spec(args=[]):
    """
{b}NAME
    {b}nuskybgd spec{o}
    Rescale the background models for a source region

{b}USAGE{o}
    nuskybgd spec infofile.json bgdparams.xcm source.reg source.pha \\
        [savefile=bgd_source.xcm]

{b}DESCRIPTION{o}
    The result from {b}nuskybgd fit{o} is first retrieved by referring to
    infofile.json and loading bgdparams.xcm in Xspec. The source spectrum is
    then appended to the spectra file list and its corresponding parameters
    values for all the background model sources are rescaled. This state is
    saved to a new *.xcm file, which the user can load in Xspec to inspect the
    quality of the background model against the source spectrum.

    {b}infofile.json{o} - The same JSON file that was used with {b}nuskybgd
        fit{o} to obtain the background model.

    {b}bgdparams.xcm{o} - The Xspec save file from {b}nuskybgd fit{o}.

    {b}source.reg{o} - Region file containing the source region.

    {b}source.pha{o} - Source spectrum file, binned appropriately and
        containing header keywords specifying its RESPFILE.

    {b}savefile{o} - (Optional) Name of the output file. By default, 'bgd_' is
        prepended to the prefix of the first source region file, e.g.
        src1.reg,src2.reg -> bgd_src1.xcm.
    """
    if env.block() is True:
        return 1

    import os
    import json
    import xspec
    import numpy as np
    import astropy.io.fits as pf
    from . import util
    from . import model as numodel

    if len(args) not in (5, 6):
        print(docformat(spec.__doc__))
        return 0

    keywords = {
        'savefile': None
    }

    for _ in args[5:]:
        arg = _.split('=')
        if arg[0] in keywords:
            keywords[arg[0]] = arg[1]

    auxildir = env.auxildir

    # Input params
    bgdinfo = numodel.check_bgdinfofile(args[1])
    if bgdinfo is False:
        print('Error: background info file %s problem.' % args[1])
        return 1

    # Check the source region files and spec files
    src_regfiles = args[3].split(',')
    src_specfiles = args[4].split(',')
    halt = False
    for _ in src_regfiles + src_specfiles:
        if not os.path.exists(_):
            halt = True
            print('Error: file %s does not exist!' % _)
    if halt:
        return 1

    # Load the Xspec save file
    if not os.path.exists(args[2]):
        print('Error: save file %s does not exist!' % args[2])
        return 1
    xspec.AllData.clear()
    xspec.Xset.restore(args[2])

    # Check order of loaded spectra vs. files listed in bgdinfo. Must have
    # save order!
    if not numodel.check_spec_order(bgdinfo):
        print('Loaded spectra files are inconsistent with bgdinfo, stopping.')
        return 1

    #####################################################
    # This part is identical to fit()
    presets = json.loads(open('%s/ratios.json' % auxildir).read())

    instrlist = numodel.get_keyword_specfiles(bgdinfo['bgfiles'],
                                              'INSTRUME', ext='SPECTRUM')

    # Compute aperture image and detector mask based weights using each
    # background region's mask
    regmask = util.mask_from_region(bgdinfo['regfiles'],
                                    bgdinfo['refimgf'])
    bgdapim, bgddetim = numodel.load_bgdimgs(bgdinfo)

    # Number of det pixels in the region mask.
    bgddetweights = numodel.calc_det_weights(bgddetim, regmask, instrlist)
    bgddetimsum = bgddetweights['sum']

    # Sum of the aperture image in the region mask.
    bgdapweights = numodel.calc_ap_weights(bgdapim, regmask, instrlist)
    bgdapimwt = bgdapweights['sum']

    refspec = numodel.get_refspec(instrlist)
    #####################################################

    # Append the corresponding values for the source region
    src_instrlist = numodel.get_keyword_specfiles(src_specfiles, 'INSTRUME', ext='SPECTRUM')
    src_regmasks = util.mask_from_region(src_regfiles, bgdinfo['refimgf'])
    src_detweights = numodel.calc_det_weights(bgddetim, src_regmasks, src_instrlist)
    src_detimsum = src_detweights['sum']
    src_apweights = numodel.calc_ap_weights(bgdapim, src_regmasks, src_instrlist)
    src_apimwt = src_apweights['sum']

    instrlist += src_instrlist
    regmask += src_regmasks
    bgddetweights['fraction'] += src_detweights['fraction']
    # The two below are references. Don't sum the ['sum'] arrays a second time.
    bgddetimsum += src_detimsum
    bgdapimwt += src_apimwt

    src_number = len(bgdinfo['bgfiles']) + 1  # First source spectrum's number
    #####################################################


    # Interact with Xspec
    numodel.addspec(src_specfiles, fresh=False)
    ##########
    # These models need to have the same number as in fit() or it won't work!
    numodel.applymodel_apbgd(presets, refspec, bgdapimwt, 2, src_number=src_number)
    numodel.applymodel_intbgd(presets, refspec, bgddetimsum, 3, src_number=src_number,
                              fix_line_ratios=bgdinfo['intbgd_fix_line_ratios'])
    numodel.applymodel_fcxb(refspec, bgddetimsum, 4, src_number=src_number)
    numodel.applymodel_intn(presets, refspec, bgddetimsum, 5, src_number=src_number)
    ##########


    if keywords['savefile'] is None:
        numodel.save_xcm(prefix='bgd_'+src_regfiles[0].replace('.reg', ''))
    else:
        numodel.save_xcm(prefix=keywords['savefile'].strip())

    return 0


def image(args=[]):
    """
{b}NAME
    {b}nuskybgd image{o}
    Create images of the nuskybgd model sources.

{b}USAGE{o}
    nuskybgd image infofile.json bgdparams.xcm emin emax [prefix=]

{b}DESCRIPTION{o}
    The result from {b}nuskybgd fit{o} is first retrieved by referring to
    infofile.json and loading bgdparams.xcm in Xspec. This is used to create
    image models of each background source defined in auxil/ratios.json are
    created, by getting model predicted counts in the background regions from
    Xspec and extrapolating to the whole field of view.

    {b}infofile.json{o} - The same JSON file that was used with {b}nuskybgd
        fit{o} to obtain the background model.

    {b}bgdparams.xcm{o} - The Xspec save file from {b}nuskybgd fit{o}.

    {b}emin{o} - Lower energy bound in keV.

    {b}emax{o} - Upper energy bound in keV.

    {b}prefix{o} - (Optional) Prefix for output files.

    The following files will be created (or overwritten):
    bgd_apbgd_[A/B].fits, bgd_intbgd_[A/B].fits, bgd_intn_[A/B].fits,
    bgd_fcxb_[A/B].fits, as well as bgd_fcxb_[spectrum_number].fits for
    each background spectrum.
    """
    import os
    import json
    import xspec
    import numpy as np
    import astropy.io.fits as pf
    from . import util
    from . import model as numodel

    if env.block() is True:
        return 1

    if len(args) not in (5, 6):
        print(docformat(image.__doc__))
        return 0

    keywords = {
        'prefix': ''
    }

    for _ in args[2:]:
        arg = _.split('=')
        if arg[0] in keywords:
            keywords[arg[0]] = arg[1]

    auxildir = env.auxildir

    # Input params
    bgdinfo = numodel.check_bgdinfofile(args[1])
    if bgdinfo is False:
        print('Error: background info file %s problem.' % args[1])
        return 1

    # Load the Xspec save file
    if not os.path.exists(args[2]):
        print('Error: save file %s does not exist!' % args[2])
        return 1
    xspec.AllData.clear()
    xspec.Xset.restore(args[2])

    # Check order of loaded spectra vs. files listed in bgdinfo. Must have
    # save order!
    if not numodel.check_spec_order(bgdinfo):
        print('Loaded spectra files are inconsistent with bgdinfo, stopping.')
        return 1

    #####################################################
    # This part is identical to fit()
    presets = json.loads(open('%s/ratios.json' % auxildir).read())

    instrlist = numodel.get_keyword_specfiles(bgdinfo['bgfiles'],
                                              'INSTRUME', ext='SPECTRUM')

    # Compute aperture image and detector mask based weights using each
    # background region's mask
    regmask = util.mask_from_region(bgdinfo['regfiles'],
                                    bgdinfo['refimgf'])
    bgdapim, bgddetim = numodel.load_bgdimgs(bgdinfo)

    # Number of det pixels in the region mask.
    bgddetweights = numodel.calc_det_weights(bgddetim, regmask, instrlist)
    bgddetimsum = bgddetweights['sum']

    # Sum of the aperture image in the region mask.
    bgdapweights = numodel.calc_ap_weights(bgdapim, regmask, instrlist)
    bgdapimwt = bgdapweights['sum']

    refspec = numodel.get_refspec(instrlist)
    
    #####################################################

    emin = float(args[3])
    emax = float(args[4])

    ignore_string = '**-%.2f %.2f-**' % (emin, emax)

    model_norms = numodel.read_xspec_model_norms()
    numodel.zero_xspec_model_norms(model_norms)


    numodel.bgimg_apbgd(bgdinfo, model_norms, bgdapweights, refspec,
        ignore=ignore_string, outprefix=keywords['prefix'])

    numodel.bgimg_intbgd(presets, refspec, bgdinfo, model_norms, bgddetweights, bgddetim,
        ignore=ignore_string, outprefix=keywords['prefix'])

    numodel.bgimg_fcxb(bgdinfo, model_norms, bgddetweights, bgddetim, regmask,
        ignore=ignore_string, outprefix=keywords['prefix'])

    numodel.bgimg_intn(presets, refspec, bgdinfo, model_norms, bgddetweights, bgddetim,
        ignore=ignore_string, outprefix=keywords['prefix'])

    return 0


def simplify(args=[]):
    """
{b}NAME
    {b}nuskybgd simplify{o}
    Remove the background region spectra from the Xspec save file

{b}USAGE{o}
    nuskybgd simplify infofile.json bgd_src.xcm [savefile=bgd_only_src.xcm]

{b}DESCRIPTION{o}
    The result from {b}nuskybgd spec{o} is loaded in Xspec, and the first few
    spectra in the list are removed (based on how many entries there are in
    bgdinfo). Only source spectra and their background models remain. All
    parameters are frozen, making it ready for the user to load the resulting
    save file and add their source model.

    {b}infofile.json{o} - The same JSON file that was used with {b}nuskybgd
        fit{o} to obtain the background model.

    {b}bgd_src.xcm{o} - The Xspec save file from {b}nuskybgd spec{o}.

    {b}savefile{o} - (Optional) Name of the output file. By default, if 'bgd_'
        is in the input xcm file name, it is replaced by 'bgd_only_' for the
        output file name. Otherwise, 'bgd_only_' is prepended to the input
        file name.
    """

    if env.block() is True:
        return 1

    import xspec
    import json
    import os
    from . import model as numodel

    if len(args) not in (3, 4):
        print(docformat(simplify.__doc__))
        return 0

    keywords = {
        'savefile': None
    }

    for _ in args[3:]:
        arg = _.split('=')
        if arg[0] in keywords:
            keywords[arg[0]] = arg[1]

    # Input params
    bgdinfo = numodel.check_bgdinfofile(args[1])
    if bgdinfo is False:
        print('Error: background info file %s problem.' % args[1])
        return 1

    # Load the Xspec save file
    if not os.path.exists(args[2]):
        print('Error: save file %s does not exist!' % args[2])
        return 1
    xspec.AllData.clear()
    xspec.Xset.restore(args[2])

    nbgd = len(bgdinfo['bgfiles'])
    nspec = xspec.AllData.nSpectra
    numodel.remove_ispec(nbgd)
    frozenpars = numodel.freeze_pars(list(range(1, nspec - nbgd + 1)))

    if keywords['savefile'] is None:
        basename = os.path.basename(args[2]).replace('.xcm', '')
        if 'bgd_' in args[2]:
            numodel.save_xcm(prefix=basename.replace('bgd_', 'bgd_only_'))
        else:
            numodel.save_xcm(prefix='bgd_only_'+basename)
    else:
        numodel.save_xcm(prefix=keywords['savefile'].strip())

    return 0


def mkinstrmap(args=[]):
    """
{b}NAME
    {b}nuskybgd mkinstrmap{o}
    Create an instrument map for a specific observation

{b}USAGE{o}
    nuskybgd mkinstrmap A01_cl.evt [usrbpix=usrbpix.fits] [prefix=prefix]

{b}DESCRIPTION{o}
    Bad pixels from CALDB will always be applied.

    {b}usrbpix{o} - Specify additional bad pixels using a FITS file that has
        BADPIX extension(s). Multiple files can be given, separated by commas
        (with no space in the argument).

    {b}prefix{o} - Add a file name prefix for the output file. This allows the
        output to be written to any specified path. If a directory is
        intended, it must end with '/'.

    {b}dryrun{o} - If 'yes', the task will stop after the bad pixel files have
        been checked, before the instrument map is calculated.
    """
    import os
    from .image import (
        apply_badpix, get_badpix_exts, make_det_mask, get_caldb_instrmap
    )
    from .util import fpm_parse, print_hr
    import astropy.io.fits as pf

    if env.block() is True:
        return 1

    ERR_BPIXFILE_NOTFOUND = """
Error: bad pixel file not found, skipping:
{file}
"""

    if len(args) == 1:
        print(docformat(mkinstrmap.__doc__))
        return 0

    obsinfofile = args[1]

    keywords = {
        'usrbpix': None,
        'prefix': '',
        'dryrun': ''
    }

    for _ in args[1:]:
        arg = _.split('=')
        if arg[0] in keywords:
            keywords[arg[0]] = arg[1]

    halt = False

    # Inputs File with observation info must exist. Typically this is the
    # event file itself, but this program needs only the information in its
    # FITS header.
    require_hdr_keywords = ['INSTRUME']
    if not os.path.exists(obsinfofile):
        print('Error: %s does not exist.' % obsinfofile)
        halt = True
        return 1  # No further checks if the main file doesn't exist

    evtfile = pf.open(obsinfofile)
    if 'EVENTS' not in evtfile:
        print('Error: EVENTS extension not found in %s.' % obsinfofile)
        halt = True
        return 1

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
        return 1

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

    if halt:
        print('mkinstrmap did not complete. See error messages.')
        return 1

    print_hr()
    print('Running nuskybgd mkinstrmap:')
    print(' '.join(args))
    print_hr()
    print('Creating instrument maps for FPM%s for observation %s' % (
        ab, obsinfofile))

    caldbbpixpath = env._CALDB.getBADPIX(
        evthdr['INSTRUME'], 'DET0', evthdr['DATE-OBS'])
    bpixfiles.append('%s/%s' % (env._CALDB_PATH, caldbbpixpath))

    print('Collecting bad pixel lists...')
    bpixexts = get_badpix_exts(bpixfiles)

    # Stop here if doing a dry run
    if keywords['dryrun'] == 'yes':
        return 0

    print('Loading CALDB instrmap...')
    instrmap, header = get_caldb_instrmap(evthdr)
    print('Loading CALDB pixpos...')
    pixmap, detnum = make_det_mask(evthdr)

    print('Applying bad pixels lists for FPM%s...' % ab)
    masked_instrmap = apply_badpix(instrmap, bpixexts[ab], pixmap, detnum)

    pf.HDUList(
        pf.PrimaryHDU(masked_instrmap, header=header)
    ).writeto(outfilename)
    print('Saved to %s' % outfilename)

    print('Done.')
    print_hr()

    return 0


def projbgd(args=[]):
    """
{b}NAME{o}
    {b}nuskybgd projbgd{o}
    Project the instrument map onto sky coordinates

{b}USAGE{o}
    nuskybgd projbgd refimg=flux.fits out=output.fits \\
        mod=[A,B] det=[1,2,3,4] chipmap=chipmap.fits aspect=aspect.fits

    Example:

    nuskybgd projbgd refimg=imA4to25keV.fits out=bgdapA.fits \\
        mod=A det=1234 chipmap=newinstrmapA.fits aspect=aspecthistA.fits

{b}DESCRIPTION{o}
    The output file bgdapA.fits has the background aperture image in sky
    coordinates. A series of files showing the projected detector masks for
    each detector specified by det= will also be created, named
    det[0-3]Aim.fits.
    """
    if env.block() is True:
        return 1

    if len(args) == 1:
        print(docformat(projbgd.__doc__))
        return 0

    import os
    import sys
    import numpy as np
    import astropy.io.fits as pf
    from .image import (
        get_det_mask, get_aspect_hist_image, get_aspect_hist_peak,
        get_aperture_image, transform_image, calc_pa_to_arot
    )
    opts = {}
    keys = ('refimg', 'out', 'det', 'mod', 'chipmap', 'aspect')
    # TODO: validate input chipmap and aspect files
    for _ in args[1:]:
        pars = _.split('=')
        if pars[0] in keys:
            opts[pars[0]] = pars[1]

    for _ in keys:
        if _ not in opts:
            print('Missing %s= keyword' % _)
            return 1

    try:
        refimg = pf.open(opts['refimg'])
    except FileNotFoundError:
        print('refimg file not found')
        return 1

    try:
        hdr = refimg[0].header
        pa = hdr['PA_PNT'] + 1.0
        required = ('CTYPE1', 'CTYPE2', 'CRPIX1', 'CRPIX2',
                    'CRVAL1', 'CRVAL2', 'CDELT1', 'CDELT2')
        for _ in required:
            if _ not in hdr:
                raise KeyError('WCS info missing from header' % _)
        if hdr['NAXIS1'] != 1000 or hdr['NAXIS2'] != 1000:
            raise ValueError('Unexpected image size (not 1000x1000)')
    except KeyError:
        print('Required information not found in header')
        return 1
    except ValueError:
        print('Unexpected refimg format')
        return 1

    if opts['mod'] not in ('A', 'B'):
        print('mod= must be set to A or B')
        return 1

    det_input = []
    for _ in '1234':
        if _ in opts['det']:
            det_input.append(int(_))
    if len(det_input) == 0:
        print('No valid detector number specified.')
        return 1

    auxildir = env.auxildir

    if opts['out'][0] != '!':
        halt = False
        if os.path.exists(opts['out']):
            print('Output file %s exists.' % opts['out'])
            halt = True
        for idet in np.arange(4):
            _ = 'det%d%sim.fits' % (idet, opts['mod'])
            if os.path.exists(_):
                print('Output file %s exists.' % _)
                halt = True
        if halt:
            return 1

    print('Using auxiliary data in %s' % auxildir)

    detmask = get_det_mask(opts['chipmap'], det_input)
    apim = get_aperture_image(opts['mod'])
    posim, pos_xoff, pos_yoff = get_aspect_hist_image(opts['aspect'])
    totalexp = np.sum(posim)

    arot = calc_pa_to_arot(pa)

    asppeakx, asppeaky = get_aspect_hist_peak(posim, pos_xoff, pos_yoff)

    xgrid, ygrid = np.meshgrid(np.arange(1000),
                               np.arange(1000), indexing='xy')

    detxa = np.int32(350 + np.around(- np.cos(arot) * (xgrid - asppeakx) +
                                     np.sin(arot) * (ygrid - asppeaky)))
    detya = np.int32(350 + np.around(np.sin(arot) * (xgrid - asppeakx) +
                                     np.cos(arot) * (ygrid - asppeaky)))
    nudge = [2.7, 0.8]
    nudge = [nudge[0] * np.cos(-arot) - nudge[1] * np.sin(-arot),
             nudge[0] * np.sin(-arot) + nudge[1] * np.cos(-arot)]

    imval = detmask / totalexp  # Normalize it, since convolving with posim
    apval = apim / totalexp        # multiplies it by sum(posim)

    refimap = np.zeros((1000, 1000), dtype=np.float64)
    refimi = np.zeros((len(det_input), 1000, 1000), dtype=np.float64)

    print('Rotating instrument map (PA = %f, %f radians)' % (pa, arot))
    for i in np.arange(1000):
        for j in np.arange(1000):
            detx = detxa[j, i]
            dety = detya[j, i]
            for idet_input in range(len(det_input)):
                if ((0 <= detx < 360) and (0 <= dety < 360) and
                        (detmask[idet_input, dety, detx] > 0)):
                    refimap[j, i] = apval[dety, detx]
                    refimi[idet_input, j, i] = imval[idet_input, dety, detx]

    bgdimap = np.zeros((1000, 1000), dtype=np.float64)
    bgdimi = np.zeros((len(det_input), 1000, 1000), dtype=np.float64)

    inx = np.where(posim > 0)
    progress_total = len(inx[0])
    progress = 0
    print('Convolving instrument map with aspect solution...')
    x_offset = pos_xoff - asppeakx - nudge[0]
    y_offset = pos_yoff - asppeaky - nudge[1]
    for x, y in zip(inx[1], inx[0]):
        bgdimap += posim[y, x] * transform_image(
            refimap, x + x_offset, y + y_offset, 0)
        for idet_input in range(len(det_input)):
            bgdimi[idet_input] += posim[y, x] * transform_image(
                refimi[idet_input], x + x_offset, y + y_offset, 0)
        progress += 1
        sys.stdout.write('\r %d/%d positions' % (progress, progress_total))
        sys.stdout.flush()

    print('\nDone.')

    outfile = pf.HDUList([pf.PrimaryHDU(bgdimap, header=hdr)])

    # Write background aperture image
    outfile.writeto(opts['out'], overwrite=True)

    # Write detector mask images
    for idet_input in range(len(det_input)):
        pf.HDUList([pf.PrimaryHDU(bgdimi[idet_input], header=hdr)]).writeto(
            'det%d%sim.fits' % (det_input[idet_input] - 1, opts['mod']),
            overwrite=True)  # For det number, input 1-4 -> output 0-3

    return 0


def aspecthist(args=[]):
    """
{b}NAME{o}
    {b}nuskybgd aspecthist{o}
    Create an aspect histogram image from pointing info after filtering by GTI

{b}USAGE{o}
    nuskybgd aspecthist aimpoints.fits gtifile=gti.fits out=aspecthist.fits

    Example:

    nuskybgd aspecthist nu90201039002A_det1.fits out=aspecthistA.fits \\
        gtifile=nu90201039002A01_gti.fits

    nuskybgd aspecthist nu90201039002B_det1.fits out=aspecthistB.fits \\
        gtifile=nu90201039002B01_gti.fits

{b}DESCRIPTION{o}
    The output file has an image representing the 2D histogram of pointing
    position with time. Gets pointing position from nu%obsid%?_det?.fits and
    good time intervals from nu%obsid%?0?_gti.fits. nu%obsid%?_det?.fits (e.g.
    nu90201039002A_det1.fits) has a table block named 'DET?_REFPOINT', with 3
    columns (TIME, X_DET?, Y_DET?) where ? is the detector number (1-4).
    nu%obsid%?0?_gti.fits (e.g. nu90201039002A01_gti.fits) has a table block
    named 'STDGTI' with two columns (START, STOP) listing intervals of good
    times.

    The image is in an extension with the name ASPECT_HISTOGRAM. Any zero
    padding has been cropped, and the x and y offsets are recorded in the
    header keywords X_OFF and Y_OFF.
    """
    import os
    import astropy.io.fits as pf
    from .util import (
        check_header_aimpoint, check_header_gti, filter_gti
    )
    from .image import (
        make_aspecthist_img, write_aspecthist_img
    )

    if len(args) != 4:
        print(docformat(aspecthist.__doc__))
        return 0

    aimpointfile = args[1]

    keywords = {
        'gtifile': None,
        'out': None
    }

    for _ in args[1:]:
        arg = _.split('=')
        if arg[0] in keywords:
            keywords[arg[0]] = arg[1]

    gtifile = keywords['gtifile']
    outfilename = keywords['out']

    # Check the arguments
    scriptname = os.path.basename(__file__)
    print('Running %s. Performing checks...' % scriptname)

    halt = False

    if not os.path.exists(aimpointfile):
        print('aimpointfile not found: %s' % aimpointfile)
        halt = True

    if gtifile is None:
        print('gtifile= missing')
        halt = True
    elif not os.path.exists(gtifile):
        print('GTI file not found: %s' % gtifile)
        halt = True

    if outfilename is None:
        print('out= missing')
        halt = True
    elif os.path.exists(outfilename) and keywords['out'][0] != '!':
        print('Output file exists: %s' % outfilename)
        halt = True

    if halt:
        print('%s did not complete. See error messages.' % scriptname)
        return 1

    # Check contents of aimpoint and gti files
    aimpointext = None
    detfh = pf.open(aimpointfile)
    for ext in detfh:
        if check_header_aimpoint(ext.header):
            aimpointext = ext
            break
    if aimpointext is None:
        print('No aspect info in the specified file %s.' % aimpointfile)
    else:
        print('Found extension %s.' % aimpointext.header['EXTNAME'])

    gtiext = None
    gtifh = pf.open(gtifile)
    for ext in gtifh:
        if check_header_gti(ext.header):
            gtiext = ext
            break
    if gtiext is None:
        print('Error: no GTI info in the specified file %s.' % gtifile)
        return 1
    else:
        print('Found extension %s.' % gtiext.header['EXTNAME'])

    # Check output file
    if os.path.exists(outfilename):
        if outfilename[0] == '!':
            outfilename = outfilename[1:]
            print('Overwriting file %s...')
        else:
            print('Output file %s exists (prefix with ! to overwrite)')
            return 1

    coords, dt = filter_gti(aimpointext, gtiext)

    asphistimg, x_min, x_max, y_min, y_max = make_aspecthist_img(coords, dt)

    write_aspecthist_img(
        asphistimg, outfilename, aimpointext,
        (x_min, x_max, y_min, y_max), overwrite=True)

    return 0


def main_nuskybgd():
    """
    Entry point for main CLI interface
    """
    import sys
    code = run(sys.argv)
