"""
Functions to handle using nuskybgd tasks from the command line. These
functions can also be used in Python, either in an interactive session or in a
script, by passing an arguments list to them in lieu of sys.argv.
"""
from . import conf

def run(args=[]):
    """
    Run nuskybgd tasks.

    Usage: nuskybgd task [arguments for task]

    Without task arguments, the usage message for that task will be printed.
    The following tasks are available:
    """
    tasks = {
        'absrmf': absrmf,
        'fitab': fitab
    }

    if len(args) == 1:
        print(run.__doc__)
        print('\n'.join(list('        %s' % _ for _ in sorted(tasks.keys()))))
        return 0

    if args[1] in tasks:
        return tasks[args[1]](args[1:])


def absrmf(args=[]):
    """
    Create RMF files that includes detector absorption (DETABS).

    Usage:

    absrmf evtfile outfile [rmffile=CALDB] [detabsfile=CALDB]

    evtfile is an event file from which the INSTRUME and DATE-OBS keywords are
    used.

    outfile will be prefixed to the output file names, and can be a file path.

    rmffile is the RMF file to multiply by absorption, set it to CALDB
    (default) to use the latest CALDB file(s).

    detabsfile is the detector absorption file to multiply the RMF with, set
    it to CALDB (default) to use the latest CALDB file(s).
    """
    if conf.block() is True:
        return 1

    import os
    from . import rmf

    if len(args) not in (3, 4, 5):
        print(absrmf.__doc__)
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


def fitab(args=[]):
    """
    Fit a multi-component background model to spectra from several background
    regions and save the best-fit model to an xcm file. If the intended save
    file exists, will retry 99 times with a number (2 to 100) appended to the
    name.

    Usage:

    fitab bgdinfo.json [savefile=bgdparams.xcm]

    fitab --help  # print a sample bgdinfo.json

    Required in bgdinfo.json:

        bgfiles - An array of spectra file names, extracted from background
            regions, grouped so that bins have gaussian statistics.

        regfiles - An array of region files for the background regions, in the
            same order as bgfiles.

        refimgf - Reference image file for creating region mask images, can be
            the background aperture image file.

        bgdapfiles - Dictionary with 'A' and 'B' keys, pointing to background
            aperture image files for each focal plane module.

        bgddetfiles - Dictionary with 'A' and 'B' keys, pointing to lists of 4
            detector mask image files for each focal plane module.
    """
    if conf.block() is True:
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
    }
}
"""

    if len(args) not in (2, 3):
        print(fitab.__doc__)
        return 0

    if args[1] == '--help':
        print(fitab.__doc__)
        print(EXAMPLE_BGDINFO)
        return 0

    keywords = {
        'savefile': None
    }

    for _ in args[2:]:
        arg = _.split('=')
        if arg[0] in keywords:
            keywords[arg[0]] = arg[1]


    auxildir = conf._AUX_DIR

    # Input params
    bgdinfofile = args[1]

    bgdinfo = numodel.check_bgdinfofile(bgdinfofile)

    if bgdinfo is False:
        print('Error: background info file %s problem.' % bgdinfofile)
        return 1

    ratios = json.loads(open('%s/ratios.json' % auxildir).read())
    # In [45]: ratiosA.keys()
    # Out[45]: dict_keys(['name', 'comment', 'models'])

    # bgfiles and regfiles must have the same ordering
    bgfiles = bgdinfo['bgfiles']
    regfiles = bgdinfo['regfiles']
    refimgf = bgdinfo['refimgf']
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

    spectra = []

    xspec.DataManager.clear(0)  # Clear any existing loaded data

    # Load each spectrum as a new data group
    for i in range(len(bgfiles)):
        spectra.append(xspec.AllData('{num}:{num} {file}'.format(
            num=i + 1,
            file=bgfiles[i])))

    # Check for valid INSTRUME keyword in spectrum header: need this to
    # determine which FPM is used.

    # Reference spectrum's index for the linked model parameters
    # 0-based; note that xspec.AllData() is 1-based
    refspec = {'A': None, 'B': None}

    for i in range(xspec.AllData.nSpectra):
        spec = xspec.AllData(i + 1)

        try:
            fpm = util.fpm_parse(spec.fileinfo('INSTRUME'))
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
            auxildir, util.fpm_parse(spec.fileinfo('INSTRUME')))  # 4:fcxb
        spec.multiresponse[5].arf = '%s/be.arf' % auxildir  # 6:grxe

    # Compute aperture image and detector mask based weights using each
    # background region's mask

    # Each background spectrum has a list of 4 values associated with each
    # CCD: number of pixels in the region mask.
    bgddetimsum = []

    # Each background spectrum has a single value that is the sum of the
    # aperture image in the region mask.
    bgdapimwt = []

    for i in range(xspec.AllData.nSpectra):
        spec = xspec.AllData(i + 1)
        fpm = util.fpm_parse(spec.fileinfo('INSTRUME'))
        regmask = util.mask_from_region(regfiles[i], refimgf)
        detnpix = [np.sum(regmask * detim) for detim in bgddetim[fpm]]
        bgddetimsum.append(detnpix)
        bgdapimwt.append(np.sum(regmask * bgdapim[fpm]))

    numodel.addmodel_apbgd(ratios, refspec, bgdapimwt, 2)
    numodel.addmodel_intbgd(ratios, refspec, bgddetimsum, 3)
    numodel.addmodel_fcxb(refspec, bgddetimsum, 4)
    numodel.addmodel_intn(ratios, refspec, bgddetimsum, 5)

    numodel.run_fit()

    if keywords['savefile'] is None:
        numodel.save_xcm()
    else:
        numodel.save_xcm(prefix=keywords['savefile'].strip())

    return 0
