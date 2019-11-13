#!/usr/bin/env python3


if __name__ == '__main__':
    """
    Create an instrument map for a specific observation. The user can specify
    additional bad pixels by giving files that have BADPIX extension(s) using
    the usrbpix= argument. Multiple files can be given, separated by commas
    (with no space in the argument). Bad pixels from CALDB will always be
    applied.

    Usage:

    nuskybgd instrmap A01_cl.evt [usrbpix=usrbpix.fits] [prefix=prefix]

    The prefix= argument lets user add a file name prefix for the output file.
    This allows the output to be written to any specified path. If a directory
    is intended, it must end with '/'.

    If dryrun=yes the task will stop after the bad pixel files have been
    checked, before the instrument map is calculated.
    """
    from .image import (
        apply_badpix, get_badpix_exts, make_det_mask, get_caldb_instrmap
    )
    from .util import fpm_parse

    ERR_CALDB_NOTFOUND = """
Error: checking for calibration files in $CALDB/data/nustar, but cannot find
caldb.indx in this CALDB location. Is the $CALDB path correct?

$CALDB = {caldbpath}

Cannot find the file

{caldbinxpath}
"""

    ERR_BPIXFILE_NOTFOUND = """
Error: bad pixel file not found, skipping:
{file}
"""

    if len(sys.argv) == 1:
        print(__doc__)
        sys.exit(0)

    obsinfofile = sys.argv[1]

    keywords = {
        'usrbpix': None,
        'prefix': '',
        'dryrun': ''
    }

    for _ in sys.argv[1:]:
        arg = _.split('=')
        if arg[0] in keywords:
            keywords[arg[0]] = arg[1]

    scriptname = os.path.basename(__file__)
    print('Running %s. Performing checks...' % scriptname)

    halt = False

    # Inputs File with observation info must exist. Typically this is the
    # event file itself, but this program needs only the information in its
    # FITS header.
    require_hdr_keywords = ['INSTRUME']
    if not os.path.exists(obsinfofile):
        print('Error: %s does not exist.' % obsinfofile)
        halt = True
        sys.exit(1)  # No further checks if the main file doesn't exist

    evtfile = pf.open(obsinfofile)
    if 'EVENTS' not in evtfile:
        print('Error: EVENTS extension not found in %s.' % obsinfofile)
        halt = True
        sys.exit(1)

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
        sys.exit(1)

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
        print('instrmap did not complete. See error messages.')
        sys.exit(1)

    print('Preliminary checks complete, proceeding...')

    print('Creating instrument maps for FPM%s for observation %s' % (
        ab, obsinfofile))

    caldb = nuskybgd.caldb.CalDB(os.environ['CALDB'])
    caldbbpixpath = caldb.getBADPIX(
        evthdr['INSTRUME'], 'DET0', evthdr['DATE-OBS'])
    bpixfiles.append('%s/%s' % (os.environ['CALDB'], caldbbpixpath))

    print('Collecting bad pixel lists...')
    bpixexts = get_badpix_exts(bpixfiles)

    # Stop here if doing a dry run
    if keywords['dryrun'] == 'yes':
        sys.exit(0)

    print('Loading CALDB instrmap...')
    instrmap, header = get_caldb_instrmap(evthdr)
    print('Loading CALDB pixpos...')
    pixmap, detnum = make_det_mask(evthdr)

    print('Applying bad pixels lists for FPM%s...' % ab)
    masked_instrmap = apply_badpix(instrmap, bpixexts[ab], pixmap, detnum)

    pf.HDUList(
        pf.PrimaryHDU(masked_instrmap, header=header)
    ).writeto(outfilename)
