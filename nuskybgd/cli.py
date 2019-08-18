"""
Functions to handle using nuskybgd tasks from the command line. These
functions can also be used in Python, either in an interactive session or in a
script, by passing an arguments list to them in lieu of sys.argv.
"""

def run(args=[]):
    """
    Run nuskybgd tasks.

    Usage: nuskybgd task [arguments for task]

    Without task arguments, the usage message for that task will be printed.
    The following tasks are available:
    """
    tasks = {
        'absrmf': absrmf
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
