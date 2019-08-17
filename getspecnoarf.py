#!/usr/bin/env python3
"""
Example usage:

getspecnoarf.py nu90201039002A01_cl.evt reg=src.reg \\
    indir=. outdir=. stemout=src \\
    attfile=../auxil/nu90201039002_att.fits.gz
"""

import sys
import os
import astropy.io.fits as pf

if len(sys.argv) == 1:
    print(__doc__)
    sys.exit(0)

keywords = {
    'indir': '.',
    'outdir': '.',
    'outprefix': None,
    'reg': None,
    'attfile': None
}

for _ in sys.argv[2:]:
    arg = _.split('=')
    if arg[0] in keywords:
        keywords[arg[0]] = arg[1]

halt = False

if not os.path.exists(keywords['attfile']):
    print('Error: attfile %s not found.' % keywords['attfile'])
    halt = True

if not os.path.exists(keywords['reg']):
    print('Error: reg file %s not found.' % keywords['reg'])
    halt = True

if keywords['outprefix'] is None:
    keywords['outprefix'] = keywords['reg'].replace('.reg', '')

# Check CALDB and CALDBCONFIG envs needed by the nustar pipeline tools
if 'CALDB' not in os.environ:
    print('You must set the $CALDB environment variable.')
    print('Point to the folder with data, docs, and software subdirectories')
    halt = True
else:
    caldbnustar = '%s/data/nustar/fpm/bcf' % os.environ['CALDB']
    if not os.path.exists(caldbnustar):
        print('NuSTAR CALDB not found at %s!' % caldbnustar)
        halt = True

if 'CALDBCONFIG' not in os.environ:
    print('You must set the $CALDBCONFIG environment variable.')
    print('This is usually at $CALDB/software/tools/caldb.config')
    halt = True
else:
    if not os.path.exists(os.environ['CALDBCONFIG']):
        print('CALDBCONFIG not found at %s!' % os.environ['CALDBCONFIG'])
        halt = True

if halt:
    sys.exit(1)

evtfile = sys.argv[1]
evtfh = pf.open(evtfile)
evthdr = evtfh['EVENTS'].header

obsid = evthdr['OBS_ID']
stemin = 'nu%s' % obsid
det = evthdr['INSTRUME'][-1]

indir = keywords['indir']
if len(indir) > 0 and indir[-1] != '/':
    indir = '%s/' % indir
stemout = keywords['outprefix']
outdir = keywords['outdir']
if len(outdir) > 0 and outdir[-1] != '/':
    outdir = '%s/' % outdir
attfile = keywords['attfile']
regfile = keywords['reg']

makearf = 'no'
makermf = 'yes'

nuproductscmd = """nuproducts \
infile={evtfile} srcregionfile={regfile} indir={indir} outdir={outdir} \
runmkarf={makearf} runmkrmf={makermf} bkgextract=no lcfile=NONE \
instrument=FPM{det} steminputs={stemin} stemout={stemout} boxsize=20 \
attfile={attfile} clobber=yes"""
os.system(nuproductscmd.format(
    evtfile=evtfile,
    regfile=regfile,
    indir=indir,
    outdir=outdir,
    makearf=makearf,
    makermf=makermf,
    det=det,
    stemin=stemin,
    stemout=stemout,
    attfile=attfile
))

# os.system("mv " + dir + "event_cl/" + outdir + "/" + regcore + "_sr.rmf " +
#           dir + "event_cl/" + outdir + "/" + regcore + "_sr_orig.rmf")
# os.system("cmprmf " + dir + "event_cl/" + outdir + "/" + regcore +
#           "_sr_orig.rmf " +
#           dir + "event_cl/" + outdir + "/" + regcore + "_sr.rmf 1e-6")
# os.system("rm -f " + dir + "event_cl/" +
#           outdir + "/" + regcore + "_sr_g30.pha")
os.system("""grppha {outdir}{stemout}_sr.pha \
\\!{outdir}{stemout}_sr_g30.pha 'group min 30 & exit'""".format(
    indir=indir, outdir=outdir, stemout=stemout))
