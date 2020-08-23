# Requirements

- Nuskybgd interacts with Xspec via PyXspec; you must install HEASOFT from source to use PyXspec.

    Download HEASOFT: https://heasarc.gsfc.nasa.gov/lheasoft/download.html

- You need a Python 3 (3.6 or newer) environment; compile PyXspec against it. See SETUP.


# Quick start guide

## 1. Environment variables

These environment variables must be set up before anything else. Modify the
paths to point to the correct location on your machine.

```bash
# bash, zsh

# NUSKYBGD
export NUSKYBGD=$HOME/astro/soft/nuskybgd
export NUSKYBGD_AUXIL=$NUSKYBGD/auxil

# HEASARC CALDB
export CALDB=/soft/astro/heasarc/CALDB
export CALDBCONFIG=$CALDB/software/tools/caldb.config

# Add nuskybgd executables to $PATH
export PATH=$NUSKYBGD/bin:$PATH

# Add nuskybgd to $PYTHONPATH for its modules to be found
export PYTHONPATH=$NUSKYBGD:$PYTHONPATH

# In csh and tcsh the export statement is setenv instead (use double quotes)
# setenv PATH "$NUSKYBGD/bin:$PATH"
```

In order to not clutter your environment variables and paths, you can wrap the
above in a file and only run the export statements when you need them. For
example, you can save the above to a file `~/nuskybgd-init.sh`. Then in your
shell rc file, add an alias:

```bash
# bash, zsh
alias nuskybgdinit="source ~/nuskybgd-init.sh"

# csh, tcsh
# alias nuskybgdinit "source ~/nuskybgd-init.sh"
```

The next time you start a shell session, the new alias will be loaded, and you
can set the environment variables with `nuskybgdinit`.


Assuming you have the 'default' layout for files, e.g. cleaned event files are
in `[target name]/[obs id]/event_cl`, the examples below puts background
modelling related files in a folder named `bgd` inside `event_cl`.

Create the `bgd/` folder if this is the first time.

---

## 2. Preparations

### 2.1 Make an image

Use mkimgs.py to create a counts image for WCS reference.

```bash
# (cd into the top level)
./mkimgs.py ./ 50002031004 3 20
```

This creates the file `imA3to20keV.fits` and `imB3to20keV.fits`, which are
counts images between 3-20 keV seen by each detector module. These images are
used later for their WCS information. All subsequent image products are
projections onto the same WCS grid as these.



### 2.2 Extract spectra from background regions

> **Note about regions**
>
> Region masking is handled by pyregion. If in doubt, test whether your region
> results in the mask as expected, directly with pyregion.
>
> In general:
>
> The mask is created by rendering each entry in the region file in sequence,
> changing pixels to `1` for an include region or changing pixels to `0` for
> an exclude region. Therefore, order matters! The final value in a given
> pixel depends on the last region that covers it.
>
> Is this behaviour the same in all software? I would not bet on it. It is
> particularly ambiguous if you are just looking at the regions in DS9, as to
> whether the include or exclude region takes precedence. To be completely
> safe, create regions such that all exclude regions come after all include
> regions.
>
> The same consideration should be given to region type. Stick to circle, box,
> and ellipse to be safe.

Create some background regions and extract spectra from both A and B modules,
e.g. I selected three background regions in DS9 and saved them in ds9 fk5
format, `bgd1.reg`, `bgd2.reg`, and `bgd3.reg`.

```bash
# New syntax for getspecnoarf.py
# In event_cl/
getspecnoarf.py nu90201039002A01_cl.evt reg=bgd/bgd1.reg \
    indir=. outdir=bgd outprefix=bgd1A \
    attfile=../auxil/nu90201039002_att.fits.gz >& bgd/bgd1A.log

getspecnoarf.py nu90201039002A01_cl.evt reg=bgd/bgd2.reg \
    indir=. outdir=bgd outprefix=bgd2A \
    attfile=../auxil/nu90201039002_att.fits.gz >& bgd/bgd2A.log

getspecnoarf.py nu90201039002A01_cl.evt reg=bgd/bgd3.reg \
    indir=. outdir=bgd outprefix=bgd3A \
    attfile=../auxil/nu90201039002_att.fits.gz >& bgd/bgd3A.log

getspecnoarf.py nu90201039002B01_cl.evt reg=bgd/bgd1.reg \
    indir=. outdir=bgd outprefix=bgd1B \
    attfile=../auxil/nu90201039002_att.fits.gz >& bgd/bgd1B.log

getspecnoarf.py nu90201039002B01_cl.evt reg=bgd/bgd2.reg \
    indir=. outdir=bgd outprefix=bgd2B \
    attfile=../auxil/nu90201039002_att.fits.gz >& bgd/bgd2B.log

getspecnoarf.py nu90201039002B01_cl.evt reg=bgd/bgd3.reg \
    indir=. outdir=bgd outprefix=bgd3B \
    attfile=../auxil/nu90201039002_att.fits.gz >& bgd/bgd3B.log
```

This creates spectral product files, among which are the ungrouped PI spectra
files named like `bgd1A_sr.pha`, grouped spectra files `bgd1A_sr_g30.pha`, and
the response matrix file `bgd1A_sr.rmf`. We will use these files for the next
step.


Extract the spectrum of an extended source in the aperture defined by src.reg.

```bash
# In event_cl/
mkdir spec

getspecnoarf.py nu90201039002A01_cl.evt reg=src1.reg \
    indir=. outdir=spec outprefix=src1A \
    attfile=../auxil/nu90201039002_att.fits.gz >& spec/src1A.log

getspecnoarf.py nu90201039002B01_cl.evt reg=src.reg \
    indir=. outdir=spec outprefix=src1B \
    attfile=../auxil/nu90201039002_att.fits.gz >& spec/src1B.log
```



#### Fix the spectral products' RESPFILE keywords (optional)

Spectral products from the old version of getspecnoarf.py wrote relative paths
in the RESPFILE keyword of the PHA file, so that the latter must be loaded
from the same relative path for XSPEC to find the RMF file. You can fix the
RESPFILE keyword in the PHA files to remove the directory name, which lets you
load the spectrum when working in the same directory as it.

```bash
# In event_cl/ or where the spectral files are
find . -iname "*.pha" -type f -exec phafix.py {} \;
```

---

## 3. Create and fit the background model

### 3.1 Instrument maps

Create image masks for the detectors. Do this for each module, `nu*A01_cl.evt`
and `nu*B01_cl.evt`.

Command line:
```bash
# In event_cl/

nuskybgd mkinstrmap nu90201039002A01_cl.evt

nuskybgd mkinstrmap nu90201039002B01_cl.evt
```

Python:
```python
from nuskybgd.cli import mkinstrmap
mkinstrmap(['mkinstrmap', 'nu90201039002A01_cl.evt'])
mkinstrmap(['mkinstrmap', 'nu90201039002B01_cl.evt'])
```

This creates the files `newinstrmapA.fits` and `newinstrmapB.fits`, which are
image masks for the detectors.


### 3.2 Aspect histogram images

Make images of the 2D histogram of the pointing position, one for each module.

Command line:
```bash
# In event_cl/

nuskybgd aspecthist nu90201039002A_det1.fits gtifile=nu90201039002A01_gti.fits \
    out=aspecthistA.fits

nuskybgd aspecthist nu90201039002B_det1.fits gtifile=nu90201039002B01_gti.fits \
    out=aspecthistB.fits
```

Python:
```python
from nuskybgd.cli import aspecthist
aspecthist(['aspecthist', 'nu90201039002A_det1.fits',
    'gtifile=nu90201039002A01_gti.fits', 'out=aspecthistA.fits'])
aspecthist(['aspecthist', 'nu90201039002B_det1.fits',
    'gtifile=nu90201039002B01_gti.fits', 'out=aspecthistB.fits'])
```

This creates the files `aspecthistA.fits` and `aspecthistB.fits` in the
directory `event_cl/`. They are images representing a 2D histogram in time of
the pointing position.


### 3.3 Background aperture images

Create images of the aperture background model and detector mask convolved
with the aspect. For each module, one image is created for the aperture
background and four images for the detector masks.

Command line:
```bash
# In event_cl/bgd/

nuskybgd projbgd refimg=../imB3to20keV.fits out=bgdapA.fits \
	mod=A det=1234 chipmap=../newinstrmapA.fits aspect=../aspecthistA.fits

nuskybgd projbgd refimg=../imB3to20keV.fits out=bgdapB.fits \
	mod=B det=1234 chipmap=../newinstrmapB.fits aspect=../aspecthistB.fits
```

Python:
```python
from nuskybgd.cli import projbgd
projbgd(['projbgd', 'refimg=../imB3to20keV.fits', 'out=bgdapA.fits',
    'mod=A', 'det=1234', 'chipmap=../newinstrmapA.fits',
    'aspect=../aspecthistA.fits'])
projbgd(['projbgd', 'refimg=../imB3to20keV.fits', 'out=bgdapB.fits',
    'mod=B', 'det=1234', 'chipmap=../newinstrmapB.fits',
    'aspect=../aspecthistB.fits'])
```

This creates the files `bgdapA.fits` and `bgdapB.fits`, which are the aperture
background images rotated and convolved with the aspect histogram images, and
`det0Aim.fits`, `det1Aim.fits`, `det2Aim.fits`, `det3Aim.fits`,
`det0Bim.fits`, `det1Bim.fits`, `det2Bim.fits`, and `det3Bim.fits`, which are
the detector masks rotated and convolved with the aspect histogram images. The
files are in the directory `event_cl/bgd/`.



### 3.4 The background model

Run `nuskybgd fit` (requires PyXspec) to create an XSPEC save file
`bgdparams.xcm`, which contains the background models with preset
normalizations.

First, create a file `bgdinfo.json` in `bgd/` with the following structure
(`nuskybgd fit --help` will print this template). The file names for the
images are default values so you may not need to modify them, but the files
for the background regions need to be updated for your data. The `"regfiles"`
list must correspond to the files in the `"bgfiles"` list.

The setting "fcxb_config" "links" are used to tie Xspec model parameters for
the FCXB model normalization. These are comprised of a list of two values for
the spectrum numbers to tie together. In the example below, all of the spectra
for FPMB are tied to the spectra for identical regions from FPMA. Note that
the Xspec spectrum number starts at 1. Without this setting, all of the FCXB
normalizations are set to free.

```json
{
    "bgfiles": [
        "bgd1A_sr_g30.pha", "bgd1B_sr_g30.pha",
        "bgd2A_sr_g30.pha", "bgd2B_sr_g30.pha",
        "bgd3A_sr_g30.pha", "bgd3B_sr_g30.pha"
    ],

    "regfiles": [
        "bgd1.reg", "bgd1.reg",
        "bgd2.reg", "bgd2.reg",
        "bgd3.reg", "bgd3.reg"
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
    }
}
```

Then, run `nuskybgd fit`, directing stdout to a log file. Check the log to see
if the task encountered any problems. Most of the logged output comes from
Xspec.

Command line:
```bash
# In event_cl/bgd/
nuskybgd fit bgdinfo.json savefile=IC342bgd >& fit.log
```

Python:
(Warning for interactive use: Xspec will flood the terminal with messages)
```python
from nuskybgd.cli import fit
fit(['fit', 'bgdinfo.json', 'savefile=IC342bgd'])
```

After successfully running `nuskybgd fit`, the save file (in this case,
`IC342bgd.xcm`) can be loaded in Xspec to recreate the same state. From there,
the user can examine and tweak with the model. They can also save a modified
version of it, preferably under a different name, to experiment with for the
next step.

```
# Start xspec, then input these commands
@IC342bgd.xcm
ignore **:**-3. 150.-**
cpd /xw
setplot energy
setplot command res y 1e-4 0.04
plot ldata delchi
```

The generated background model then needs to be fitted to the background
spectra.

The first thing the user should check is whether the preset normalizations
provide a close fit to the background spectra. If the background regions do
not have extended emission, the user may be able to obtain a good background
model by simply running `fit` in Xspec.

> The save file contains settings for Xspec fit equivalent to the following:
> `statistic chi; method leven 30000 1e-4; ignore **:**-3. 150.-**`
> The user should change this setting as required.

On the other hand, if the background regions do contain extended emission,
simply running `fit` may not be enough, because the excess emission from the
real source cannot be accounted for by the preset model components. In that
case, the user must use discretion to account for the additional source(s).
Due to the number of free parameters in the nuskybgd model and possible
degeneracies with any added source models, it may be necessary to hand-tune
the fit and not fit everything simultaneously.

> Do not remove any of the nuskybgd generated models or change their names.

When the model is OK, write the current state to a save file under a different
name to the nuskybgd generated save file. This save file contains model
parameters for the background components. The subsequent tasks will look for
these using the model component names.

```
# Creates a save file mymodel.xcm
save all mymodel
```

The saved xcm file should not contain any general XSPEC/Tcl scripts because
PyXspec will not properly execute it.

---

## 4. Applying the fitted background


### 4.1 Spectral analysis --- calculate the background model for a source region

To calculate the amount of background in some source region, create a region
file and extract the spectrum from it. For this example, the region file is
named `src1.reg` and the grouped spectrum file is `src1_g30.pha`. Use
`nuskybgd spec` on the best fit model save file.

```bash
nuskybgd spec bgdinfo.json mymodel.xcm src1.reg src1_g30.pha
```

See the help message of `nuskybgd spec` for what it does. This step creates a
new Xspec save file, `bgd_src1.xcm` (for custom output file name, use the
optional argument `savefile=` to specify it.)

This save file is very similar to the input, with the addition of the source
spectrum, and the background model parameter values have been calculated for
the source spectrum region. The background model parameters are free, and the
user can verify the quality of the background subtraction and tweak the model.

When the background model is satisfactory, create an Xspec save file with the
command `save all bgd_src1_mymodel` and use this for the following step.
Otherwise, just use the output from `nuskybgd spec` (as below).

```bash
nuskybgd simplify bgdinfo.json bgd_src1.xcm
```

This final step removes the spectra from the background regions from the saved
Xspec state, leaving only the source spectrum. The background model parameters
are frozen, and the user can load this save file in Xspec and proceed to
adding their source model. The default output file name for the example
command above is `bgd_only_src1.xcm`, but this can be customized using the
`savefile=` argument.


### 4.2 Generate background images

Images of the model background sources can be created from the Xspec save file
obtained in the previous step. For example, to generate images of the model
counts between 3 and 20 keV,

```bash
nuskybgd image bgdinfo.json mymodel.xcm 3 20
```

This results in two files for each background source (A and B), as well as
images of the FCXB in each background spectrum region.


---

## Add detabs to rmf

Creates det0A.rmf, det1A.rmf, etc... in bgd/

```
cd /Users/qw/astro/nustar/IC342_X1/90201039002/event_cl
absrmf.py nu90201039002A01_cl.evt bgd/det
absrmf.py nu90201039002B01_cl.evt bgd/det
```



## Create fake spectra of apbgd and fcxb components

```
cd bgd/
imrefspec.py AB 0123
```


---


TODO

Lookup nucoord source to figure out coord transform, make instrmap use the
same method.

hard background modelling examples
ophiuchus -- cluster fills fov
A2146 -- gain shifts


FCXB choice for users whether to tie some regions/how

Create hardcoded table model of the total background to read in to Xspec?





---


