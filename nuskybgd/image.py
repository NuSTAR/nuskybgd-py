# (c) 2019 Qian Wang
# This file is part of nuskybgd released under MIT License.
# See LICENSE file for full details.

import numpy as np
import astropy.io.fits as pf
from scipy.ndimage import affine_transform
from scipy.ndimage import shift
from numpy.linalg import inv


def make_aspecthist_img(coords, dt):

    asphistimg = np.zeros((1000, 1000), dtype=np.float64)
    x_min, y_min, x_max, y_max = 999, 999, 0, 0
    for i in range(len(dt)):
        x, y = int(np.floor(coords[i][0])), int(np.floor(coords[i][1]))
        asphistimg[y, x] += dt[i]
        if y < y_min:
            y_min = y
        elif y > y_max:
            y_max = y
        if x < x_min:
            x_min = x
        elif x > x_max:
            x_max = x

    print('Image subspace: ', x_min, x_max, y_min, y_max)

    return asphistimg, x_min, x_max, y_min, y_max


def write_aspecthist_img(asphistimg, outfilename, aimpointext, offsets,
                         overwrite=False):
    import datetime

    x_min, x_max, y_min, y_max = offsets

    # Write some useful info into the header
    aspecthistfh = pf.HDUList(
        pf.PrimaryHDU(asphistimg[y_min:y_max + 1, x_min:x_max + 1])
    )
    aspecthistfh[0].header['EXTNAME'] = 'ASPECT_HISTOGRAM'
    aspecthistfh[0].header['X_OFF'] = (x_min, 'x offset in pixels')
    aspecthistfh[0].header['Y_OFF'] = (y_min, 'y offset in pixels')
    aspecthistfh[0].header['EXPOSURE'] = (
        np.float32(np.sum(asphistimg[y_min:y_max + 1, x_min:x_max + 1])),
        'seconds, total exposure time'
    )
    aspecthistfh[0].header['COMMENT'] = (
        'Add the x/y offset to image coordinates to recover aimpoint '
        'in detector coordinates.'
    )
    aspecthistfh[0].header['DATE'] = (
        datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%S'),
        'File creation date (YYYY-MM-DDThh:mm:ss UTC)'
    )
    for hdrkey in ('TELESCOP', 'INSTRUME', 'OBS_ID', 'OBJECT', 'TARG_ID',
                   'RA_OBJ', 'DEC_OBJ', 'RA_NOM', 'DEC_NOM',
                   'RA_PNT', 'DEC_PNT', 'TIMESYS', 'MJDREFI', 'MJDREFF',
                   'CLOCKAPP', 'TSTART', 'TSTOP', 'DATE-OBS', 'DATE-END'):
        if hdrkey in aimpointext.header:
            aspecthistfh[0].header[hdrkey] = (
                aimpointext.header[hdrkey],
                aimpointext.header.comments[hdrkey]
            )

    aspecthistfh.writeto(outfilename, overwrite=overwrite)


def get_det_mask(instmap_mask_file, detnum):
    """
    Return an image of the instrument map mask for the specified detector(s).
    Active area is set to 1, the rest 0. Require instmap_mask_file to be str.
    Require detnum to be int or [int].

    Usage: get_det_mask(instmap_mask_file, detnum)

    Output: mask[det, y, x] is a 3D numpy array, in which mask[i] are detector
    masks for each detector.

    TODO: add input validation.
    """
    if isinstance(detnum, int):
        detnum = [detnum]

    mskfile = pf.open(instmap_mask_file)
    shape = mskfile['INSTRMAP'].data.shape
    mask = np.zeros((len(detnum), shape[0], shape[1]), dtype=np.float64)
    for i in np.arange(len(detnum)):
        inx = np.where(mskfile['INSTRMAP'].data == detnum[i])
        mask[i, inx[0], inx[1]] = 1.
    return mask


def get_aspect_hist_image(asphistimgfile):
    """
    Return an image of the aspect histogram and x and y offsets if any (from
    the X_OFF and Y_OFF keywords in header).
    """
    fh = pf.open(asphistimgfile)
    # Check if the image has x, y offsets
    if ('X_OFF' in fh[0].header and 'Y_OFF' in fh[0].header):
        x_off = fh[0].header['X_OFF']
        y_off = fh[0].header['Y_OFF']
    else:
        x_off, y_off = 0, 0
    return fh[0].data, x_off, y_off


def get_aspect_hist_peak(asphistimg, xoff=0, yoff=0):
    """
    Return the (x, y) array indices of the maximum in the aspect histogram
    image.

    TODO: what to do if there are two such pixels
    """
    peak = np.argmax(asphistimg)
    return (peak % asphistimg.shape[1] + xoff,
            peak // asphistimg.shape[1] + yoff)


def get_aperture_image(detector):
    from . import conf
    import os

    auxildir = os.environ[conf._AUX_ENV]
    if detector not in ('A', 'B'):
        print('Detector must have value A or B.')
        return False

    # Look up a preset params file to determine shift and rotation, rescale
    # factor
    with open('%s/nomapbgdparams.dat' % auxildir, 'r') as pfile:
        _content = pfile.read().splitlines()
        # Extract the first part before the # on each line
        paperbgd = [float(_.split('#')[0]) for _ in _content]

    if detector == 'A':
        aperture_db = pf.open('%s/detA_det1.img.gz' % auxildir)[0].data
        shift_x = paperbgd[4] / 0.12096
        shift_y = paperbgd[5] / 0.12096
        rot_angle = paperbgd[6]
        rescale = 3.1865E-05 * paperbgd[0] * 0.12096 * 0.12096
    else:
        aperture_db = pf.open('%s/detB_det1.img.gz' % auxildir)[0].data
        shift_x = paperbgd[7] / 0.12096
        shift_y = paperbgd[8] / 0.12096
        rot_angle = paperbgd[9]
        rescale = 3.1865E-05 * paperbgd[1] * 0.12096 * 0.12096

    aperture = rescale * transform_image(
        aperture_db, shift_x, shift_y, rot_angle)

    aperture += 0.00495 * 0.7

    # Crop to 360 x 360
    return 1. * aperture[108:468, 108:468]


def transform_image(image, shiftx, shifty, angle, order=1):
    """
    Apply shift and rotation to the image. The translation is applied first,
    then the rotation.

    Usage: transform_image(image, dx, dy, angle)

    The rotation angle is in radians and the translation is in pixel.

    The transformation is implemented as sequence of affine transformations.
    The scipy module takes a matrix of the form (ndim + 1, ndim + 1), where it
    assumes that the transformation is specified using homogeneous
    coordinates. This matrix has the 2x2 rotation matrix in the top left
    corner, and the linear shifts in the top right. They are applied in this
    order:

    1  0  shiftx
    0  1  shifty
    0  0  1
    (translation by shift amounts)

    1  0  -(X-1)/2
    0  1  -(Y-1)/2
    0  0  1
    (translation to center rotation on the IDL rot center)

    cos  sin 0
    -sin cos 0
    0    0   1
    (clockwise rotation)

    1  0  +(X-1)/2
    0  1  +(Y-1)/2
    0  0  1
    (undo translation for center of rotation)
    """

    if shiftx == 0 and shifty == 0 and angle == 0:
        return image

    elif angle == 0:
        return shift(image, [shifty, shiftx],
                     order=order, mode='wrap', prefilter=False)

    else:
        # The original IDL implementation performs the linear translation
        # first (wrapping at borders), then rotates the image clockwise by an
        # angle in degrees. The center of the rotation is (X-1)/2, (Y-1)/2. In
        # both steps, bilinear interpolation is used.

        # Numpy array coordinates are (y, x). This swaps dx and dy in the
        # translation part, and the position of the -sin element in the
        # rotation part, compared to the standard version for (x, y, 1).
        # Beware that the coordinate transforms are calculated in the
        # conventional (x, y) sense, but are written in numpy's (y, x) order
        # when implementing them in the transformation matrix.
        cx, sx = np.cos(angle), np.sin(angle)

        # Center of rotation
        rot_x = 0.5 * (image.shape[1] - 1)
        rot_y = 0.5 * (image.shape[0] - 1)

        dx = cx * (shiftx - rot_x) + sx * (shifty - rot_y) + rot_x
        dy = -sx * (shiftx - rot_x) + cx * (shifty - rot_y) + rot_y

        tx = np.array([[cx, -sx, dy],
                       [sx, cx, dx],
                       [0, 0, 1]], dtype=np.float64)

        # print(cx, sx)
        # print(tx)

        # The prefilter option, which is turned on by default, applies an
        # image sharpening. It must not be applied. The mode is set to 'wrap'
        # to emulate the behavior of the original implementation of image
        # shift. The spline interpolation order is 1 for bilinear, 3 for
        # bicubic (bilinear is the original behavior).
        return affine_transform(image, inv(tx),
                                order=order, mode='wrap', prefilter=False)


def crop_image(image):
    """
    Crop an image with 0 padding so that there is no empty margin. Return the
    cropped image, and the x, y offsets of the new image coordinates relative
    to the uncropped image. If there is no non-zero value in the image,
    returns False.

    Usage: crop_image(image)

    Return: (cropped image, [offset_x, offset_y])
    """
    _ = np.abs(image)
    filled_x = np.where(np.sum(_, axis=0) != 0)[0]
    filled_y = np.where(np.sum(_, axis=1) != 0)[0]

    if len(filled_x) == 0:
        # Empty image
        print('The image has no non-zero pixels')
        return False

    min_x, max_x = np.min(filled_x), np.max(filled_x)
    min_y, max_y = np.min(filled_y), np.max(filled_y)

    return image[min_y:max_y + 1, min_x:max_x + 1], [min_x, min_y]


def get_evtfile_pa(evtfile):
    """
    Get position angle from events file header.
    """
    hdus = pf.open(evtfile)
    pa = hdus['EVENTS'].header['PA_PNT'] + 1.0
    return pa


def calc_pa_to_arot(pa, old_pa=False):
    if not old_pa:
        arot = - pa * np.pi / 180.
    else:
        arot = (pa - 90.) * np.pi / 180.
    return arot


def apply_badpix(instrmap, bpixexts, pixmap, detnum):
    DETNAM = {'DET0': 0, 'DET1': 1, 'DET2': 2, 'DET3': 3}

    output = 1. * instrmap

    for ext in bpixexts:
        hh = ext.header
        if ('TSTART' in hh and 'TSTOP' in hh):
            tobs = hh['TSTOP'] - hh['TSTART']
        else:
            tobs = None

        if not ('EXTNAME' in ext.header and
                'BADPIX' in ext.header['EXTNAME'] and
                'DETNAM' in ext.header and
                ext.header['DETNAM'] in DETNAM):
            continue

        idet = DETNAM[ext.header['DETNAM']]
        # print(idet)
        badpix = ext.data

        x = badpix['RAWX']
        y = badpix['RAWY']
        flags = badpix['BADFLAG']
        if tobs is not None:
            dt = badpix['TIME_STOP'] - badpix['TIME']

        for i in np.arange(len(x)):
            if ((tobs is not None and dt[i] > 0.8 * tobs) or
                    flags[i][-2] is True):
                ii = np.where((pixmap == x[i] + y[i] * 32) *
                              (detnum == idet))
                output[ii] = 0

            ii = np.where(output > 0)
            output[ii] = detnum[ii] + 1

    return output


def get_badpix_exts(bpixfiles):
    DETNAM = {'DET0': 0, 'DET1': 1, 'DET2': 2, 'DET3': 3}

    # Collect list of bad pixel extensions, group into FPMA and B
    bpixlist = {'A': [], 'B': []}

    for _ in bpixfiles:
        ff = pf.open(_)
        # Check FPM
        phdr = ff[0].header
        try:
            if phdr['TELESCOP'] != 'NuSTAR':
                print('%s: TELESCOP != NuSTAR, skipping.' % _)
                continue
        except KeyError:
            print('%s: no TELESCOP keyword, skipping.' % _)
            continue

        try:
            if phdr['INSTRUME'] == 'FPMA':
                ab = 'A'
            elif phdr['INSTRUME'] == 'FPMB':
                ab = 'B'
            else:
                print('%s: unknown INSTRUME: %s' % (_, phdr['INSTRUME']))
                continue
        except KeyError:
            print('%s: no INSTRUME keyword, skipping.' % _)
            continue

        for ext in ff[1:]:
            try:
                ehdr = ext.header
                if not (ehdr['XTENSION'] == 'BINTABLE' and
                        'BADPIX' in ehdr['EXTNAME']):
                    continue
                if ehdr['DETNAM'] in DETNAM:
                    bpixlist[ab].append(ext)
                    print('Added %s %s %s to list.' % (_, ab, ehdr['DETNAM']))
                else:
                    print('%s: unknown DETNAM: %s' % (_, ehdr['DETNAM']))
                    continue
            except KeyError:
                # Extension does not have XTENSION, EXTNAME or DETNAM
                continue

    return bpixlist


def make_det_mask(evthdr):
    import os
    from . import caldb
    from . import conf

    # grade weighting from NGC253 002 obs.
    GRADE_WT = [1.00000, 0.124902, 0.117130, 0.114720, 0.118038,
                0.0114296, 0.0101738, 0.0113617, 0.0122017, 0.0157910,
                0.0144079, 0.0145691, 0.0149934, 0.00165462, 0.00194312,
                0.00156128, 0.00143400, 0.00210433, 0.00180735, 0.00140006,
                0.00169704, 0.00189220, 0.00160371, 0.00150188, 0.00168007,
                0.000296983, 0.000364864]

    fpm = evthdr['INSTRUME']
    obsutctime = evthdr['DATE-OBS']
    cal = caldb.CalDB(os.environ[conf._CALDB_ENV])
    pixpospath = cal.getPIXPOS(fpm, 'DET0', obsutctime)

    pixposf = pf.open('%s/%s' % (os.environ[conf._CALDB_ENV], pixpospath))

    pixmap = np.full((360, 360), -1, dtype=np.int32)
    detnum = np.full((360, 360), -1, dtype=np.int32)
    allpdf = np.zeros((360, 360), dtype=np.float64)

    for ext in pixposf:
        if (('EXTNAME' not in ext.header) or
            ('PIXPOS' not in ext.header['EXTNAME']) or
                ('DETNAM' not in ext.header)):
            continue
        idet = int(ext.header['DETNAM'].replace('DET', ''))
        pixpos = ext.data

        # Store references to columns
        pp_det1x = pixpos['REF_DET1X']
        pp_det1y = pixpos['REF_DET1Y']
        pp_rawx = pixpos['RAWX']
        pp_rawy = pixpos['RAWY']
        pp_grade = pixpos['GRADE']
        pp_pdf = pixpos['PDF']

        for ix in np.arange(32):
            for iy in np.arange(32):
                # Get array indices where all of the following are True
                ii = np.where((pp_det1x != -1) *
                              (pp_rawx == ix) *
                              (pp_rawy == iy) *
                              (pp_grade <= 26))[0]

                thispdf = np.zeros((360, 360), dtype=np.float64)

                for i in ii:
                    if not np.isnan(pp_pdf[i]).any():
                        # No nan value in PDF
                        ref_x = pp_det1x[i]
                        ref_y = pp_det1y[i]
                        thispdf[ref_y:ref_y + 7, ref_x:ref_x + 7] += (
                            pp_pdf[i] *
                            GRADE_WT[pp_grade[i]])

                ii = np.where(thispdf > allpdf)
                if len(ii) > 0:
                    allpdf[ii] = thispdf[ii]
                    pixmap[ii] = ix + iy * 32
                    detnum[ii] = idet

    pixmap = shift(pixmap, [-1, -1], mode='wrap', prefilter=False, order=1)
    detnum = shift(detnum, [-1, -1], mode='wrap', prefilter=False, order=1)

    return pixmap, detnum


def get_caldb_instrmap(evthdr):
    """
    Get instrument map image from CALDB, shifted by (-1, -1).

    Returns (image, FITS header).
    """
    import os
    from . import conf
    from . import caldb
    from . import util

    fpm = evthdr['INSTRUME']
    obsutctime = evthdr['DATE-OBS']

    cal = caldb.CalDB(os.environ[conf._CALDB_ENV])
    imappath = cal.getINSTRMAP(fpm, obsutctime)
    imapf = pf.open('%s/%s' % (os.environ[conf._CALDB_ENV], imappath))
    try:
        imap = imapf['INSTRMAP'].data
    except KeyError:
        print('Error: INSTRMAP extension not found.')
        return False
    hdr = imapf['INSTRMAP'].header
    hdr['HISTORY'] = 'Created inst. map from %s' % os.path.basename(imappath)
    util.hdu_timestamp_write(imapf['INSTRMAP'])

    return imap, hdr
