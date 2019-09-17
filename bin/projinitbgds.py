#!/usr/bin/env python3
"""
Project the instrument map onto sky coordinates to get the projected detector
mask images and background aperture images.

Usage: projinitbgds refimg=flux.fits out=output.fits mod=[A,B] det=[1,2,3,4] \\
           chipmap=chipmap.fits aspect=aspect.fits

Example:

projinitbgds refimg=imA4to25keV.fits out=bgdapA.fits mod=A det=1234 \\
    chipmap=newinstrmapA.fits aspect=projobsA.fits

The output file bgdapA.fits has the background aperture image in sky
coordinates. A series of files showing the projected detector masks for each
detector specified by det= will also be created, named det[0-3]Aim.fits.
"""
import sys
import os
import astropy.io.fits as pf
import numpy as np
from scipy.ndimage import affine_transform
from scipy.ndimage import shift
from numpy.linalg import inv


_NUSKYBGD_AUX_NAME = 'NUSKYBGD_AUXIL'


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

    auxildir = os.environ[_NUSKYBGD_AUX_NAME]
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


if __name__ == '__main__':

    args = {}
    keys = ('refimg', 'out', 'det', 'mod', 'chipmap', 'aspect')
    # TODO: validate input chipmap and aspect files
    for _ in sys.argv[1:]:
        pars = _.split('=')
        if pars[0] in keys:
            args[pars[0]] = pars[1]

    for _ in keys:
        if _ not in args:
            print('Missing %s= argument' % _)
            print(__doc__)
            sys.exit(1)

    try:
        refimg = pf.open(args['refimg'])
    except FileNotFoundError:
        print('refimg file not found')
        sys.exit(1)

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
        sys.exit(1)
    except ValueError:
        print('Unexpected refimg format')
        sys.exit(1)

    if args['mod'] not in ('A', 'B'):
        print('mod= must be set to A or B')
        sys.exit(1)

    det_input = []
    for _ in '1234':
        if _ in args['det']:
            det_input.append(int(_))
    if len(det_input) == 0:
        print('No valid detector number specified.')
        sys.exit(1)

    if _NUSKYBGD_AUX_NAME not in os.environ:
        print('Environment variable $%s must be set.' % _NUSKYBGD_AUX_NAME)
        sys.exit(1)

    auxildir = os.environ[_NUSKYBGD_AUX_NAME]

    if not os.path.exists(auxildir):
        print('Error: NUSKYBGD_AUXIL path does not exists (%s)' % auxildir)

    if args['out'][0] != '!':
        halt = False
        if os.path.exists(args['out']):
            print('Output file %s exists.' % args['out'])
            halt = True
        for idet in np.arange(4):
            _ = 'det%d%sim.fits' % (idet, args['mod'])
            if os.path.exists(_):
                print('Output file %s exists.' % _)
                halt = True
        if halt:
            sys.exit(1)

    print('Using auxiliary data in %s' % auxildir)

    detmask = get_det_mask(args['chipmap'], det_input)
    apim = get_aperture_image(args['mod'])
    posim, pos_xoff, pos_yoff = get_aspect_hist_image(args['aspect'])
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
    outfile.writeto(args['out'], overwrite=True)

    # Write detector mask images
    for idet_input in range(len(det_input)):
        pf.HDUList([pf.PrimaryHDU(bgdimi[idet_input], header=hdr)]).writeto(
            'det%d%sim.fits' % (det_input[idet_input] - 1, args['mod']),
            overwrite=True)  # For det number, input 1-4 -> output 0-3

    sys.exit(0)
