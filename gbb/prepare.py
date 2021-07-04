# -*- coding: utf-8 -*-
"""Utility functions for preprocessing epi time series before use in GBB. These
functions call afni and fsl programs which need to be installed and accessible
from the command line.
"""

# python standard library inputs
import os
import sys
import subprocess

# external inputs
import numpy as np
import nibabel as nb
from sh import gunzip
from nighres.intensity import phase_unwrapping
from scipy.ndimage import gaussian_filter

# local inputs
from .io import get_filename, load_volume, save_volume
from .plot import plot_moco

__all__ = ['preprocess_epi', 'motion_correction', 'allineate', 'susan_filter']


def preprocess_epi(file_magn, file_phase, phase_clip=0.25, std_max=0.25,
                   sigma=10.0):
    """Computation of contrast-enhanced mean from magnitude and phase data of a
    functional time series. The magnitude time series is realigned using afni's
    `3dvolreg` and the phase time series is unwrapped. Then, the realignment
    parameters are applied to the unwrapped phase time series. Motion outliers
    are removed from realigned time series. The mean across time is computed for
    both time series and the phase data is rescaled to the range [0, 1].
    Locations in the mean phase image are smoothed out (gaussian filter) which
    show high temporal variance. This tries to decrease the impact of
    singularities to the final image. The enhanced mean epi is then created by
    multiplying the mean magnitude image by the mean phase image.

    Parameters
    ----------
    file_magn : str
        File name of 4d nifti volume (magnitude data).
    file_phase : str
        File name of 4d nifti volume (phase data).
    phase_clip : float
        Clip mean phase data.
    std_max : float
        Temporal variance threshold to smooth out singularities from the mean
        phase image.
    sigma : float
        Standard deviation of gaussian kernel to smooth out singularities from
        the mean phase image.

    Returns
    -------
    None

    """

    # get file parts of magnitude and phase data
    path_m, name_m, ext_m = get_filename(file_magn)
    path_p, name_p, ext_p = get_filename(file_phase)

    # file names for intermediate and final files
    file_phase_unwrap = os.path.join(path_p, name_p + "_unwrap" + ext_p)
    file_umagn = os.path.join(path_m, "u" + name_m + ext_m)
    file_uphase = os.path.join(path_p, "u" + name_p + ext_p)
    file_umagn_m = os.path.join(path_m, "mean_u" + name_m + ext_m)
    file_uphase_m = os.path.join(path_p, "mean_u" + name_p + ext_p)
    file_umagn_w = os.path.join(path_m, "mean_u" + name_m + "_enhanced" + ext_m)

    # motion correction of magnitude data
    res_moco = motion_correction(file_magn)

    # unwrap phase data
    arr_p, affine_p, header_p = load_volume(file_phase)
    arr_p_unwrap = np.zeros_like(arr_p)
    for i in range(np.shape(arr_p)[3]):
        print("Unwrap volume " + str(i))
        tmp_in = nb.Nifti1Image(arr_p[:, :, :, i], affine_p, header_p)
        tmp_out = phase_unwrapping(tmp_in, mask=None, nquadrants=3,
                                   tv_flattening=True, tv_scale=0.5,
                                   save_data=False, overwrite=False,
                                   output_dir=False, file_name=False)

        arr_p_unwrap[:, :, :, i] = tmp_out["result"].get_fdata()

    # save unwrapped time series
    save_volume(file_phase_unwrap, arr_p_unwrap, affine_p, header_p)

    # apply motion correction to unwrapped phase time series
    allineate(file_phase_unwrap, file_uphase, res_moco["matrix"])

    # load realigned timeseries
    arr_m, affine, header_m = load_volume(file_umagn)
    arr_p, _, header_p = load_volume(file_uphase)

    # remove outliers from magnitude and phase data
    t = np.loadtxt(res_moco["outlier"]).astype(int)
    arr_m = arr_m[:, :, :, t == 0]
    arr_p = arr_p[:, :, :, t == 0]

    # get magnitude and phase mean (plus variance)
    arr_m_mean = np.mean(arr_m, axis=3)
    arr_p_mean = np.mean(arr_p, axis=3)
    arr_p_std = np.std(arr_p, axis=3)

    # threshold phase data
    frac_max = phase_clip * np.max(arr_p_mean)
    arr_p_mean[arr_p_mean > frac_max] = frac_max
    arr_p_mean[arr_p_mean < -frac_max] = -frac_max

    # rescale and invert
    phase_range = np.max(arr_p_mean) - np.min(arr_p_mean)
    arr_p_mean = (arr_p_mean - np.min(arr_p_mean)) / phase_range
    arr_p_mean = 1 - arr_p_mean

    save_volume(file_umagn_m, arr_m_mean, affine, header_m)
    save_volume(file_uphase_m, arr_p_mean, affine, header_p)

    # create binary mask by thresholding
    arr_mask = arr_p_std.copy()
    arr_mask[arr_mask < std_max] = 0
    arr_mask[arr_mask != 0] = 1

    # filter phase
    arr_p_mean = _fill_with_gaussian(arr_p_mean, arr_mask, sigma)

    # weight magnitude image (enhanced contrast)
    arr_m_mean *= arr_p_mean
    save_volume(file_umagn_w, arr_m_mean, affine, header_m)

    # delete unaligned unwrapped phase timeseries
    os.remove(file_phase_unwrap)


def motion_correction(file_in):
    """Simple wrapper of afni's `3dvolreg` program to do motion correction of
    functional time series. Transformation matrices (target to volume) and
    motion parameters for each volume are saved to disk. Motion parameters are
    saved as 6 ASCII formatted columns (roll, pitch, yaw, dS, dL, dP) with the
    following convention:

    * roll: rotation about the I-S axis (degrees ccw)
    * pitch: rotation about the R-L axis (degrees ccw)
    * yaw: rotation about the A-P axis (degrees ccw)
    * dS: displacement in the superior direction (mm)
    * dL: displacement in the left direction (mm)
    * dP: displacement in the posterior direction (mm)

    Plots are saved in svg format to illustrate translational and rotational
    movements. Furthermore, a regressor of no interest is saved which indicates
    all volumes which were classified as motion outliers. Motion outliers are
    those volumes which exceed a defined maximal translation and/or rotation on
    a volume-to-volume (short) adn/or volume-to-volume0 (long) basis. The
    threshold for translational movements is taken as the smallest voxel edge
    length. The threshold for rotational movements is hard coded.

    Parameters
    ----------
    file_in : str
        File name of functional time series.

    Returns
    -------
    dict
        Dictionary collecting the output under the following keys

        * file_out (str) : File name of realigned functional time series.
        * params (str) : File name of saved motion parameters.
        * matrix (str) : File name of saved matrix transformations.
        * outlier (str) : File name of saved motion outliers.

    """

    path_in, name_in, ext_in = get_filename(file_in)
    path_moco = os.path.join(path_in, 'diagnosis')

    # make subfolders
    if not os.path.exists(path_moco):
        os.makedirs(path_moco)

    file_out = os.path.join(path_in, "u"+name_in+ext_in)
    file_params = os.path.join(path_moco, "moco_params.1D")
    file_matrix = os.path.join(path_moco, "moco_matrix.1D")
    file_outlier = os.path.join(path_moco, "outlier.1D")

    # run afni 3dvolreg
    # ---
    # Fourier: interpolation method
    # clipit: clips the values in each volume to be in the same range
    # zpad: zero pad around the edges by n voxels
    # float: for output to be written in floating point format
    # base: set the target volume to the nth volume
    # twopass: do two passes of the registration algorithm

    mc = '3dvolreg'
    mc += ' -Fourier -clipit -zpad 4 -float -base 0 -twopass'
    mc += ' -prefix ' + file_out
    mc += ' -1Dfile ' + file_params
    mc += ' -1Dmatrix_save ' + file_matrix
    mc += ' ' + file_in

    dash = '-' * 40
    print(dash)
    print('{}'.format('Run afni 3dvolreg'))
    print(dash)
    print(mc)

    try:
        subprocess.check_output(mc, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        print(e.output)

    # plot motion parameters
    plot_moco(file_params, os.path.join(path_moco, "plot_trans.svg"), "trans")
    plot_moco(file_params,  os.path.join(path_moco, "plot_rot.svg"), "rot")

    mp = np.loadtxt(file_params)
    _, _, header = load_volume(file_in)

    # define outlier thresholds
    trans_max = np.min(header["pixdim"][1:4])
    rot_max = 1.0

    # get euclidean distances for trans and rot
    trans = np.sqrt(mp[:, 3]**2+mp[:, 4]**2+mp[:, 5]**2)
    rot = np.sqrt(mp[:, 0]**2+mp[:, 1]**2+mp[:, 2]**2)

    # search for short (volume-to-volume) and long (volume-to-volume0) motion
    # and classify as outliers
    nt = header["dim"][4]
    outlier = np.zeros(nt)
    for i in range(nt-1):
        if trans[i+1] > trans_max:
            outlier[i+1] = 1
        elif rot[i+1] > rot_max:
            outlier[i+1] = 1
        elif trans[i+1] - trans[i] > trans_max / 2:
            outlier[i+1] = 1
        elif rot[i+1] - rot[i] > rot_max / 2:
            outlier[i+1] = 1

    # save outliers as regressor of no interest in ASCII format
    np.savetxt(file_outlier, outlier, fmt='%i')

    return {'res':  file_out,
            'params': file_params,
            'matrix': file_matrix,
            'outlier': file_outlier
            }


def allineate(file_in, file_out, file_matrix):
    """Simple wrapper of afni's `3dAllineate` program to apply a transformation
    matrix (1Dmatrix) to an input volume.

    Parameters
    ----------
    file_in : str
        File name of nifti volume.
    file_out : str
        File name of aligned nifti volume.
    file_matrix : str
        File name of transformation matrix (afni).

    Returns
    -------
    None

    """

    # run afni 3dAllineate
    # ---
    # final: interpolation method
    # nopad: do not use zero padding on the source image
    # float: for output to be written in floating point format

    aa = '3dAllineate'
    aa += ' -final linear -nopad -float'
    aa += ' -source ' + file_in
    aa += ' -prefix ' + file_out
    aa += ' -1Dmatrix_apply ' + file_matrix

    dash = '-' * 40
    print(dash)
    print('{}'.format('Run afni 3dAllineate'))
    print(dash)
    print(aa)

    try:
        subprocess.check_output(aa, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        print(e.output)


def susan_filter(file_in, file_out, sigma=3):
    """Simple wrapper of fsl's `susan` noise filtering program to apply an
    edge-preserving smoothing filter to an input volume.

    Parameters
    ----------
    file_in : str
        File name of nifti volume.
    file_out : str
        File name of filtered nifti volume.
    sigma : float
        spatial size (half-width) of smoothing in mm.

    Returns
    -------
    None

    """

    try:
        subprocess.run('susan', stdout=subprocess.DEVNULL)
    except FileNotFoundError:
        sys.exit("Could not find command 'susan'. Make sure FSL is installed "
                 "and can be accessed from the command line.")

    # run fsl susan
    # ---
    # bt_threshold: brightness threshold which should be greater than noise
    # level and less then contrast of edges to be preserved
    # dim: dimensionality (2 or 3)
    # use_median: use local median filter for single-point noise (0 or 1)
    # n_usans: find smoothing area (USAN) from secondary images (0, 1, or 2)

    arr = nb.load(file_in).get_fdata()
    arr_mean = np.mean(arr)
    arr_noise = np.mean(arr[:5, :5, :5])

    bt_threshold = 0.5 * (arr_mean - arr_noise)
    dim = 3
    use_median = 1
    n_usans = 0

    sf = 'susan'
    sf += " " + file_in
    sf += " " + str(bt_threshold)
    sf += " " + str(sigma)
    sf += " " + str(dim)
    sf += " " + str(use_median)
    sf += " " + str(n_usans)
    sf += " " + file_out

    dash = '-' * 40
    print(dash)
    print('{}'.format('Run susan noise filtering'))
    print(dash)
    print(sf)
    print("{:<22s}{:<.3f}".format("Brightness threshold: ", bt_threshold))
    print("{:<22s}{:<.1f}".format("Sigma: ", sigma))

    try:
        subprocess.check_output(sf, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        print(e.output)

    if os.path.isfile(file_out + ".gz"):
        gunzip(file_out + ".gz")


def _fill_with_gaussian(arr, arr_mask, sigma_gaussian=10.0):
    """Within a binary mask, all voxels are replaced by its gaussian filtered
    values. This is done to smoothed out edges and holes in the image array.

    Parameters
    ----------
    arr : (N,M,O) np.ndarray
        Image array.
    arr_mask : (N,M,O) np.ndarray
        Binary mask array.
    sigma_gaussian : float, optional
        Sigma for gaussian filter.

    Returns
    -------
    arr : (N,M,O) np.ndarray
        Filtered image array.

    """

    # apply gaussian filter to phase data
    arr_gaussian = gaussian_filter(arr, sigma_gaussian, order=0, output=None,
                                   mode='reflect', cval=0.0, truncate=4.0)

    # replace data
    arr[arr_mask == 1] = arr_gaussian[arr_mask == 1]

    return arr
