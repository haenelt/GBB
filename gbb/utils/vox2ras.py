# -*- coding: utf-8 -*-

# python standard library inputs
import subprocess

# external inputs
import numpy as np
import nibabel as nb
from numpy.linalg import inv


def vox2ras(file_in):
    """ VOX to RAS

    This function reads an input volume and computes the transformation between 
    voxel space and freesurfer vertex RAS (right-anterior-superior) coordinate 
    system from the header information. Transformation for both directions are 
    returned.    

    Parameters
    ----------
    file_in : str
        Filename of nifti image.

    Returns
    -------
    vox2ras_tkr : ndarray
        Transformation matrix from voxel to ras space.
    ras2vox_tkr : ndarray
        Transformation matrix from ras to voxel space.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 18-12-2019
    Last modified: 05-10-2020
    
    """
    
    # get affine vox2ras-tkr and ras2vox-tkr transformation to reference volume
    freesurfer_command = ['mri_info', file_in, '--{}'.format("ras2vox-tkr")]
    transformation = subprocess.check_output(freesurfer_command, 
                                             stderr=subprocess.DEVNULL).decode()
    
    # ignore if warning is stated in first line
    if transformation[:7] == "WARNING":
        i = 0
        while True:
            if transformation[i] == "\n":
                transformation = transformation[i+1:]
                break
            else:
                i += 1
    
    num_transformation = [[float(x) for x in line.split()] 
                          for line in transformation.split('\n') if len(line)>0]
    
    # get final transformation matriced as numpy array
    ras2vox_tkr = np.array(num_transformation)
    vox2ras_tkr = inv(np.array(num_transformation))
    
    return vox2ras_tkr, ras2vox_tkr


def vox2ras2(file_in):
    """
    https://neurostars.org/t/get-voxel-to-ras-transformation-from-nifti-file/4549/3

    Parameters
    ----------
    file_in

    Returns
    -------

    """

    nii = nb.load(file_in)
    mgh = nb.MGHImage(nii.dataobj, nii.affine)
    vox2ras_tkr = mgh.header.get_vox2ras_tkr()

    # get final transformation matriced as numpy array
    ras2vox_tkr = inv(np.array(vox2ras_tkr))


    return vox2ras_tkr, ras2vox_tkr
