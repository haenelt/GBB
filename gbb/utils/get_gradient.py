# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import cv2
import numpy as np
import nibabel as nb
from scipy.ndimage.filters import gaussian_filter

# local inputs
from gbb.utils.line_ras2vox import line_ras2vox


def get_gradient(input_vol, ras2vox_tkr, line_dir=2, sigma=0, kernel_size=3, 
                 write_output=None, path_output="", name_output="gradient"):
    """ Get gradient
    
    This function computes the second order gradient along one axis 
    (line_dir: 0,1,2) or the sum over all directions (line_dir: 3).    

    Parameters
    ----------
    input_vol : str
        Filename of input nifti.
    ras2vox_tkr : ndarray
        Transformation from ras to voxel space.
    line_dir : int, optional
        Direction of gradient calculation in ras space. The default is 2.
    sigma : float, optional
        Gaussian blurring kernel of gradient map (if set > 0). The default is 0.
    kernel_size : float, optional
        Kernel size for gradient calculation. The default is 3.
    write_output : bool, optional
        Write output image. The default is None.
    path_output : str, optional
        Path where output is written. The default is "".
    name_output : str, optional
        Basename of output file. The default is "gradient".

    Returns
    -------
    res : ndarray
        Gradient array.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 30-10-2019
    Last modified: 08-10-2020
    
    """
    
    # get line direction in voxel space
    if line_dir != 3:
        line_dir = line_ras2vox(line_dir, ras2vox_tkr)

    # load input
    vol = nb.load(input_vol)
    vol_array = vol.get_fdata()

    # get gradient
    res = np.zeros_like(vol_array)
    if line_dir == 0:
        for i in range(np.shape(vol_array)[2]):
            res[:,:,i] = cv2.Sobel(vol_array[:,:,i],cv2.CV_64F,0,1,kernel_size)
            res[:,:,i] = cv2.Sobel(res[:,:,i],cv2.CV_64F,0,1,kernel_size)
    
    elif line_dir == 1:
        for i in range(np.shape(vol_array)[2]):
            res[:,:,i] = cv2.Sobel(vol_array[:,:,i],cv2.CV_64F,1,0,kernel_size)
            res[:,:,i] = cv2.Sobel(res[:,:,i],cv2.CV_64F,1,0,kernel_size)        
    
    elif line_dir == 2:
        for i in range(np.shape(vol_array)[0]):
            res[i,:,:] = cv2.Sobel(vol_array[i,:,:],cv2.CV_64F,1,0,kernel_size)
            res[i,:,:] = cv2.Sobel(res[i,:,:],cv2.CV_64F,1,0,kernel_size) 
    
    elif line_dir == 3:
        res_x = np.zeros_like(res)
        for i in range(np.shape(vol_array)[2]):
            res_x[:,:,i] = cv2.Sobel(vol_array[:,:,i],cv2.CV_64F,0,1,kernel_size)
            res_x[:,:,i] = cv2.Sobel(res_x[:,:,i],cv2.CV_64F,0,1,kernel_size)
        
        res_y = np.zeros_like(res)
        for i in range(np.shape(vol_array)[2]):
            res_y[:,:,i] = cv2.Sobel(vol_array[:,:,i],cv2.CV_64F,1,0,kernel_size)
            res_y[:,:,i] = cv2.Sobel(res_y[:,:,i],cv2.CV_64F,1,0,kernel_size)   
        
        res_z = np.zeros_like(res)
        for i in range(np.shape(vol_array)[0]):
            res_z[i,:,:] = cv2.Sobel(vol_array[i,:,:],cv2.CV_64F,1,0,kernel_size)
            res_z[i,:,:] = cv2.Sobel(res_z[i,:,:],cv2.CV_64F,1,0,kernel_size) 

        res = res_x + res_y + res_z

    # gaussian blurring (optional)
    if sigma:
        res = gaussian_filter(res, 
                              sigma, 
                              order = 0, 
                              output = None, 
                              mode = 'reflect', 
                              cval = 0.0, 
                              truncate = 4.0)
    
    # write gradient image
    if write_output:
        output = nb.Nifti1Image(res, vol.affine, vol.header)
        nb.save(output,os.path.join(path_output,name_output+".nii"))
    
    return res
