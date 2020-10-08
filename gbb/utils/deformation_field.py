# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb
from nibabel.affines import apply_affine
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter

    
def deformation_field(vtx_old, vtx_new, input_vol, ras2vox_tkr, sigma=1, 
                      path_output="", name_output="", write_output=True):
    """
    This function computes a deformation field from an array of shifted 
    vertices. Each voxel in the deformation field corresponds to a shift in 
    voxel space along one direction. Optionally, a gaussian filter can be 
    applied to the final deformation field.
    Inputs:
        *vtx_old (arr): original array of vertices.
        *vtx_new (arr): new array of vertices.
        *input_vol (str): filename of reference volume.
        *ras2vox_tkr (arr): ras to voxel space transformation.
        *sigma (float): sigma for gaussian filtering of deformation field.
        *path_output (str): path where output is written.
        *name_output (str): basename of output file without file extension.
        *write_output (bool): write output file.
    Outputs:
        *arr_deform (arr): deformation field.
        
    created by Daniel Haenelt
    Date created: 28-12-2019       
    Last modified: 08-10-2020
    """
    
    # load reference volume
    vol = nb.load(input_vol)
    vol.header["dim"][4] = 3
    
    # get volume dimension
    x_dim = vol.header["dim"][1]
    y_dim = vol.header["dim"][2]
    z_dim = vol.header["dim"][3]
    
    # initialize deformation array
    arr_deform = np.zeros(vol.header["dim"][1:5])
    
    # get vertices to voxel space
    vtx_old = apply_affine(ras2vox_tkr, vtx_old)
    vtx_new = apply_affine(ras2vox_tkr, vtx_new)
    
    # get vertex shifts in voxel space
    vtx_shift = vtx_new - vtx_old
    
    # get volume coordinates in voxel space
    xf = np.arange(0,x_dim)
    yf = np.arange(0,y_dim)
    zf = np.arange(0,z_dim)
    
    y_plane, x_plane, z_plane = np.meshgrid(yf,xf,zf)
    x_plane = x_plane.flatten()
    y_plane = y_plane.flatten()
    z_plane = z_plane.flatten()
    
    # grid interpolation
    nn = griddata(vtx_old, vtx_shift, np.column_stack([x_plane,y_plane,z_plane]), method='nearest')
    li = griddata(vtx_old, vtx_shift, np.column_stack([x_plane,y_plane,z_plane]), method='linear')
    
    # fill all nans from the linear interpolation with nearest neighbor values
    li[np.isnan(li)] = nn[np.isnan(li)]
    
    # reshape to deformation volume
    for i in range(3):
        
        # reshape
        arr_deform[:,:,:,i] = np.reshape(li[:,i],(x_dim,y_dim,z_dim))
    
        # apply gaussian filter
        if sigma:
            arr_deform[:,:,:,i] = gaussian_filter(arr_deform[:,:,:,i], sigma)
        
    # write output    
    if write_output:
        output = nb.Nifti1Image(arr_deform, vol.affine, vol.header)
        nb.save(output, os.path.join(path_output,name_output+".nii"))
    
    return arr_deform
