# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
from nibabel.freesurfer.io import write_geometry
from nibabel.affines import apply_affine

# local inputs
from gbb.interpolation.linear_interpolation3d import linear_interpolation3d
from gbb.utils.remove_vertex import remove_vertex
    
    
def apply_deformation(vtx, fac, arr_deform, vox2ras_tkr, ras2vox_tkr, 
                      path_output="", name_output="", write_output=True):
    """ Apply deformation

    This function applies a coordinate shift to an array of vertices. The 
    coordinate shift is taken from a deformation field where each voxel 
    corresponds to a shift along one direction in voxel space. Vertices outside 
    the deformation field are removed from the mesh.    

    Parameters
    ----------
    vtx : ndarray
        Array of vertices.
    fac : ndarray
        Array of faces.
    arr_deform : ndarray
        4D volume array containing shifts in voxel space.
    vox2ras_tkr : ndarray
        Transformation from voxel to ras space.
    ras2vox_tkr : ndarray
        Transformation from ras to voxel space.
    path_output : str, optional
        Path where output is written. The default is "".
    name_output : str, optional
        Name of output file. The default is "".
    write_output : bool, optional
        Write output file. The default is True.

    Returns
    -------
    vtx : ndarray
        Deformed vertices.
    fac : ndarray
        Corresponding faces.
    ind_keep : ndarray
        Input vertex indices to keep after cleaning.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 28-12-2019
    Last modified: 23-10-2020
    
    """
       
    # get array dimensions
    xdim = np.shape(arr_deform)[0]
    ydim = np.shape(arr_deform)[1]
    zdim = np.shape(arr_deform)[2]    
    
    # get vertices to voxel space
    vtx = apply_affine(ras2vox_tkr, vtx)
    
    # only keep vertex indices within the slab
    vtx[vtx[:,0] < 0,0] = np.nan
    vtx[vtx[:,1] < 0,1] = np.nan
    vtx[vtx[:,2] < 0,2] = np.nan
    
    vtx[vtx[:,0] > xdim - 1,0] = np.nan
    vtx[vtx[:,1] > ydim - 1,1] = np.nan
    vtx[vtx[:,2] > zdim - 1,2] = np.nan
    
    # remove vertices outside the slab
    if np.any(np.isnan(vtx)):
        ind_keep = np.arange(len(vtx))
        ind_keep = ind_keep[~np.isnan(np.sum(vtx, axis=1))]
        
        vtx, fac, ind_keep = remove_vertex(vtx, fac, ind_keep)
    
    # sample shifts from deformation map
    vtx_shift = np.zeros((len(vtx),3))
    for i in range(3):
        vtx_shift[:,i] = linear_interpolation3d(vtx[:,0], vtx[:,1], vtx[:,2], 
                                                arr_deform[:,:,:,i])
        vtx_shift[np.isnan(vtx_shift[:,i]),i] = 0
    
    # shift coordinates to new location
    vtx += vtx_shift
    
    # get new coordinates in ras space
    vtx = apply_affine(vox2ras_tkr, vtx)
    
    if write_output:
        file_out = os.path.join(path_output, name_output)
        write_geometry(file_out, vtx, fac)
        np.savetxt(file_out+"_ind.txt", ind_keep, fmt="%d")
    
    return vtx, fac, ind_keep
