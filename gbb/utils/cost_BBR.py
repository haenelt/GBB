# -*- coding: utf-8 -*-
import numpy as np
from nibabel.affines import apply_affine
from gbb.interpolation import linear_interpolation3d
    

def cost_BBR(vtx, vtx_n, arr_vol, ras2vox, vol_max, Q0=0, M=0.5, h=1, t2s=True):
    """
    This function computes the cost function defined in the original BBR paper 
    (Greve and Fischl, 2009). Vertices should be in freesurfer ras_tkr space. 
    First, vertex coordinates on both sides normal to the surface in 1 mm 
    distance are computed. Based on a tranformation matrix, vertex coordinates 
    are transformed then into the voxel space of the input array. Vertices are 
    excluded with coordinates exceeding the array limits. GM and WM values are 
    sampled using linear interpolation. A percent contrast measure and the cost 
    function value are then calculated.
    Inputs:
        *vtx (arr): array of vertices.
        *vtx_n (arr): array of unit vertices in normal direction.
        *vol_array (arr): 3D array of image volume.
        *ras2vox (arr): transformation matrix to voxel space.
        *vol_max (arr): array of maximum voxel coordinates in x-, y-, and z-direction.
        *Q0 (float): offset parameter in percent contrast measure.
        *M (float): slope parameter in percent contrast measure.
        *h (float): weight for each vertex in percent contrast measure.
        *t2s (bool): t2star weighing, i.e., gm > wm intensity.
    Outputs:
        *J (float): cost function value.
        
    created by Daniel Haenelt
    Date created: 21-12-2019
    Last modified: 03-10-2020
    """
       
    # sort offset in two groups according to normal direction
    gm_pts = vtx + vtx_n
    wm_pts = vtx - vtx_n
    
    # ras2vox transformation
    gm_pts = apply_affine(ras2vox, gm_pts)
    wm_pts = apply_affine(ras2vox, wm_pts)
    
    # get location of outlier coordinates
    outlier = np.zeros(len(gm_pts), dtype=np.int8())
    outlier[gm_pts[:,0] > vol_max[0] - 1] = 1
    outlier[gm_pts[:,1] > vol_max[1] - 1] = 1
    outlier[gm_pts[:,2] > vol_max[2] - 1] = 1
    outlier[wm_pts[:,0] > vol_max[0] - 1] = 1
    outlier[wm_pts[:,1] > vol_max[1] - 1] = 1
    outlier[wm_pts[:,2] > vol_max[2] - 1] = 1
    outlier[gm_pts[:,0] < 0] = 1
    outlier[gm_pts[:,1] < 0] = 1
    outlier[gm_pts[:,2] < 0] = 1
    outlier[wm_pts[:,0] < 0] = 1
    outlier[wm_pts[:,1] < 0] = 1
    outlier[wm_pts[:,2] < 0] = 1
    
    # get rid of outliers
    gm_pts = gm_pts[outlier == 0]
    wm_pts = wm_pts[outlier == 0]
    
    # get values in GM and WM
    gm_val = linear_interpolation3d(gm_pts[:,0], gm_pts[:,1], gm_pts[:,2], arr_vol)
    wm_val = linear_interpolation3d(wm_pts[:,0], wm_pts[:,1], wm_pts[:,2], arr_vol)
    
    # percent contrast measure
    if t2s == True:
        Q = 100 * ( wm_val - gm_val ) / ( 0.5 * ( gm_val + wm_val ) )
    else:
        Q = 100 * ( gm_val - wm_val ) / ( 0.5 * ( gm_val + wm_val ) )
    
    # cost value
    J = 1 / len(Q) * np.sum( h * ( 1 + np.tanh(M * ( Q - Q0 )) ) )
    
    return J
