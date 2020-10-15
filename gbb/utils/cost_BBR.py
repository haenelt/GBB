# -*- coding: utf-8 -*-

# external inputs
import numpy as np
from nibabel.affines import apply_affine

# local inputs
from gbb.interpolation import linear_interpolation3d
    

def cost_BBR(vtx, vtx_n, arr_vol, ras2vox_tkr, Q0=0, M=0.5, h=1, s=1, t2s=True):
    """ Cost BBR

    This function computes the cost function defined in the original BBR paper 
    [1]. Vertices should be in freesurfer ras_tkr space. First, vertex 
    coordinates on both sides normal to the surface are computed. Based on a 
    tranformation matrix, vertex coordinates are transformed then into voxel 
    space of the input array. Vertices are excluded with coordinates exceeding 
    the array limits. GM and WM values are sampled using linear interpolation. A 
    percent contrast measure and the cost function value are calculated.    

    Parameters
    ----------
    vtx : ndarray
        Array of vertices.
    vtx_n : ndarray
        Array of unit vertices in normal direction.
    arr_vol : ndarray
        3D array of image volume.
    ras2vox_tkr : ndarray
        Transformation from ras to voxel space.
    Q0 : float, optional
        Offset parameter in percent contrast measure. The default is 0.
    M : float, optional
        Slope parameter in percent contrast measure. The default is 0.5.
    h : float, optional
        Weight for each vertex in percent contrast measure. The default is 1.
    s : float, optional
        Distance scaling factor for sampling normal to the surface. The default 
        is 1.
    t2s : bool, optional
        T2star weighing, i.e., gm > wm intensity. The default is True.

    Returns
    -------
    J : float
        Cost function value.

    References
    -------
    .. [1] Greve, DN, Fischl, B, Accurate and robust brain image alignment using 
    boundary-based registration, Neuroimage 48(1), 63--72 (2009).

    Notes
    -------
    created by Daniel Haenelt
    Date created: 21-12-2019
    Last modified: 08-10-2020
    
    """
    
    # maximum voxel coordinates in x-, y-, and z-direction
    vol_max = np.shape(arr_vol)
       
    # sort offset in two groups according to normal direction
    gm_pts = vtx + s * vtx_n
    wm_pts = vtx - s * vtx_n
    
    # ras2vox transformation
    gm_pts = apply_affine(ras2vox_tkr, gm_pts)
    wm_pts = apply_affine(ras2vox_tkr, wm_pts)
    
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
