# -*- coding: utf-8 -*-

# python standard library inputs
import sys

# external inputs
import numpy as np
import matplotlib.pyplot as plt
from nibabel.affines import apply_affine

# local inputs
from gbb.interpolation import linear_interpolation3d
    
    
def get_shift(vtx, fac, vtx_norm, ind, arr_grad, vox2ras_tkr, ras2vox_tkr, 
              arr_vein=None, line_dir=2, line_length=3, t2s=True, 
              show_plot=True):
    """ Get shift

    This function computes the vertex shift in one direction (line_dir: 0,1,2) 
    or along the surface normal (line_dir: 3) towards the highest GM/WM 
    gradient. First, a line with a defined direction is computed for a chosen 
    vertex point and values of the second order gradient are sampled onto the 
    line. If the line direction is along one axis, lines are omitted if the 
    normal projection along that axis has a length below a set threshold. This 
    is done to exlude shifts for surface locations which are almost 
    perpendicularly oriented to the shift direction. Starting from the current 
    vertex point, the nearest zero crossing is found (if it exists). If a t2s 
    contrast is considered, the intensity values are going from dark to bright 
    for WM -> GM. In this case, we expect a transition from positive to negative 
    for the zero crossing in the second order gradient. Optionally, the found 
    shift is only considered if no vein is found in shift direction along the 
    line. Vertex coordinates are in ras space.    

    Parameters
    ----------
    vtx : ndarray
        Array of vertex points.
    fac : ndarray
        Array of corresponding faces.
    vtx_norm : ndarray
        Array of corresponding vertex normals.
    ind : int
        Current vertex index.
    arr_grad : ndarray
        3D array of second order gradient values along one axis.
    vox2ras_tkr : ndarray
        Transformation from voxel to ras space.
    ras2vox_tkr : ndarray
        Transformation from ras to voxel space.
    arr_vein : ndarray, optional
        3D array with masked veins. The default is None.
    line_dir : int, optional
        Line direction in ras conventions (0,1,2) or normal direction (3). The 
        default is 2.
    line_length : float, optional
        Length of vertex shift in one direction in mm. The default is 3.
    t2s : bool, optional
        WM darker than GM. The default is True.
    show_plot : bool, optional
        Show line plot in command window. The default is True.

    Returns
    -------
    shift_curr : ndarray
        Shift of vertex in ras coordinates.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 21-12-2019
    Last modified: 05-10-2020
    
    """
    
    # initialize
    shift_curr = []
    zero_found = []
    
    # fix parameters
    n_line = 1000 # number of line points
    line_threshold = 0.05 # if direction is along one axis, omit line if length is below threshold
    
    # maximum voxel coordinates in x-, y-, and z-direction
    vol_max = np.shape(arr_grad)
    
    # get current vertex and normal
    vtx_curr = vtx[ind,:].copy()
    normal_curr = vtx_norm[ind,:].copy()
    
    if line_dir == 0:
        normal_curr[0] = normal_curr[0] / np.abs(normal_curr[0])
        normal_curr[1] = 0
        normal_curr[2] = 0
        
        if np.abs(normal_curr[0]) < line_threshold:
            return []
        
    elif line_dir == 1:
        normal_curr[0] = 0
        normal_curr[1] = normal_curr[1] / np.abs(normal_curr[1])
        normal_curr[2] = 0
        
        if np.abs(normal_curr[1]) < line_threshold:
            return []
        
    elif line_dir == 2:
        normal_curr[0] = 0
        normal_curr[1] = 0
        normal_curr[2] = normal_curr[2] / np.abs(normal_curr[2])
        
        if np.abs(normal_curr[2]) < line_threshold:
            return []
        
    elif line_dir > 3 or line_dir < 0:
        sys.exit("error: choose a valid line direction!")
    
    # get line a -> b (WM -> GM)
    pt_start = vtx_curr.copy()
    pt_end = vtx_curr.copy()
    pt_start += line_length * normal_curr
    pt_end -= line_length * normal_curr
    
    line_curr = np.linspace((0,0,0),(1,1,1,),n_line,dtype=np.float)
    line_curr = (pt_end-pt_start) * line_curr + pt_start
    line_curr = apply_affine(ras2vox_tkr, line_curr)
    
    # remove coordinates exceeding the array limits
    outlier = np.zeros(len(line_curr), dtype=np.int8)
    outlier[line_curr[:,0] < 0] = 1
    outlier[line_curr[:,1] < 0] = 1
    outlier[line_curr[:,2] < 0] = 1
    outlier[line_curr[:,0] >= vol_max[0] - 1] = 1
    outlier[line_curr[:,1] >= vol_max[1] - 1] = 1
    outlier[line_curr[:,2] >= vol_max[2] - 1] = 1
    
    line_curr = line_curr[outlier == 0,:]
    n_line = len(line_curr) # update line length
    
    # get gradient and vein data along line
    grad_curr = linear_interpolation3d(line_curr[:,0],line_curr[:,1],line_curr[:,2],arr_grad)
    
    if arr_vein is not None:
        vein_curr = linear_interpolation3d(line_curr[:,0],line_curr[:,1],line_curr[:,2],arr_vein)
    
    # show line plot
    if show_plot:
        plt.clf()
        plt.plot(grad_curr)
        plt.xlabel("WM -> GM")
        plt.ylabel("Second order gradient")
    
    # mid-point of line
    i = np.floor(n_line / 2).astype(int)
    
    # get point of zero crossing (closest to mid-point) if grad_curr exists, i.e., if it does not 
    # contain only outliers
    j = 0
    switch = 0
    if len(grad_curr):
    
        # look for veins in line
        if arr_vein is not None:
            vein_up = np.sum(vein_curr[i:])
            vein_down = np.sum(vein_curr[:i])
        else:
            vein_up = 0
            vein_down = 0
    
        # start search
        while i > 0 and i < n_line - 2:
            
            # if vein up -> i nur in negativer richtung
            # if vein down -> i nur in positiver richtung
            # if not -> look in both directions (switch)
            # if both -> return empty array
            if vein_up and not vein_down:
                i -= 1
            elif not vein_up and vein_down:
                i += 1
            elif not vein_up and not vein_down:
                j += 1
                if switch == 0:
                    switch = 1
                    i += j
                else:
                    switch = 0
                    i -= j
            else:
                break
            
            if grad_curr[i] < 0 and grad_curr[i+1] > 0 and t2s:
                zero_found = 1
                break
            elif grad_curr[i] > 0 and grad_curr[i+1] < 0 and not t2s:
                zero_found = 1
                break
    
    # only consider shift if not positioned within a vein
    if zero_found:
        zero_curr = apply_affine(vox2ras_tkr, line_curr[i,:])
        shift_curr = vtx_curr - zero_curr
    
    return shift_curr
