# -*- coding: utf-8 -*-

# python standard library inputs
import os
import sys

# external inputs
import numpy as np
from numpy.linalg import norm
from nibabel.freesurfer.io import read_geometry, write_geometry

# local inputs
from gbb.neighbor.nn_2d import nn_2d
from gbb.utils.update_mesh import update_mesh
from gbb.utils.smooth_surface import smooth_surface


def run_anchor(vtx, fac, adjm, anchor, n_neighbor=20, smooth_iter=0):
    """
    This function loads a list of control points and shifts a mesh with its 
    local neighborhood to match these control points. Control points are assumed 
    to be in ras coordinates. Optionally, the output mesh can be smoothed.    
    Inputs.
        *vtx (arr): array of vertex points.
        *fac (arr): array of corresponding faces.
        *anchor (list): list of control points.
        *n_neighbor (int): neighborhood size.
        *smooth_iter (int): number of smoothing iterations of final mesh.
    Outputs:
        *vtx (arr): shifted array of vertex points.
        *ind_control (arr): indices of closest vertices to control points.

    created by Daniel Haenelt
    Date created: 12-05-2020 
    Last modified: 05-10-2020  
    """
    
    print("start mesh initialization (anchoring)")
    
    # check vein array
    if anchor is None:
        sys.exit("error: no control points found for anchoring!")
    
    # number of anchor points
    n_anchor = len(anchor)
    
    # loop through control points
    ind_control = np.zeros(n_anchor).astype(int)
    for i in range(n_anchor):
    
        # print current status
        print("apply control point: "+str(i+1)+"/"+str(n_anchor))
    
        # get nearest vertex
        vtx_dist = norm(vtx-anchor[i,:], axis=1)
        ind_min = np.where(vtx_dist == np.min(vtx_dist))[0][0]
        ind_control[i] = ind_min
        
        # get neighborhood
        nn_ind = nn_2d(ind_min, adjm, n_neighbor)
        
        # get shift
        vtx_shift = vtx[ind_min,:] - anchor[i,:]
        
        # update mesh
        vtx = update_mesh(vtx, vtx_shift, ind_min, nn_ind, 1)
    
    # smooth output
    if smooth_iter:
        tmp = np.random.randint(0, 10, 5)
        tmp_string = ''.join(str(i) for i in tmp)
        surf_temp = os.path.join(os.getcwd(),"surf_"+tmp_string)
        write_geometry(surf_temp, vtx, fac)
        smooth_surface(surf_temp, surf_temp, smooth_iter)
        vtx, _ = read_geometry(surf_temp)
        os.remove(surf_temp)
        
    return vtx, ind_control
