# -*- coding: utf-8 -*-

# external inputs
import numpy as np


def update_mesh(vtx, vtx_shift, source_ind, nn_ind, l_rate=0.9):
    """
    This function updates the vertices based on a calculated vertex shift within 
    a defined neighboorhood and specified learning rate. If the neighborhood 
    only contains the source vertex, no distance weighting is computed. The 
    learning rate should be a value between 0 and 1. Vertex coordinates are in 
    ras space.
    Inputs:
        *vtx (arr): array of vertices.
        *vtx_shift (arr): vertex shift in mm.
        *source_ind (int): source index.
        *nn_ind (list): neighborhood indices.
        *l_rate (float): learning rate.
    Outputs:
        *vtx_new (arr): updated array of vertices.
        
    created by Daniel Haenelt
    Date created: 22-12-2019      
    Last modified: 05-10-2020
    """
    
    if set(nn_ind) == set([source_ind]):
        s = [1]
    else:
        # get distance to source vertex
        sx = ( vtx[nn_ind,0] - vtx[source_ind,0] ) ** 2
        sy = ( vtx[nn_ind,1] - vtx[source_ind,1] ) ** 2
        sz = ( vtx[nn_ind,2] - vtx[source_ind,2] ) ** 2

        # scale shift depending on learning rate and neighborhood distance
        s = np.sqrt( sx + sy + sz )
        s = ( np.max(s) - s ) / np.max(s)

    vtx_new = vtx.copy()
    for i in range(len(nn_ind)):
        vtx_new[nn_ind[i],:] = vtx_new[nn_ind[i],:] - l_rate * s[i] * vtx_shift
    
    return vtx_new
