# -*- coding: utf-8 -*-

# external inputs
import numpy as np


def nn_3d(vtx0, vtx, r_size):
    """
    This function computes all nearest neighbors found within a sphere with a 
    radius defined in ras coordinates. Note that the defined neighborhood does 
    not have to be fully connected in this case. Vertex coordinates are in ras 
    space.
    Inputs:
        *vtx0 (arr): vertex point.
        *vtx (arr): array of vertices.
        *r_size (float): radius of sphere in ras coordinates.
    Outputs:
        *nn (arr): array of neighbor indices.
        *r[nn] (arr): euclidean distance to neighbors.
        
    created by Daniel Haenelt
    Date created: 22-12-2019     
    Last modified: 05-10-2020
    """
    
    rx = ( vtx[:,0] - vtx0[0] ) ** 2
    ry = ( vtx[:,1] - vtx0[1] ) ** 2
    rz = ( vtx[:,2] - vtx0[2] ) ** 2
    
    r = np.sqrt( rx + ry + rz )
    
    nn = np.where(r < r_size)[0]
    
    return nn, r[nn]
