# -*- coding: utf-8 -*-

# external inputs
import numpy as np


def nn_3d(vtx0, vtx, r_size):
    """ NN 3D

    This function computes all nearest neighbors found within a sphere with a 
    radius defined in ras coordinates. Note that the defined neighborhood does 
    not have to be fully connected in this case. Vertex coordinates are in ras 
    space.    

    Parameters
    ----------
    vtx0 : ndarray
        Vertex point.
    vtx : ndarray
        Array of vertices.
    r_size : float
        Radius of sphere in ras coordinates.

    Returns
    -------
    nn : ndarray
        Array of neighbor indices.
    r[nn] : ndarray
        Euclidean distance to neighbors.

    Notes
    -------
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
