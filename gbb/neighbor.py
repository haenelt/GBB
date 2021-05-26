# -*- coding: utf-8 -*-
"""Functions for finding nearest vertex neighbors."""

# external inputs
import numpy as np

__all__ = ['nn_2d', 'nn_3d']


def nn_2d(ind, adjm, n_iter):
    """Neighbor 2D.

    This function finds all nearest neighbors of a source vertex. By applying 
    the function iteratively, the output is used as source for the next 
    iteration to grow the neighborhood.

    Parameters
    ----------
    ind : int
        Vertex index.
    adjm : obj
        Sparse adjacency matrix (in csr_matrix format).
    n_iter : int
        Number of iterations.

    Returns
    -------
    nn : (N,) np.ndarray
        Array of neighbor indices.

    """
    
    # get first order neighbours
    nn = adjm[ind, :].indices

    i = 0
    while i < n_iter - 1:
        nn_temp = nn.copy()
        for j in range(len(nn_temp)):
            nn = np.append(nn, adjm[nn_temp[j],:].indices)
            nn = np.unique(nn)

        i += 1
        
    return nn


def nn_3d(vtx0, vtx, r_size):
    """Neighbor 3D.

    This function finds all vertex neighbors within a sphere. The radius is 
    defined in tkr-ras coordinates. Note that the defined neighborhood does 
    not have to be fully connected in this case. Vertex coordinates are in 
    tkr-ras space.    

    Parameters
    ----------
    vtx0 : (3,) np.ndarray
        Center coordinate.
    vtx : (N,) np.ndarray
        Vertex array.
    r_size : float
        Radius of sphere in tkr-ras coordinates.

    Returns
    -------
    nn : ndarray
        Array of neighbor indices.
    
    """
    
    rx = ( vtx[:,0] - vtx0[0] ) ** 2
    ry = ( vtx[:,1] - vtx0[1] ) ** 2
    rz = ( vtx[:,2] - vtx0[2] ) ** 2
    
    r = np.sqrt( rx + ry + rz )
    
    nn = np.where(r < r_size)[0]
    
    return nn
