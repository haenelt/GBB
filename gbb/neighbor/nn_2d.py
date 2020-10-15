# -*- coding: utf-8 -*-

# external inputs
import numpy as np


def nn_2d(ind, adjm, n_iter):
    """ NN 2D

    This function finds all nearest neigbors on the surface mesh which show a 
    connection within the number of iterations steps.    

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
    nn : ndarray
        Array of neighbor indices.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 21-12-2019     
    Last modified: 05-10-2020

    """
    
    # get first order neighbours
    nn = adjm[ind,:].indices

    i = 0
    while i < n_iter - 1:
        nn_temp = nn.copy()
        for j in range(len(nn_temp)):
            nn = np.append(nn, adjm[nn_temp[j],:].indices)
            nn = np.unique(nn)

        i += 1
        
    return nn
