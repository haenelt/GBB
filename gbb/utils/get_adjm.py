# -*- coding: utf-8 -*-

# external inputs
import numpy as np
from scipy.sparse import csr_matrix


def get_adjm(vtx, fac):
    """
    This function computes a sparse adjacency matrix for a triangular surface
    mesh. The matrix has the size (nvertex,nvertex). Each matrix entry with 
    value 1 stands for an edge of the surface mesh.
    Inputs:
        *vtx (arr): array of vertices.
        *fac (arr): array of faces.
    Outputs:
        *sparse_adjm (obj): sparse adjacency matrix.
        
    created by Daniel Haenelt
    Date created: 20-12-2019
    Last modified: 09-10-2020
    """

    # get number of vertices and faces
    nvtx = len(vtx)
    nfac   = len(fac)
    
    # initialise
    row = []
    col = []
    
    # get rows and columns of edges
    row.extend( [fac[i,0] for i in range(nfac)] )
    col.extend( [fac[i,1] for i in range(nfac)] )
    
    row.extend( [fac[i,1] for i in range(nfac)] )
    col.extend( [fac[i,2] for i in range(nfac)] )
    
    row.extend( [fac[i,2] for i in range(nfac)] )
    col.extend( [fac[i,0] for i in range(nfac)] )
        
    # make sure that all edges are symmetric
    row.extend( [fac[i,1] for i in range(nfac)] )
    col.extend( [fac[i,0] for i in range(nfac)] )
    
    row.extend( [fac[i,2] for i in range(nfac)] )
    col.extend( [fac[i,1] for i in range(nfac)] )
    
    row.extend( [fac[i,0] for i in range(nfac)] )
    col.extend( [fac[i,2] for i in range(nfac)] )
    
    # adjacency entries get value 1
    data = np.ones(len(row), dtype=np.int8)
    
    # write sparse adjacency matrix
    sparse_adjm = csr_matrix((data, (row, col)), shape=(nvtx,nvtx))
    
    return sparse_adjm
