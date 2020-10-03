# -*- coding: utf-8 -*-
import numpy as np
from scipy.sparse import csr_matrix
from nibabel.freesurfer.io import read_geometry
    
    
def get_adjm(input_surf):
    """
    This function computes a sparse adjacency matrix from a surface mesh in 
    freesurfer format. The matrix has the size (nvertex,nvertex). Each matrix 
    entry with value 1 stands for an edge of the surface mesh.
    Inputs:
        *input_surf (str): filename of surface mesh.
    Outputs:
        *sparse_adjm (obj): sparse adjacency matrix.
        
    created by Daniel Haenelt
    Date created: 20-12-2019
    Last modified: 03-10-2020
    """

    # load surface mesh
    vtx, fac = read_geometry(input_surf)

    # get number of vertices and faces
    nvertex = len(vtx[:,0])
    nface   = len(fac[:,0])
    
    # initialise
    row = []
    col = []
    
    # get rows and columns of edges
    row.extend( [fac[i,0] for i in range(nface)] )
    col.extend( [fac[i,1] for i in range(nface)] )
    
    row.extend( [fac[i,1] for i in range(nface)] )
    col.extend( [fac[i,2] for i in range(nface)] )
    
    row.extend( [fac[i,2] for i in range(nface)] )
    col.extend( [fac[i,0] for i in range(nface)] )
        
    # make sure that all edges are symmetric
    row.extend( [fac[i,1] for i in range(nface)] )
    col.extend( [fac[i,0] for i in range(nface)] )
    
    row.extend( [fac[i,2] for i in range(nface)] )
    col.extend( [fac[i,1] for i in range(nface)] )
    
    row.extend( [fac[i,0] for i in range(nface)] )
    col.extend( [fac[i,2] for i in range(nface)] )
    
    # adjacency entries get value 1
    data = np.ones(len(row), dtype=np.int8)
    
    # write sparse adjacency matrix
    sparse_adjm = csr_matrix((data, (row, col)), shape=(nvertex,nvertex))
    
    return sparse_adjm
