# -*- coding: utf-8 -*-

# external inputs
import numpy as np

# local inputs
from gbb.utils.get_adjm import get_adjm
from gbb.neighbor.nn_2d import nn_2d


def remove_vertex(vtx, fac, ind_keep):
    """
    This function removes vertices from an array of vertex points and updates 
    the corresponding face array.
    Inputs:
        *vtx (arr): array of vertices.
        *fac (arr): array of faces.
        *ind_keep (list): index list of vertices to keep.
    Outputs:
        *vtx (arr): remaining vertices.
        *fac (arr): remaining faces.
        
    created by Daniel Haenelt
    Date created: 22-06-2020
    Last modified: 08-10-2020
    """

    # get indices which will be removed
    ind_tmp = np.arange(len(vtx))
    ind_remove = list(set(ind_tmp) - set(ind_keep))

    # get new vertices
    vtx = vtx[ind_keep,:]

    # get new faces
    fac_keep = np.zeros(len(fac))
    fac_keep += np.in1d(fac[:,0], ind_keep)
    fac_keep += np.in1d(fac[:,1], ind_keep)
    fac_keep += np.in1d(fac[:,2], ind_keep)
    fac = fac[fac_keep == 3,:]
    
    # reindex faces
    loop_status = 0
    loop_length = len(ind_remove)
    for i in range(loop_length):
        
        # print status
        counter = np.floor(i / loop_length * 100)
        if counter != loop_status:
            print("sort faces: "+str(counter)+" %")
            loop_status = counter
        
        tmp = fac[fac >= ind_remove[i]] - 1
        fac[fac >= ind_remove[i]] = tmp

    # get adjacency matrix
    adjm = get_adjm(vtx, fac)

    # remove singularities (vertices without faces)
    loop_status = 0
    loop_length = len(vtx)
    for i in range(loop_length):
        
        # print status
        counter = np.floor(i / loop_length * 100)
        if counter != loop_status:
            print("clean faces: "+str(counter)+" %")
            loop_status = counter
        
        n_neighbor = nn_2d(i, adjm, 0)
        if not len(n_neighbor):
            
            # remove vertex
            vtx = np.delete(vtx, i, 0)
            
            # sort faces
            tmp = fac[fac >= i] - 1
            fac[fac >= i] = tmp
    
    return vtx, fac
