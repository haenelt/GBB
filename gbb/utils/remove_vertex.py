# -*- coding: utf-8 -*-

# python standard library inputs
import itertools

# external inputs
import numpy as np


def remove_vertex(vtx, fac, ind_keep):
    """ Remove vertex
    
    This function removes vertices from an array of vertex points and updates 
    the corresponding face array.    

    Parameters
    ----------
    vtx : ndarray
        Array of vertices.
    fac : ndarray
        Array of faces.
    ind_keep : list
        Index list of vertices to keep.

    Returns
    -------
    vtx : ndarray
        Remaining vertices.
    fac : ndarray
        Remaining faces.
    ind_keep : ndarray
        Updated index list of vertices to keep after vertex cleaning.
    
    Notes
    -------
    created by Daniel Haenelt
    Date created: 22-06-2020
    Last modified: 19-10-2020

    """
    
    # get indices which will be removed
    ind_tmp = np.arange(len(vtx))
    ind_remove = list(set(ind_tmp) - set(ind_keep))
    ind_remove = sorted(ind_remove, reverse=True)

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

    # get indices which will be cleaned
    ind_vtx = np.arange(len(vtx))
    ind_fac = list(itertools.chain(*fac))
    ind_fac = list(set(ind_fac))
    ind_remove = list(set(ind_vtx) - set(ind_fac))
    ind_remove = sorted(ind_remove, reverse=True)

    # remove singularities (vertices without faces)
    loop_status = 0
    loop_length = len(ind_remove)
    for i in range(loop_length):
        
        # print status
        counter = np.floor(i / loop_length * 100)
        if counter != loop_status:
            print("clean faces: "+str(counter)+" %")
            loop_status = counter
                   
        # remove vertex and index
        vtx = np.delete(vtx, ind_remove[i], 0)
        ind_keep = np.delete(ind_keep, ind_remove[i], 0)    
            
        # sort faces
        tmp = fac[fac >= ind_remove[i]] - 1
        fac[fac >= ind_remove[i]] = tmp
    
    return vtx, fac, ind_keep
