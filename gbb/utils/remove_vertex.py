# -*- coding: utf-8 -*-

# external inputs
import numpy as np


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
    Last modified: 05-10-2020
    """

    # get new vertices
    vtx = vtx[ind_keep,:]

    # get new faces
    fac_keep = np.zeros(len(fac[:,0]))
    fac_keep += np.in1d(fac[:,0],ind_keep)
    fac_keep += np.in1d(fac[:,1],ind_keep)
    fac_keep += np.in1d(fac[:,2],ind_keep)
    fac_temp = fac[fac_keep == 3,:]
    fac = fac[fac_keep == 3,:]

    # sort new faces
    c_step = 0
    n_step = [10,20,30,40,50,60,70,80,90,100]
    for i in range(len(ind_keep)):
        temp = np.where(ind_keep[i] == fac_temp)
        fac[temp] = i
    
        # print status
        counter = np.floor(i / len(ind_keep) * 100).astype(int)
        if counter == n_step[c_step]:
            print("sort faces: "+str(counter)+" %")
            c_step += 1
    
    # remove singularities (vertices without faces)
    fac_counter = 0
    fac_old = fac.copy()
    n_singularity = np.zeros(len(vtx))
    c_step = 0
    for i in range(len(vtx)):
        row, col = np.where(fac_old == i)
 
        n_singularity[i] = len(row)
        if not n_singularity[i]:    
            fac_temp = fac.copy()
            fac_temp[fac_temp >= fac_counter] = -1
            fac_temp[fac_temp != -1] = 0
            fac += fac_temp
            fac_counter -= 1
 
        # update face counter
        fac_counter += 1
 
        # print status
        counter = np.floor(i / len(vtx) * 100).astype(int)
        if counter == n_step[c_step]:
            print("clean vertices: "+str(counter)+" %")
            c_step += 1

    # vertices and indices without singularities
    vtx = vtx[n_singularity != 0]
    
    return vtx, fac
