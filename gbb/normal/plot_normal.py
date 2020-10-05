# -*- coding: utf-8 -*-

# external inputs
import numpy as np
from nibabel.freesurfer.io import write_geometry

# local inputs
from gbb.neighbor.nn_2d import nn_2d
from gbb.normal import get_normal
    

def plot_normal(vtx, fac, adjm, file_out, step_size=100, shape="line"):
    """
    This function generates lines to visualize outward directed surface normals 
    of an input surface mesh.
    Inputs:
        *vtx (arr): array of vertex points.
        *fac (arr): corresponding face array.
        *adjm (obj): adjacency matrix.
        *file_out (str): filename of output surface mesh.
        *step_size (int): subset of vertices.
        *shape (str): line, triangle, prism
        
    created by Daniel Haenelt
    Date created: 13-12-2019
    Last modified: 05-10-2020
    """
    
    # get surface normals
    normal = get_normal(vtx, fac)
    
    # array containing a list of considered vertices
    t = np.arange(0,len(vtx),step_size)
    
    # initialise faces for specific shape
    if shape == "prism":
        fac_new = [[0, 1, 2],
                   [3, 4, 5],
                   [0, 1, 4],
                   [0, 3, 4],
                   [1, 2, 5],
                   [1, 4, 5],
                   [0, 2, 5],
                   [0, 3, 5]]
        fac_iter = 6
    elif shape == "triangle":
        fac_new = [[0,1,2]]
        fac_iter = 3
    elif shape == "line":
        fac_new = [[0,1,0]]
        fac_iter = 2
    
    vtx_res = []
    fac_res = []
    for i in range(len(t)):
        
        # get index from nearest neighbour of a given vertex
        nn = nn_2d(t[i], adjm, 0)
        nn = nn[:2]
        
        # get all vertex points for specific shape
        if shape == "prism":
            A = list(vtx[t[i]])
            B = list(vtx[nn[0]])
            C = list(vtx[nn[1]])
            D = list(vtx[t[i]] - normal[t[i]])
            E = list(vtx[nn[0]] - normal[nn[0]])
            F = list(vtx[nn[1]] - normal[nn[1]])
            vtx_new = [A,B,C,D,E,F]
        elif shape == "triangle":
            A = list(vtx[t[i]])
            B = list(vtx[nn[0]])
            C = list(vtx[t[i]] - normal[t[i]])
            vtx_new = [A,B,C]
        elif shape == "line":
            A = list(vtx[t[i]])
            B = list(vtx[t[i]] - normal[t[i]])
            vtx_new = [A,B]
    
        # update faces
        if i > 0:
            for j in range(len(fac_new)):
                fac_new[j] = [x+fac_iter for x in fac_new[j]]
        
        # update resulting vertex and face list
        vtx_res.extend(vtx_new)
        fac_res.extend(fac_new)
    
    # vertices and faces as array
    vtx_res = np.array(vtx_res)
    fac_res = np.array(fac_res)
    
    # write output geometry
    write_geometry(file_out, vtx_res, fac_res)
    