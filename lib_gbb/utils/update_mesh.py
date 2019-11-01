def update_mesh(input_surf, vtx_n, vtx_shift, n_iter=2, dir=2, sigmoid_range=5, write_output=False):
    """
    This function applies a shift to a vertex in a specified direction and updates vertex positions
    in the local neighbourhood. The neighbourhood is shifted in dependence of its euclidean distance
    to the shifted center vertex. A sigmoid filter is applied to dampen the shift towards the
    periphery.
    Inputs:
        *input_surf: filename of surface mesh.
        *vtx_n: vertex number of shifted vertex.
        *vtx_shift: vertex shift in mm.
        *n_iter: number of iterations to collect neighbourhood vertices.
        *dir: direction of vertex shift.
        *sigmoid_range: strength of sigmoid filter.
        *write_output: write updated surface mesh (optional).
    Outputs:
        *vtx_res: updated vertex coordinates.
        
    created by Daniel Haenelt
    Date created: 01-11-2019      
    Last modified: 01-11-2019
    """
    import os
    import numpy as np
    import math
    import datetime
    from nibabel.freesurfer.io import read_geometry, write_geometry
    from numpy.linalg import norm

    def sigmoid(x):
        return 1 / (1 + math.exp(-x))

    # load surface mesh
    vtx, fac = read_geometry(input_surf)

    # get n_iter neighbourhood
    neighbour_n = np.unique(np.concatenate(fac[np.where(fac==vtx_n)[0],:]))
    if n_iter > 1:
        for i in range(n_iter):
            temp_n = []
            for j in range(len(neighbour_n)):
                temp_n = np.concatenate((temp_n,np.unique(np.concatenate(fac[np.where(fac==neighbour_n[j])[0],:]))))

            neighbour_n = np.concatenate((neighbour_n,temp_n))
            neighbour_n = np.unique(neighbour_n).astype(int)

    # remove vtx_n from neighbourhood
    neighbour_n = neighbour_n[neighbour_n != vtx_n]

    # do one shift
    vtx_old = vtx[vtx_n].copy()
    vtx_new = vtx[vtx_n].copy()
    vtx_new[dir] += vtx_shift

    # initialise new array
    vtx_res = vtx.copy()
    vtx_res[vtx_n] = vtx_new

    # get change in euclidean distance for each neighbour
    d = np.zeros(len(neighbour_n))
    for i in range(len(neighbour_n)):
        d[i] = norm(vtx[neighbour_n[i]] - vtx_old)

    # normalise distances in neighbourhood and apply sigmoid filter
    d = 2*( d - np.min(d) ) / (np.max(d)-np.min(d)) - 1 
    d = -1 * sigmoid_range * d

    # update shifts in neighbourhood
    for i in range(len(neighbour_n)):
        vtx_res[neighbour_n[i],dir] = vtx_res[neighbour_n[i],dir] + sigmoid(d[i]) * vtx_shift

    # write updated mesh
    if write_output:
        write_geometry(os.path.join(os.getcwd(),"upate_mesh_"+datetime.datetime.now().strftime("%Y%m%d%H%M%S")),vtx_res,fac)
    
    return vtx_res