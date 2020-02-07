def update_mesh(vtx, vtx_shift, source_ind, nn_ind, shift_dir=2, l_rate=0.9):
    """
    This function updates the vertices based on a calculated vertex shift within a defined
    neighboorhood and specified learning rate. The learning rate should be a value between 0 and 1.
    Vertex coordinates are in ras space.
    Inputs:
        *vtx: array of vertices.
        *vtx_shift
        *source_ind: source index.
        *nn_ind: neighborhood indices.
        *l_rate: learning rate.
    Outputs:
        *vtx_new: updated array of vertices.
        
    created by Daniel Haenelt
    Date created: 22-12-2019      
    Last modified: 22-12-2019
    """
    import numpy as np

    # get distance to source vertex
    sx = ( vtx[nn_ind,0] - vtx[source_ind,0] ) ** 2
    sy = ( vtx[nn_ind,1] - vtx[source_ind,1] ) ** 2
    sz = ( vtx[nn_ind,2] - vtx[source_ind,2] ) ** 2

    # scale shift depending on learning rate and neighborhood distance
    s = np.sqrt( sx + sy + sz )
    s = ( np.max(s) - s )/ np.max(s)

    vtx_new = vtx.copy()
    for i in range(len(nn_ind)):
        vtx_new[nn_ind[i],shift_dir] = vtx_new[nn_ind[i],shift_dir] - l_rate * s[i] * vtx_shift
    
    return vtx_new