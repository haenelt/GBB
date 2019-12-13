def get_direction(vtx_source, vtx_ref, diff_dir=2, diff_threshold=0.05):
    """
    This function finds the nearest reference vertex to a source point, i.e., the reference point
    with minimum euclidean distance. The distance along one axis is then computed between both 
    points by taking the difference reference point minus source point. The sign of the distance is 
    used to infer the direction between white and pial surface.
    Inputs:
        *vtx_source: source vertex point (white surface).
        *vtx_ref: array of reference vertices (pial surface).
        *diff_dir: axis for distance calculation (0,1,2).
        *diff_threshold: threshold to ignore vertices in parallel to the considered axis.
    Outputs:
        *r_dist: normalized direction to reference mesh.
        
    created by Daniel Haenelt
    Date created: 13-12-2019        
    Last modified: 13-12-2019
    """
    import numpy as np
    
    # get euclidean distance
    rx_temp = ( vtx_ref[:,0] - vtx_source[0] ) ** 2
    ry_temp = ( vtx_ref[:,1] - vtx_source[1] ) ** 2
    rz_temp = ( vtx_ref[:,2] - vtx_source[2] ) ** 2
    r_temp = np.sqrt(rx_temp + ry_temp + rz_temp)
       
    # reference index number closest to current source vertex
    ind_min = np.where(r_temp == np.min(r_temp))[0][0]
       
    # distance between points along specified axis
    r_dist = vtx_ref[ind_min,diff_dir] - vtx_source[diff_dir]
    
    # get direction
    if diff_threshold > 0:
        
        if r_dist > diff_threshold:
            r_dist = 1
        elif np.abs(r_dist) <= diff_threshold:
            r_dist = 0
        else:
            r_dist = -1
    
    else:
        
        if r_dist > 0:
            r_dist = 1
        else:
            r_dist = -1

    return r_dist