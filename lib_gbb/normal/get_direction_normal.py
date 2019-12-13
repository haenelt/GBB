def get_direction_normal(ind_source, vtx_source, fac_source, diff_dir=2, diff_threshold=0.05):
    """
    This function gets the length of the vertex normal along one axis. The sign is used tp infer the
    direction between white and pial surface.
    Inputs:
        *ind_source: index of source vertex point (white surface).
        *vtx_source: array of source vertices (white surface).
        *fac_source: array of source faces (white surface).
        *diff_dir: axis for distance calculation (0,1,2).
        *diff_threshold: threshold to ignore vertices in parallel to the considered axis.
    Outputs:
        *r_dist: normalized direction to reference mesh (pial surface).
        
    created by Daniel Haenelt
    Date created: 13-12-2019      
    Last modified: 13-12-2019
    """
    import numpy as np
    from lib_gbb.normal import get_normal_surf
    
    # get surface normals per vertex
    norm = get_normal_surf(vtx_source, fac_source)
    
    # norm length along one axis
    r_dist = norm[ind_source, diff_dir]
    
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