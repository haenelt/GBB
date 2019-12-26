def get_normal_direction(vtx_source, fac_source, diff_dir=2, diff_threshold=0.05):
    """
    This function gets the length of the vertex normal along one axis. The sign is used to infer the
    direction between white and pial surface. Vertices are in ras space.
    Inputs:
        *vtx_source: array of source vertices (white surface).
        *fac_source: array of source faces (white surface).
        *diff_dir: axis for distance calculation in ras.
        *diff_threshold: threshold to ignore vertices in parallel to the considered axis.
    Outputs:
        *r_dist: normalized directions to reference mesh (pial surface).
        *norm: normal per vertex.
        
    created by Daniel Haenelt
    Date created: 13-12-2019      
    Last modified: 18-12-2019
    """
    import numpy as np
    from lib_gbb.normal import get_normal
    
    # get surface normals per vertex
    norm = get_normal(vtx_source, fac_source)
    
    # get distance along one axis
    r_dist = norm[:, diff_dir]
    
    # get directions
    if diff_threshold > 0:
        r_dist[r_dist > diff_threshold] = 1
        r_dist[r_dist < -diff_threshold] = -1
        r_dist[np.abs(r_dist) != 1] = 0
    else:
        r_dist[r_dist > 0] = 1
        r_dist[r_dist != 1] = -1

    return r_dist, norm