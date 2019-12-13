def get_direction_all(input_source, input_ref, diff_dir=2):
    """
    This function finds the nearest reference vertex to a each source vertex, i.e., the reference 
    point with minimum euclidean distance. The distance along one axis is then computed by taking 
    the difference reference point minus source point. The sign of the distance is used to infer the 
    direction between white and pial surface. Both surface meshs do not have to coincide in their 
    faces or number of vertices. The distance for each source vertex is written to a file.
    Inputs:
        *input_source: filename of source mesh (white surface).
        *input_ref: filename of reference mesh (pial surface).
        *diff_dir: axis for distance calculation (0,1,2).
        
    created by Daniel Haenelt
    Date created: 13-12-2019         
    Last modified: 13-12-2019
    """
    import numpy as np
    import nibabel as nb
    from nibabel.freesurfer.io import read_geometry

    # load geometry
    vtx_source, _ = read_geometry(input_source)
    vtx_ref, _ = read_geometry(input_ref)

    r_dist = np.zeros((len(vtx_source),1,1))
    for i in range(len(vtx_source)):
        
        # get euclidean distance
        rx_temp = ( vtx_ref[:,0] - vtx_source[i,0] ) ** 2
        ry_temp = ( vtx_ref[:,1] - vtx_source[i,1] ) ** 2
        rz_temp = ( vtx_ref[:,2] - vtx_source[i,2] ) ** 2
        r_temp = np.sqrt(rx_temp + ry_temp + rz_temp)
        
        # reference index number closest to current source vertex
        ind_min = np.where(r_temp == np.min(r_temp))[0][0]
        
        # distance between points along specified axis
        r_dist[i] = vtx_ref[ind_min,diff_dir] - vtx_source[i,diff_dir]
        
        # normalize distances to [-1, 1]
        r_dist = 2 * ( r_dist - np.min(r_dist) ) / ( np.max(r_dist) - np.min(r_dist) ) - 1
        
    # write output
    output = nb.Nifti1Image(r_dist, np.eye(4), nb.Nifti1Header())
    nb.save(output, input_source+"_ndist_dir"+str(diff_dir))
