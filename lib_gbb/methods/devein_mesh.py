def devein_mesh(vtx, fac, arr_vein, arr_ignore, n, adjm, ras2vox, n_neighbor=20, shift_dir=2, 
                smooth_iter=30, max_iterations=1000):
    """
    This function finds vertex points which are located within marked veins and shift these and 
    their neighborhood until the mesh is free of trapped vertices (or the maximum number of
    iterations is reached). Shifts will be applied only along one axis and in inward direction. 
    Optionally, the output mesh can be smoothed.    
    Inputs.
        *vtx: array of vertex points.
        *fac: array of corresponding faces.
        *arr_vein: vein mask.
        *arr_ignore: binary mask where deveining is omitted.
        *n: array of vertex normal directions.
        *adjm: adjacecy matrix.
        *ras2vox: transformation matrix from ras to voxel space.
        *n_neighbor: neighborhood size.
        *shift_dir: shift direction in ras conventions.
        *smooth_iter: number of smoothing iterations of final surface mesh.
        *max_iterations: maximum number of deveining iterations.
    Outputs:
        *vtx: shifted array of vertex points.
        *counter: number of deveining iterations.
    
    created by Daniel Haenelt
    Date created: 06-02-2020             
    Last modified: 13-05-2020  
    """
    import os
    import numpy as np
    from numpy.linalg import norm
    from nibabel.freesurfer.io import read_geometry, write_geometry
    from nibabel.affines import apply_affine
    from lib.surface.smooth_surface import smooth_surface
    from lib_gbb.utils.update_mesh import update_mesh
    from lib_gbb.neighbor.nn_2d import nn_2d

    # load arrays
    arr_vein = np.round(arr_vein).astype(int)
    
    if arr_ignore is not None:
        arr_ignore = np.round(arr_ignore).astype(int)
        arr_vein[arr_ignore == 1] = 0
   
    # centroid
    vtx_c = np.mean(vtx, axis=0)

    # get nearest voxel coordinates
    vtx_vox = apply_affine(ras2vox, vtx)
    vtx_vox = np.round(vtx_vox).astype(int)   

    # get vertices trapped in veins
    vein_mask = np.zeros(len(vtx))
    vein_mask = arr_vein[vtx_vox[:,0], vtx_vox[:,1], vtx_vox[:,2]]
    n_veins = len(vein_mask[vein_mask == 1])

    print("start mesh initialization (deveining)")

    # loop through vertices
    counter = 0
    while n_veins > 0 and counter < max_iterations:
    
        # print current status
        print("i: "+str(counter)+", # of trapped vertices: "+str(len(vein_mask[vein_mask == 1])))
    
        # select random vein vertex
        vein_ind = np.where(vein_mask == 1)[0]
        curr_ind = vein_ind[np.random.randint(len(vein_ind))]
        nn_ind = nn_2d(curr_ind, adjm, n_neighbor)
    
        # get shift as weighted average
        vtx_shift = np.zeros((len(nn_ind),3))
        vtx_shift[:,shift_dir] = vtx[nn_ind,shift_dir]
        vtx_shift[:,shift_dir] = vtx[curr_ind,shift_dir] - vtx_shift[:,shift_dir]
        vtx_shift = np.mean(vtx_shift, axis=0)
        vtx_shift = np.abs(vtx_shift)
        
        # do only inward shifts
        if n[curr_ind] < 0:
            vtx_shift = -1 * vtx_shift
        
        # check inward shift by comparing distance to centroid
        vtx_dist_noshift = norm(vtx[curr_ind,:] - vtx_c)
        vtx_dist_shift = norm(vtx[curr_ind,:] - vtx_shift - vtx_c)
        
        # update mesh if valid inward shift        
        if vtx_dist_shift - vtx_dist_noshift > 0 or n[curr_ind] == 0:
            vtx_temp_vox = apply_affine(ras2vox, vtx[curr_ind])
            vtx_temp_vox = np.round(vtx_temp_vox).astype(int)
            arr_vein[vtx_temp_vox[0],vtx_temp_vox[1],vtx_temp_vox[2]] = 0
        else:
            vtx = update_mesh(vtx, vtx_shift, curr_ind, nn_ind, 1)
            vtx_vox = apply_affine(ras2vox, vtx)
            vtx_vox = np.round(vtx_vox).astype(int)    
    
        # get all vertices within vein
        vein_mask = np.zeros(len(vtx))
        vein_mask = arr_vein[vtx_vox[:,0], vtx_vox[:,1], vtx_vox[:,2]]
        n_veins = len(vein_mask[vein_mask == 1])
        
        counter += 1
    
    if counter < max_iterations:
        print("Deveining converged!")
    
    # smooth output
    if smooth_iter:
        tmp = np.random.randint(0, 10, 5)
        tmp_string = ''.join(str(i) for i in tmp)
        surf_temp = os.path.join(os.getcwd(),"surf_"+tmp_string)
        write_geometry(surf_temp, vtx, fac)
        smooth_surface(surf_temp, surf_temp, smooth_iter)
        vtx, _ = read_geometry(surf_temp)
        os.remove(surf_temp)
        
    return vtx, counter