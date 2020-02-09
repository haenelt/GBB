def devein_mesh(surf_in, ref_in, vein_in, surf_out=None, n_neighbor=20, shift_dir=2, smooth_iter=30, 
                max_iterations=1000):
    """
    This function finds vertex points which are located within marked veins and shift these and 
    their neighborhood until the mesh is free of trapped vertices (or the maximum number of
    iterations is reached). Shifts will be applied only along one axis and in inward direction. 
    Optionally, the output mesh can be smoothed.    
    Inputs.
        *surf_in: filename of input surface.
        *ref_in: filename of reference volume.
        *vein_in: filename of vein mask.
        *surf_out: filename of output surface.
        *n_neighbor: neighborhood size.
        *shift_dir: shift direction in ras conventions.
        *smooth_iter: number of smoothing iterations of final surface mesh.
        *max_iterations: maximum number of deveining iterations.
    Outputs:
        *vtx: shifted array of vertex points.
    
    created by Daniel Haenelt
    Date created: 06-02-2020             
    Last modified: 09-02-2020  
    """
    import os
    import numpy as np
    import nibabel as nb
    from nibabel.freesurfer.io import read_geometry, write_geometry
    from nibabel.affines import apply_affine
    from lib.surface.vox2ras import vox2ras
    from lib.surface.smooth_surface import smooth_surface
    from lib_gbb.utils.get_adjm import get_adjm
    from lib_gbb.utils.update_mesh import update_mesh
    from lib_gbb.neighbor.nn_2d import nn_2d
    from lib_gbb.normal.get_normal_direction import get_normal_direction

    # load data
    vein = np.round(nb.load(vein_in).get_fdata())
    vtx, fac = read_geometry(surf_in)
    adjm = get_adjm(surf_in)
    _, ras2vox_tkr = vox2ras(ref_in)

    # get nearest voxel coordinates
    vtx_vox = apply_affine(ras2vox_tkr, vtx)
    vtx_vox = np.round(vtx_vox).astype(int)   

    # get vertices trapped in veins
    vein_mask = np.zeros(len(vtx))
    vein_mask = vein[vtx_vox[:,0], vtx_vox[:,1], vtx_vox[:,2]]
    n_veins = len(vein_mask[vein_mask == 1])

    # get surface normals
    norm, _ = get_normal_direction(vtx, fac, shift_dir, 0.05)

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
        if norm[curr_ind] < 0:
            vtx_shift = -1 * vtx_shift
        elif norm[curr_ind] == 0:
            vein[vtx_vox[curr_ind,0],vtx_vox[curr_ind,1],vtx_vox[curr_ind,2]] = 0
            continue
    
        # update mesh
        vtx = update_mesh(vtx, vtx_shift, curr_ind, nn_ind, 1)
        
        vtx_vox = apply_affine(ras2vox_tkr, vtx)
        vtx_vox = np.round(vtx_vox).astype(int)    
    
        # get all vertices within vein
        vein_mask = np.zeros(len(vtx))
        vein_mask = vein[vtx_vox[:,0], vtx_vox[:,1], vtx_vox[:,2]]
        n_veins = len(vein_mask[vein_mask == 1])
        
        counter += 1
    
    if counter < max_iterations:
        print("Deveining converged!")
    
    # smooth output
    if smooth_iter:
        surf_temp = os.path.join(os.path.dirname(surf_out),"surf_temp")
        write_geometry(surf_temp, vtx, fac)
        smooth_surface(surf_temp, surf_temp, smooth_iter)
        vtx, _ = read_geometry(surf_temp)
        os.remove(surf_temp)
     
    # write output
    if surf_out:
        write_geometry(surf_out, vtx, fac) 
        
    return vtx