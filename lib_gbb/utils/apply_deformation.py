def apply_deformation(vtx, fac, arr_deform, vox2ras_tkr, ras2vox_tkr, path_output="", 
                      name_output="", write_output=True):
    """
    This function applies the coordinate shift to an array of vertices. The coordinate shift is 
    taken from a deformation field where each voxel corresponds to a shift along one direction in 
    voxel space.
    Inputs:
        *vtx: array of vertices.
        *fac: array of faces.
        *vox2ras_tkr: voxel to ras space transformation.
        *ras2vox_tkr: ras to voxel space transformation.
        *arr_deform: 4D volume array containing shifts in voxel space.
        *path_output: path where output is written.
        *name_output: name of output file.
        *write_output: write output file (boolean).
    Outputs:
        *vtx: deformed vertices.
        *fac: original faces.
        
    created by Daniel Haenelt
    Date created: 28-12-2019
    Last modified: 14-05-2020
    """
    import os
    import numpy as np
    from nibabel.freesurfer.io import write_geometry
    from nibabel.affines import apply_affine
    from lib_gbb.interpolation.linear_interpolation3d import linear_interpolation3d
    
    # get vertices to voxel space
    vtx = apply_affine(ras2vox_tkr, vtx)
    
    # sample shifts from deformation map
    vtx_shift = np.zeros((len(vtx),3))
    for i in range(3):
        vtx_shift[:,i] = linear_interpolation3d(vtx[:,0], vtx[:,1], vtx[:,2], arr_deform[:,:,:,i])
    
    # shift coordinates to new location
    vtx += vtx_shift
    
    # get new coordinates in ras space
    vtx = apply_affine(vox2ras_tkr, vtx)
    
    if write_output:
        write_geometry(os.path.join(path_output,name_output), vtx, fac)
    
    return vtx, fac