def get_ignore(vtx, ignore_array, ras2vox_tkr, write_output=False, path_output=False):
    """
    This function gets vertex coordinates and vertex indices of all points which are within a binary 
    mask.
    Inputs:
        *vtx: array of vertex points.
        *ignore_array: 3D array with binary mask.
        *ras2vox_tkr: ras to voxel transformation matrix.
        *write_output: write output file (boolean).
        *path_output: path where output is written.
    Outputs:
        *vtx_ignore: array of vertex indices within mask.
        *ind_ignore: array of vertex points within mask.
        
    created by Daniel Haenelt
    Date created: 11-05-2020         
    Last modified: 11-05-2020
    """
    import os
    import numpy as np
    import nibabel as nb
    from nibabel.affines import apply_affine
    from lib_gbb.interpolation import nn_interpolation3d
    
    # transform vertices to voxel space
    vtx = apply_affine(ras2vox_tkr, vtx)
    
    # sample binary mask
    mask = nn_interpolation3d(vtx[:,0],vtx[:,1],vtx[:,2],ignore_array)
    
    # get all vertices within mask
    ind_ignore = np.where(mask == 1)[0].astype(int)
    vtx_ignore = vtx[ind_ignore,:]
    
    # write output
    if write_output:
        header = nb.freesurfer.mghformat.MGHHeader()
        output = nb.freesurfer.mghformat.MGHImage(ind_ignore, np.eye(4), header)
        nb.save(output, os.path.join(path_output, "ignore_mask.mgh"))
    
    return vtx_ignore, ind_ignore