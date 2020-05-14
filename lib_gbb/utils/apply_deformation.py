def apply_deformation(vtx, fac, vox2ras_tkr, ras2vox_tkr, line_dir, deformation, path_output="", 
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
        *line_dir: shift direction in ras space.
        *deformation: 3D volume array containing shifts in voxel space.
        *path_output: path where output is written.
        *name_output: name of output file.
        *write_output: write output file (boolean).
    Outputs:
        *vtx: deformed vertices.
        *fac: original faces.
        
    created by Daniel Haenelt
    Date created: 28-12-2019
    Last modified: 28-12-2019
    """
    import os
    from nibabel.freesurfer.io import write_geometry
    from nibabel.affines import apply_affine
    from lib_gbb.interpolation.linear_interpolation3d import linear_interpolation3d
    from lib_gbb.utils.line_ras2vox import line_ras2vox

    # get line direction in voxel space
    line_dir = line_ras2vox(line_dir, ras2vox_tkr)

    # transform vertices
    vtx = apply_affine(ras2vox_tkr, vtx)
    
    # sample shifts from deformation map
    shift = linear_interpolation3d(vtx[:,0], vtx[:,1], vtx[:,2], deformation)

    # shift coordinates to new location along line direction
    vtx[:,line_dir] = vtx[:,line_dir] + shift
    
    # get new coordinates in ras space
    vtx = apply_affine(vox2ras_tkr, vtx)
    
    if write_output:
        write_geometry(os.path.join(path_output,name_output), vtx, fac)
    
    return vtx, fac