def write_shift(vtx_new, vtx_old, line_dir, path_output, name_output):
    """
    This function writes an overlay with final vertex shifts. The difference between new and old
    coordinates is computed. Since they differ only along one axis, we can just sum across columns.
    This gets rid of the zeros in all other axes.
    Inputs:
        *vtx_new: source array of vertices.
        *vtx_old: target array of vertices.
        *line_dir: shift direction.
        *path_output: path where output is written.
        *name_output: basename of written output file without file extension.
        
    created by Daniel Haenelt
    Date created: 22-12-2019         
    Last modified: 28-12-2019
    """
    import os
    import numpy as np
    import nibabel as nb

    # distance between old and new points
    vtx_shift = vtx_new - vtx_old
    
    # get shift along one direction
    vtx_shift = vtx_shift[:,line_dir]

    # write output
    header = nb.freesurfer.mghformat.MGHHeader()
    output = nb.freesurfer.mghformat.MGHImage(vtx_shift, np.eye(4), header)
    nb.save(output, os.path.join(path_output,name_output+".mgh"))