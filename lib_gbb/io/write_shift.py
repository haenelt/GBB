def write_shift(vtx_old, vtx_new, line_dir, path_output, name_output):
    """
    This function writes an overlay with final vertex shifts. The difference between new and old
    coordinates is computed. If line_dir is 3, one overlay for each axis is written. 
    Inputs:
        *vtx_old: target array of vertices.
        *vtx_new: source array of vertices.
        *line_dir: shift direction (0,1,2,3).
        *path_output: path where output is written.
        *name_output: basename of written output file without file extension.
        
    created by Daniel Haenelt
    Date created: 22-12-2019         
    Last modified: 14-05-2020
    """
    import os
    import sys
    import numpy as np
    import nibabel as nb

    # suffix
    suffix = ["x","y","z"]
    
    # define header
    header = nb.freesurfer.mghformat.MGHHeader()
    
    # distance between old and new points
    vtx_shift = vtx_new - vtx_old
    
    # write output
    if line_dir > 3 or line_dir < 0:
        sys.exit("Choose a valid line direction!")
    elif line_dir == 3:
        for i in range(3):
            output = nb.freesurfer.mghformat.MGHImage(vtx_shift[:,i], np.eye(4), header)
            nb.save(output, os.path.join(path_output,name_output+"_"+suffix+".mgh"))
    else:
        output = nb.freesurfer.mghformat.MGHImage(vtx_shift[:,line_dir], np.eye(4), header)
        nb.save(output, os.path.join(path_output,name_output+".mgh"))