def write_shift(vtx_new, vtx_old, path_output, name_output):
    """
    This function writes an overlay in freesurfer mgh format showing the euclidean distance between 
    source and target surface for each vertex. Vertex coordinates are in ras space.
    Inputs:
        *vtx_new: source array of vertices.
        *vtx_old: target array of vertices.
        *path_output: path where output is written.
        *name_output: basename of written output file.
        
    created by Daniel Haenelt
    Date created: 22-12-2019         
    Last modified: 22-12-2019
    """
    import os
    import numpy as np
    import nibabel as nb

    # euclidean distance between old and new points
    rx = ( vtx_old[:,0] - vtx_new[:,0] ) ** 2
    ry = ( vtx_old[:,1] - vtx_new[:,1] ) ** 2
    rz = ( vtx_old[:,2] - vtx_new[:,2] ) ** 2
    r = np.sqrt( rx + ry + rz )

    # write output
    header = nb.freesurfer.mghformat.MGHHeader()
    output = nb.freesurfer.mghformat.MGHImage(r, np.eye(4), header)
    nb.save(output, os.path.join(path_output,name_output+"_shift.mgh"))
