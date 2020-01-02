def plot_normal_direction(input_surf, diff_dir=2, diff_threshold=0.05):
    """
    This function plots the direction from the white surface towards the pial surface based on the 
    white surface vertex normals along one axis.
    Inputs:
        *input_surf: filename of source mesh (white surface).
        *diff_dir: axis for distance calculation in ras space.
        *diff_threshold: threshold to ignore normals perpendicular to the considered axis.
        
    created by Daniel Haenelt
    Date created: 13-12-2019         
    Last modified: 20-12-2019
    """
    import numpy as np
    import nibabel as nb
    from nibabel.freesurfer.io import read_geometry
    from lib_gbb.normal import get_normal_direction

    # load geometry
    vtx_source, fac_source = read_geometry(input_surf)

    # get directions 
    r_dist, _ = get_normal_direction(vtx_source, fac_source, diff_dir, diff_threshold)

    # write output
    header = nb.freesurfer.mghformat.MGHHeader()
    output = nb.freesurfer.mghformat.MGHImage(r_dist, np.eye(4), header)
    nb.save(output, input_surf+"_plot_normal_dir"+str(diff_dir)+".mgh")