def plot_direction_normal(input_source, diff_dir=2, diff_threshold=0.05):
    """
    This function plots the direction from the white surface towards the pial surface based on the 
    white surface vertex normals along one axis.
    Inputs:
        *input_source: filename of source mesh (white surface).
        *diff_dir: axis for distance calculation (0,1,2).
        *diff_threshold: threshold to ignore vertices in parallel to the considered axis.
        
    created by Daniel Haenelt
    Date created: 13-12-2019         
    Last modified: 13-12-2019
    """
    import numpy as np
    import nibabel as nb
    from nibabel.freesurfer.io import read_geometry
    from lib_gbb.normal import get_direction_normal

    # load geometry
    vtx_source, fac_source = read_geometry(input_source)

    r_dist = np.zeros((len(vtx_source),1,1))
    for i in range(len(vtx_source)):
        r_dist[i] = get_direction_normal(i, vtx_source, fac_source, diff_dir, diff_threshold)

        
    # write output
    output = nb.Nifti1Image(r_dist, np.eye(4), nb.Nifti1Header())
    nb.save(output, input_source+"_plot_normal_dir"+str(diff_dir)+".mgh")