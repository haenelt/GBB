def plot_direction_euclidean(input_source, input_ref, diff_dir=2, diff_threshold=0.05):
    """
    This function plots the direction of from the white surface towards the pial surface based on
    the minimum euclidean distance between white and pial vertex points along one axis.
    Inputs:
        *input_source: filename of source mesh (white surface).
        *input_ref: filename of reference mesh (pial surface).
        *diff_dir: axis for distance calculation (0,1,2).
        *diff_threshold: threshold to ignore vertices in parallel to the considered axis.
        
    created by Daniel Haenelt
    Date created: 13-12-2019         
    Last modified: 13-12-2019
    """
    import numpy as np
    import nibabel as nb
    from nibabel.freesurfer.io import read_geometry
    from lib_gbb.normal import get_direction_euclidean

    # load geometry
    vtx_source, _ = read_geometry(input_source)
    vtx_ref, _ = read_geometry(input_ref)

    r_dist = np.zeros((len(vtx_source),1,1))
    for i in range(len(vtx_source)):
        r_dist[i] = get_direction_euclidean(vtx_source[i,:], vtx_ref, diff_dir, diff_threshold)

        
    # write output
    output = nb.Nifti1Image(r_dist, np.eye(4), nb.Nifti1Header())
    nb.save(output, input_source+"_plot_euclidean_dir"+str(diff_dir)+".mgh")