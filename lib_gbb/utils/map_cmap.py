def map_cmap(input_surf, input_ref, hemi, path_output, write_txt=True, cleanup=True):
    """
    This function computes for each vertex the corresponding nearest voxel location. Freesurfer has
    to be included to run this function.
    Inputs:
        *input_surf: filename of input surface.
        *input_ref: reference volume in the same space as the input surface (nifti format).
        *hemi: hemisphere.
        *path_output: path where output is written.
        *cleanup: remove intermediate files.
    Outputs:
        *vtx2vox: list containing voxel coordinates for each vertex.
    
    created by Daniel Haenelt
    Date created: 30-10-2019
    Last modified: 12-11-2019
    """
    import os
    import shutil as sh
    import numpy as np
    import nibabel as nb
    from nibabel.freesurfer.io import read_geometry
    from nipype.interfaces.freesurfer.preprocess import MRIConvert
    from nipype.interfaces.freesurfer import SampleToSurface
    from lib.cmap import generate_coordinate_mapping

    # set freesurfer path environment
    os.environ["SUBJECTS_DIR"] = path_output

    # freesurfer subject
    tmp = np.random.randint(0, 10, 5)
    tmp_string = ''.join(str(i) for i in tmp)
    sub = "tmp_"+tmp_string

    # make output folder
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # mimic freesurfer folder structure (with some additional folder for intermediate files)
    path_sub = os.path.join(path_output,sub)
    path_mri = os.path.join(path_sub,"mri")
    path_surf = os.path.join(path_sub,"surf")

    os.makedirs(path_sub)
    os.makedirs(path_mri)
    os.makedirs(path_surf)

    # copy input surface to mimicked freesurfer folde
    sh.copyfile(input_surf, os.path.join(path_surf,hemi+".source"))

    # create orig.mgz from reference volume
    mc = MRIConvert()
    mc.inputs.in_file = input_ref
    mc.inputs.out_file = os.path.join(path_mri,"orig.mgz")
    mc.inputs.in_type = "nii"
    mc.inputs.out_type = "mgz"
    mc.run()    

    # read surface geometry
    vtx, _ = read_geometry(input_surf)

    # generate cmap and save single dimensions
    cmap_img = generate_coordinate_mapping(input_ref, 
                                           0, 
                                           path_output = "", 
                                           suffix = "", 
                                           time = False, 
                                           write_output = False)
    cmap_img.header["dim"][0] = 3
    cmap_img.header["dim"][4] = 1
    cmap_array = cmap_img.get_fdata()

    components = ["x", "y", "z"]
    for i in range(len(components)):
        output = nb.Nifti1Image(cmap_array[:,:,:,i], cmap_img.affine, cmap_img.header)
        nb.save(output,os.path.join(path_mri,components[i]+".nii"))

    # map cmap to surf and get voxel locations for each vertex
    vtx2vox = np.zeros([len(vtx),3])
    for i in range(len(components)):

        # mri_vol2surf
        sampler = SampleToSurface()
        sampler.inputs.subject_id = sub
        sampler.inputs.reg_header = True
        sampler.inputs.hemi = hemi
        sampler.inputs.source_file = os.path.join(path_mri,components[i]+".nii")
        sampler.inputs.surface = "source"
        sampler.inputs.sampling_method = "point"
        sampler.inputs.sampling_range = 0
        sampler.inputs.sampling_units = "mm"
        sampler.inputs.interp_method = "trilinear" # or trilinear
        sampler.inputs.out_type = "mgh"
        sampler.inputs.out_file = os.path.join(path_surf,hemi+"."+components[i]+"_sampled.mgh")
        sampler.run()
    
        data_img = nb.load(os.path.join(path_surf,hemi+"."+components[i]+"_sampled.mgh"))
        vtx2vox[:,i] = np.squeeze(data_img.get_fdata())

    if write_txt:
        np.savetxt(os.path.join(path_output,"vtx2vox.txt"), vtx2vox, "%.5f")
    
    # delete intermediate files
    if cleanup:
        sh.rmtree(path_sub, ignore_errors=True)
    
    return vtx2vox