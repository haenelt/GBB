def deformation_field(vtx_old, vtx_new, input_vol, line_ras, vox2ras_tkr, ras2vox_tkr, sigma=1, 
                      path_output="", name_output="", write_output=True):
    """
    This function computes a deformation field from a array of shifted vertices. Each voxel in the 
    deformation field corresponds to a shift in voxel space along one direction. Optionally, a 
    gaussian filter can be applied to the final deformation field.
    Inputs:
        *vtx_old: original array of vertices.
        *vtx_new: new array of vertices.
        *input_vol: filename of reference volume.
        *line_ras: shift direction in ras space.
        *vox2ras_tkr: voxel to ras space transformation.
        *ras2vox_tkr: ras to voxel space transformation.
        *sigma: sigma for gaussian filtering of deformation field.
        *path_output: path where output is written.
        *name_output: basename of output file without file extension.
        *write_output: write output file (boolean).
    Outputs:
        *deform_li: deformation field.
        
    created by Daniel Haenelt
    Date created: 28-12-2019       
    Last modified: 13-05-2020
    """
    import os
    import numpy as np
    import nibabel as nb
    from nibabel.affines import apply_affine
    from scipy.interpolate import griddata
    from scipy.ndimage import gaussian_filter
    from lib_gbb.utils.line_ras2vox import line_ras2vox

    # transform vertices to voxel space
    vtx_old = apply_affine(ras2vox_tkr, vtx_old)
    vtx_new = apply_affine(ras2vox_tkr, vtx_new)

    # line direction in voxel space
    line_vox = line_ras2vox(line_ras, ras2vox_tkr)

    # get vertex shift along one axis
    vtx_shift = vtx_new[:,line_vox] - vtx_old[:,line_vox]

    # load volume
    vol = nb.load(input_vol)
    arr_vol = vol.get_fdata()

    # initialize arrays
    deform_li = np.zeros_like(arr_vol)
    deform_nn = np.zeros_like(arr_vol)

    # get coordinate grid of one slice
    if line_vox != 2:
        xf = np.arange(0,np.shape(arr_vol)[0])
        yf = np.arange(0,np.shape(arr_vol)[1])
    else:
        xf = np.arange(0,np.shape(arr_vol)[1])
        yf = np.arange(0,np.shape(arr_vol)[2])
    
    y_plane, x_plane = np.meshgrid(yf,xf)

    # interpolate index values to grid
    x_plane_reshape = x_plane.reshape(len(xf)*len(yf),)
    y_plane_reshape = y_plane.reshape(len(xf)*len(yf),)

    # get slices in voxel space
    vtx_slice = np.round(vtx_old[:,line_ras]).astype(int)

    # get deformation field
    n_iter = np.unique(vtx_slice)
    for i in n_iter:
    
        # get one slice
        vtx_old_temp = vtx_old[vtx_slice == i]
        vtx_shift_temp = vtx_shift[vtx_slice == i]
        
        # delete slice column
        vtx_old_temp = np.delete(vtx_old_temp, line_ras, axis=1)

        # grid interpolation of index data
        coord_plane = np.transpose(np.array((x_plane_reshape, y_plane_reshape)))

        # grid interpolation
        nn = griddata(vtx_old_temp, vtx_shift_temp, coord_plane, method='nearest')
        li = griddata(vtx_old_temp, vtx_shift_temp, coord_plane, method='linear')

        if line_vox != 2:
            deform_nn[:,:,i] = nn.reshape(len(yf),len(xf)) 
            deform_li[:,:,i] = li.reshape(len(yf),len(xf)) 
        else:
            deform_nn[i,:,:] = nn.reshape(len(yf),len(xf)) 
            deform_li[i,:,:] = li.reshape(len(yf),len(xf))             
    
    # fill all nans from the linear interpolation with nearest neighbor values
    deform_li[np.isnan(deform_li)] = deform_nn[np.isnan(deform_li)]
    
    # apply gaussian filtering
    if sigma:
        deform_li = gaussian_filter(deform_li, sigma)

    # write output    
    if write_output:
        output = nb.Nifti1Image(deform_li, vol.affine, vol.header)
        nb.save(output,os.path.join(path_output,name_output+".nii"))
    
    return deform_li