def deformation_field(vtx_old, vtx_new, input_vol, vox2ras_tkr, ras2vox_tkr, sigma=1, 
                      path_output="", name_output="", write_output=True):
    """
    This function computes a deformation field from a array of shifted vertices. Each voxel in the 
    deformation field corresponds to a shift in voxel space along one direction. Optionally, a 
    gaussian filter can be applied to the final deformation field.
    Inputs:
        *vtx_old: original array of vertices.
        *vtx_new: new array of vertices.
        *input_vol: filename of reference volume.
        *vox2ras_tkr: voxel to ras space transformation.
        *ras2vox_tkr: ras to voxel space transformation.
        *sigma: sigma for gaussian filtering of deformation field.
        *path_output: path where output is written.
        *name_output: basename of output file without file extension.
        *write_output: write output file (boolean).
    Outputs:
        *arr_deform: deformation field.
        
    created by Daniel Haenelt
    Date created: 28-12-2019       
    Last modified: 14-05-2020
    """
    import os
    import numpy as np
    import nibabel as nb
    from nibabel.affines import apply_affine
    from scipy.interpolate import griddata
    from scipy.ndimage import gaussian_filter
    from lib_gbb.utils.line_ras2vox import line_ras2vox
    
    # line directions
    line_dir = [0,1,2]
    line_vox = []
    
    # load reference volume
    vol = nb.load(input_vol)
    vol.header["dim"][4] = 3
    
    # get volume dimension
    x_dim = vol.header["dim"][1]
    y_dim = vol.header["dim"][2]
    z_dim = vol.header["dim"][3]
    
    # initialize deformation array
    arr_deform = np.zeros(vol.header["dim"][1:5])
    
    # get vertices to voxel space
    vtx_old = apply_affine(ras2vox_tkr, vtx_old)
    vtx_new = apply_affine(ras2vox_tkr, vtx_new)
    
    # get line directions to voxel space
    line_vox = [line_ras2vox(line_dir[i], ras2vox_tkr) for i in line_dir]
    
    # get vertex shifts in voxel space
    vtx_shift = np.zeros((len(vtx_old),3))
    for i in range(len(line_vox)):
        vtx_shift[:,line_vox[i]] = vtx_new[:,line_vox[i]] - vtx_old[:,line_vox[i]]
    
    # get volume coordinates in voxel space
    xf = np.arange(0,vol.header["dim"][1])
    yf = np.arange(0,vol.header["dim"][2])
    zf = np.arange(0,vol.header["dim"][3])
    
    x_plane, y_plane, z_plane = np.meshgrid(xf,yf,zf)
    x_plane = x_plane.flatten()
    y_plane = y_plane.flatten()
    z_plane = z_plane.flatten()
    
    vox_coords = np.column_stack([y_plane,x_plane,z_plane])
    
    # grid interpolation
    nn = griddata(vtx_old, vtx_shift, vox_coords, method='nearest')
    li = griddata(vtx_old, vtx_shift, vox_coords, method='linear')
    
    # fill all nans from the linear interpolation with nearest neighbor values
    li[np.isnan(li)] = nn[np.isnan(li)]
    
    # reshape to deformation volume
    for i in range(len(line_vox)):
        
        # reshape
        arr_deform[:,:,:,i] = np.reshape(nn[:,i],(y_dim,x_dim,z_dim))
    
        # apply gaussian filter
        if sigma:
            arr_deform[:,:,:,i] = gaussian_filter(arr_deform[:,:,:,i], sigma)
        
    # write output    
    if write_output:
        output = nb.Nifti1Image(arr_deform, vol.affine, vol.header)
        nb.save(output, os.path.join(path_output,name_output+".nii"))
    
    return arr_deform