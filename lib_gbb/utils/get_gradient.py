def get_gradient(input_vol, ras2vox, line_dir=2, sigma=0, kernel_size=3, write_output=None):
    """
    This function computes the second order gradient in one direction (phase-encoding direction) of 
    a mean bold image.
    Inputs:
        *input_vol: filename of input nifti.
        *ras2vox: transformation matrix from ras to voxel space.
        *line_dir: direction of gradient calculation in ras space.
        *sigma: gaussian blurring kernel of gradient map (if set > 0).
        *kernel_size: kernel size for gradient calculation.
        *write_output: write output image (boolean).
    Outputs:
        *res: gradient array.
        
    created by Daniel Haenelt
    Date created: 30-10-2019     
    Last modified: 23-12-2019
    """
    import sys
    import os
    import cv2
    import numpy as np
    import nibabel as nb
    from nibabel.affines import apply_affine
    from scipy.ndimage.filters import gaussian_filter

    # get line direction in ras space
    if line_dir == 0:
        pt_ras = [1,0,0]
    elif line_dir == 1:
        pt_ras = [0,1,0]
    elif line_dir == 2:
        pt_ras = [0,0,1]
    else:
        sys.exit("Invalid axis direction in gradient calculation!")

    # get unit vector in voxel space
    pt0_vox = apply_affine(ras2vox, [0,0,0])
    pt_vox = apply_affine(ras2vox, pt_ras)

    pt_vox = np.abs( pt_vox - pt0_vox )
    pt_vox = pt_vox / np.max(pt_vox)
    pt_vox = pt_vox.astype(int)

    # updated line direction in voxel space
    line_dir = np.where(pt_vox == 1)[0][0]

    # load input
    vol = nb.load(input_vol)
    vol_array = vol.get_fdata()

    # get gradient
    res = np.zeros_like(vol_array)
    if line_dir == 0:
        for i in range(np.shape(vol_array)[2]):
            res[:,:,i] = cv2.Sobel(vol_array[:,:,i],cv2.CV_64F,0,1,kernel_size)
            res[:,:,i] = cv2.Sobel(res[:,:,i],cv2.CV_64F,0,1,kernel_size)
    
    elif line_dir == 1:
        for i in range(np.shape(vol_array)[2]):
            res[:,:,i] = cv2.Sobel(vol_array[:,:,i],cv2.CV_64F,1,0,kernel_size)
            res[:,:,i] = cv2.Sobel(res[:,:,i],cv2.CV_64F,1,0,kernel_size)        
    
    else:
        for i in range(np.shape(vol_array)[0]):
            res[i,:,:] = cv2.Sobel(vol_array[i,:,:],cv2.CV_64F,1,0,kernel_size)
            res[i,:,:] = cv2.Sobel(res[i,:,:],cv2.CV_64F,1,0,kernel_size) 
    
    # gaussian blurring (optional)
    if sigma:
        res = gaussian_filter(res, 
                              sigma, 
                              order = 0, 
                              output = None, 
                              mode = 'reflect', 
                              cval = 0.0, 
                              truncate = 4.0)
    
    # write gradient image
    if write_output:
        output = nb.Nifti1Image(res, vol.affine, vol.header)
        nb.save(output,os.path.join(os.path.dirname(input_vol),"gradient.nii"))
    
    return res