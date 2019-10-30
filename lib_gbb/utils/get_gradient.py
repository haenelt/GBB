def get_gradient(input_vol, path_output, input_vein="", sigma_gaussian=3, grad_dir=1, kernel_size=3):
    """
    This function computes the gradient in one direction (phase-encoding direction) of a mean bold
    image. Optionally, a vein mask can be applied.
    Inputs:
        *input_vol: filename of input nifti.
        *path_output: path where output is written.
        *input_vein: binary vein mask (1: veins, 0: background).
        *sigma_gaussian: gaussian blurring kernel for vein masking.
        *grad_dir: direction of gradient calculation.
        *kernel_size: kernel size for gradient calculation.
        
    created by Daniel Haenelt
    Date created: 30-10-2019     
    Last modified: 30-10-2019
    """
    import os
    import numpy as np
    import nibabel as nb
    import cv2
    from scipy.ndimage.filters import gaussian_filter

    # load input
    vol = nb.load(input_vol)
    vol_array = vol.get_fdata()

    # mask veins if a vein mask is given
    if input_vein:
        vein_array = np.round(nb.load(input_vein).get_fdata()).astype(int)
        vol_array_blurred = gaussian_filter(vol_array, 
                                            sigma_gaussian, 
                                            order = 0, 
                                            output = None, 
                                            mode = 'reflect', 
                                            cval = 0.0, 
                                            truncate = 4.0)
        
        vol_array[vein_array == 1] = vol_array_blurred[vein_array == 1]

    res = np.zeros_like(vol_array)
    for i in range(np.shape(vol_array)[2]):
        if grad_dir == 0:
            res[:,:,i] = cv2.Sobel(vol_array[:,:,i],cv2.CV_64F,0,1,kernel_size)
        elif grad_dir == 1:
            res[:,:,i] = cv2.Sobel(vol_array[:,:,i],cv2.CV_64F,1,0,kernel_size)
        else:
            print("invalid gradient direction!")
    
    # write gradient image
    output = nb.Nifti1Image(res,vol.affine,vol.header)
    nb.save(output,os.path.join(path_output,"gradient.nii"))
    
    output = nb.Nifti1Image(vol_array,vol.affine,vol.header)
    nb.save(output,os.path.join(path_output,"mean_data_blurred.nii"))