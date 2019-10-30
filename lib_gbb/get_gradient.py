# get line
import numpy as np
import nibabel as nb
import cv2


# get gradient from input volume
# mask out veins

input_vol = "/home/raid2/haenelt/projects/GBB/test_data/mean_data.nii"
vol = nb.load(input_vol)
vol_array = vol.get_fdata()

#vol_grad = np.gradient(vol_array[:,:,30],axis = 1)

res = np.zeros_like(vol_array)
for i in range(np.shape(vol_array)[2]):
    res[:,:,i] = sobx = cv2.Sobel(vol_array[:,:,i],cv2.CV_64F,1,0,ksize=3)

output = nb.Nifti1Image(res,vol.affine,vol.header)
nb.save(output,"/data/pt_01880/test2.nii")

#%%