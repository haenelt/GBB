import random
import numpy as np
import nibabel as nb
from lib_gbb.utils import get_gradient
from nibabel.freesurfer.io import read_geometry
from lib_gbb.normal import get_normal_direction
from lib_gbb.utils.get_adjm import get_adjm

input_surf = "/home/daniel/projects/GBB/test_data/lh.layer10_def"
input_ref = "/home/daniel/projects/GBB/test_data/mean_data.nii"
input_vein = "/home/daniel/projects/GBB/test_data/vein.nii"

""" do not edit below """

print("compute adjacency matrix")
adjm = get_adjm(input_surf)

#%% get neighbour

nn = None
ind = 0
n_iter = 10


i = 0
while i < n_iter:
    
    if i == 0:
        nn = adjm[ind,:].indices
    else:
        nn_temp = nn.copy()
        x = [np.append(nn,adjm[nn_temp[j],:].indices) for j in range(len(nn_temp))]
    

#%%

# get current vertex point
n_coords = np.arange(0,len(vtx),1)
n_vertex = random.choice(n_coords)

#%%

ind = 27000
mm_line = 3 # in mm
n_line = 1000
show_plot = True
delta_dir = 2
vol_xmax = 100
vol_ymax = 100
vol_zmax = 100

# get transformation
vox2ras, ras2vox = vox2ras(input_ref)

# get vol and gradient
vol_array = nb.load(input_ref).get_fdata()
grad_array = get_gradient(input_ref, input_vein, sigma_gaussian=1, grad_dir=1, kernel_size=3)

# get surf
vtx, fac = read_geometry(input_geometry)

# get normals
n = get_normal_direction(vtx, fac, diff_dir