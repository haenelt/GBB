import random
import numpy as np
import nibabel as nb
from lib_gbb.utils import get_gradient
from nibabel.freesurfer.io import read_geometry
from lib_gbb.normal import get_normal_direction
from lib_gbb.utils.get_adjm import get_adjm
from lib.surface.vox2ras import vox2ras
from lib_gbb.utils.get_shift import get_shift

input_surf = "/home/daniel/projects/GBB/test_data/lh.layer10_def"
input_ref = "/home/daniel/projects/GBB/test_data/mean_data.nii"
input_vein = "/home/daniel/projects/GBB/test_data/vein.nii"

# parameters
delta_dir = 2

""" do not edit below """

print("load data")
vtx, fac = read_geometry(input_surf)
vol_array = nb.load(input_ref).get_fdata()
grad_array = get_gradient(input_ref, input_vein, sigma_gaussian=1, grad_dir=1, kernel_size=3)

print("compute adjacency matrix")
adjm = get_adjm(input_surf)

# get transformation
print("get transformations")
vox2ras_tkr, ras2vox_tkr = vox2ras(input_ref)

# get normals
print("compute normals")
n = get_normal_direction(vtx, fac, delta_dir)

# 

#%%


# get current vertex point
n_coords = np.arange(0,len(vtx),1)
n_vertex = random.choice(n_coords)

#%%

ind = 27000
mm_line = 3 # in mm
n_line = 1000
show_plot = True
vol_xmax = 100
vol_ymax = 100
vol_zmax = 100