import random
import numpy as np
import nibabel as nb
from nibabel.freesurfer.io import write_geometry
from lib_gbb.utils import get_gradient
from nibabel.freesurfer.io import read_geometry
from lib_gbb.normal import get_normal_direction
from lib_gbb.utils.get_adjm import get_adjm
from lib.surface.vox2ras import vox2ras
from lib_gbb.utils.get_shift import get_shift
from lib_gbb.utils.cost_BBR import cost_BBR
from lib_gbb.neighbor.nn_3d import nn_3d
from lib_gbb.utils.update_mesh import update_mesh

input_surf = "/home/daniel/projects/GBB/test_data/lh.layer10_def"
input_ref = "/home/daniel/projects/GBB/test_data/mean_data.nii"
input_vein = "/home/daniel/projects/GBB/test_data/vein.nii"

# parameters
delta_dir = 2
mm_line = 5 # in mm
n_line = 1000
show_plot = True
t2s = True
r_size = 10 # in mm

""" do not edit below """

print("load data")
vtx, fac = read_geometry(input_surf)
vol_array = nb.load(input_ref).get_fdata()
grad_array = get_gradient(input_ref, input_vein, sigma_gaussian=1, grad_dir=1, kernel_size=3)

# get array size
vol_max = np.shape(vol_array)

print("compute adjacency matrix")
adjm = get_adjm(input_surf)

# get transformation
print("get transformations")
vox2ras_tkr, ras2vox_tkr = vox2ras(input_ref)

# get normals
print("compute normals")
n = get_normal_direction(vtx, fac, delta_dir)

# get cost function
cost_array = []
cost_check = np.zeros(50)
cost_val = cost_BBR(vtx, fac, vol_array, n, ras2vox_tkr, vol_max, delta_size=2, delta_dir=2)
cost_array.append(cost_val)
cost_check = np.roll(cost_check, 1)
cost_check[0] = cost_val

#%%

# get current vertex point
n_coords = np.arange(0,len(vtx),1)
n_vertex = 154121#random.choice(n_coords)

# get shift
vtx_shift = get_shift(vtx, fac, n, n_vertex, grad_array, vox2ras_tkr, ras2vox_tkr, vol_max, mm_line, 
                      n_line, delta_dir, t2s, show_plot=True)

# get neighborhood
nn_ind = nn_3d(vtx[n_vertex,:], vtx, r_size)

#%% update mesh
if shift_curr:
    vtx = update_mesh(vtx,vtx_shift,n_vertex,nn_ind,l_rate=1.0)

# exit criterion
# mache loop
# dynamic cost function plot
# write out intermediate surfaces
# parameter changes: r_size (neighborhood), scale factor (learning rate)
# final output: (1) surface, (2) final shifts, (3) accuracy to GM/WM border (metric for each vertex convergence)
# count iterations, count number of Nones










