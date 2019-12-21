import os
import numpy as np
import nibabel as nb
import random
from lib_gbb.utils import get_gradient
from nibabel.freesurfer.io import read_geometry
from lib_gbb.utils import get_line
from lib_gbb.utils import update_mesh

input_geometry = "/home/daniel/projects/GBB/test_data/lh.layer10_def"
input_ref = "/home/daniel/projects/GBB/test_data/mean_data.nii"
input_vein = "/home/daniel/projects/GBB/test_data/vein.nii" 
path_output = "/home/daniel/Schreibtisch"
line_size = 2 # in mm
step_size = 0.4 # in mm
line_direction = 1 # line in F/H direction (phase-encoding direction)

""" do not edit below """

#%%

# get gradient (second order) of input volume
get_gradient(input_ref, path_output, input_vein, sigma_gaussian=1)

# load data
surf = read_geometry(input_geometry)
vol = nb.load(input_ref).get_fdata()
gradient = nb.load(os.path.join(path_output,"gradient.nii")).get_fdata()

# get line
n_coords = np.arange(0,len(surf[0]),1)
n_vertex = 83883 # random.choice(n_coords) 
vtx_coords = surf[0][n_vertex]
vox_coords = np.loadtxt(os.path.join(path_output,"vtx2vox.txt"))[n_vertex,:]

# get mean white matter response
vox_coords_all = np.loadtxt(os.path.join(path_output,"vtx2vox.txt")).astype(int)
mean_white = np.median(vol[vox_coords_all[:,0],vox_coords_all[:,1],vox_coords_all[:,2]])

# here is the get line function
line_coords, line_values = get_line(vox_coords, 
                                    gradient, 
                                    line_direction, 
                                    line_size, 
                                    step_size, 
                                    interpolation="linear", 
                                    show_plot=True)

# check shift
exit_loop = 0
i = 0
loc_new = []
while exit_loop == 0:
    if line_values[i+1] > 0 and line_values[i] < 0:
        loc_new = i
        break
    elif line_values[i+1] < 0 and line_values[i] > 0:
        loc_new = i
        break
    elif line_values[i] == 0:
        loc_new = i
        break
    
    if i >= len(line_values) - 2:
        exit_loop = 1
    else:
        i += 1

line_coords_res = line_coords[loc_new]
print(line_coords_res)

# get updated mesh
os.chdir(path_output)
update_mesh(input_geometry, n_vertex, vox_coords[1] - line_coords_res, n_iter=2, dir=2, sigmoid_range=5, write_output=True)
