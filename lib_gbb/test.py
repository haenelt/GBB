from lib_gbb.utils import map_cmap

input_surf = "/home/daniel/projects/GBB/test_data/lh.layer10_def"
input_ref = "/home/daniel/projects/GBB/test_data/mean_data.nii"
hemi = "lh"
path_output = "/home/daniel/Schreibtisch/test"

map_cmap(input_surf, input_ref, hemi, path_output)

#%%
from lib_gbb.utils import get_gradient

input_vol = "/home/daniel/projects/GBB/test_data/mean_data.nii"
input_vein = "/home/daniel/projects/GBB/test_data/vein.nii" 
path_output = "/home/daniel/projects/GBB/test_data"

get_gradient(input_vol, path_output, input_vein=input_vein)

#%%

# get line
import numpy as np
import nibabel as nb
import random
from nibabel.freesurfer.io import read_geometry
from lib_gbb.utils import get_line
    
input_geometry = "/home/daniel/projects/GBB/test_data/lh.layer10_def"
input_vol = "/home/daniel/projects/GBB/test_data/mean_data.nii"
input_gradient = "/home/daniel/projects/GBB/test_data/gradient.nii"
input_vtx2vox = "/home/daniel/Schreibtisch/test/vtx2vox.txt"
line_size = 2 # in mm
step_size = 0.4 # in mm
trans_white = 100
edge_threshold = 5
direction = 1
interpolation = "linear"
    
surf = read_geometry(input_geometry)
vol = nb.load(input_vol).get_fdata()
gradient = nb.load(input_gradient).get_fdata()
 
n_coords = np.arange(0,len(surf[0]),1)
n_vertex = random.choice(n_coords) 
vtx_coords = surf[0][n_vertex]
vox_coords = np.loadtxt(input_vtx2vox)[n_vertex,:]

# get mean white matter response
vox_coords_all = np.loadtxt(input_vtx2vox).astype(int)
mean_white = np.median(vol[vox_coords_all[:,0],vox_coords_all[:,1],vox_coords_all[:,2]])

# here is the get line function
line_coords, line_values = get_line(vox_coords, 
                                    gradient, 
                                    direction, 
                                    line_size, 
                                    step_size, 
                                    interpolation="linear", 
                                    show_plot=True)

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