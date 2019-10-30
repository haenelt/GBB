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
from scipy.interpolate import splev, splrep
from lib_gbb.interpolation import linear_interpolation2d, nn_interpolation2d
import matplotlib.pyplot as plt
    
input_geometry = "/home/daniel/projects/GBB/test_data/lh.layer10_def"
input_vol = "/home/daniel/projects/GBB/test_data/mean_data.nii"
input_gradient = "/home/daniel/projects/GBB/test_data/gradient.nii"
input_vtx2vox = "/home/daniel/Schreibtisch/test/vtx2vox.txt"
line_size = 2 # in mm
step_size = 0.5 # in mm
direction = 1
interpolation = "linear"
    
surf = read_geometry(input_geometry)
vol = nb.load(input_vol).get_fdata()
gradient = nb.load(input_gradient).get_fdata()
 
n_coords = np.arange(0,len(surf[0]),1)
n_vertex = random.choice(n_coords) 
vtx_coords = surf[0][n_vertex]
vox_coords = np.loadtxt(input_vtx2vox)[n_vertex,:]

line_coords = np.arange(vox_coords[direction] - line_size, vox_coords[direction] + line_size, step_size) + 1

# exclude coordinates at the edges
line_coords = line_coords[line_coords > 0]
line_coords = line_coords[line_coords < np.shape(vol)[direction]] -1 

# get mean white matter response
vox_coords_all = np.loadtxt(input_vtx2vox).astype(int)
mean_white = np.median(vol[vox_coords_all[:,0],vox_coords_all[:,1],vox_coords_all[:,2]])
    
# get line
line = []
x = vox_coords[0]
z = np.round(vox_coords[2]).astype(int)
for i in range(len(line_coords)):
        
    y = line_coords[i]
    if interpolation == "linear":
        line.append(linear_interpolation2d(x, y, gradient[:,:,z]))
    elif interpolation == "nearest":
        line.append(nn_interpolation2d(x, y, gradient[:,:,z]))
    else:
        print("choose a valid interpolation method!")        

line_coords_new = np.linspace(line_coords.min(),line_coords.max(),1000) #300 represents number of points to make between T.min and T.max
line_spl = splrep(line_coords, line)
line_smooth = splev(line_coords_new, line_spl)
plt.plot(line_coords, line, 'o', line_coords_new, line_smooth)
plt.show()

# decide if in gm or wm
if vol[int(vox_coords[0]),int(vox_coords[1]),int(vox_coords[2])] > mean_white: 
    test_loc = line_coords_new[int(np.where(line_smooth==np.min(line_smooth))[0][0])] # in gm
else:
    test_loc = line_coords_new[int(np.where(line_smooth==np.max(line_smooth))[0][0])] # in wm

print(test_loc)