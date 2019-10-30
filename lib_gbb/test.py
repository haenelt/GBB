from lib_gbb.utils import map_cmap

input_surf = "/home/raid2/haenelt/projects/GBB/test_data/lh.layer10_def"
input_ref = "/home/raid2/haenelt/projects/GBB/test_data/mean_data.nii"
hemi = "lh"
path_output = "/data/pt_01880/test"

map_cmap(input_surf, input_ref, hemi, path_output)

#%%

# get line
import numpy as np
import nibabel as nb
from nibabel.freesurfer.io import read_geometry
from lib_gbb.interpolation import linear_interpolation2d
import matplotlib.pyplot as plt
    
input_geometry = "/home/raid2/haenelt/projects/GBB/test_data/lh.layer10_def"
input_vol = "/home/raid2/haenelt/projects/GBB/test_data/mean_data.nii"
input_vtx2vox = "/data/pt_01880/test/vtx2vox.txt"
n_vertex = 5755
line_size = 1 # in mm
step_size = 0.5 # in mm
direction = 1
interpolation = "linear"
    
surf = read_geometry(input_geometry)
vol = nb.load(input_vol).get_fdata()
    
vtx_coords = surf[0][n_vertex]
vox_coords = np.loadtxt(input_vtx2vox)[n_vertex,:]

line_coords = np.arange(vox_coords[direction] - line_size, vox_coords[direction] + line_size, step_size) + 1

# exclude coordinates at the edges
line_coords = line_coords[line_coords > 0]
line_coords = line_coords[line_coords < np.shape(vol)[direction]] -1 

# get mean white matter response
vox_coords_all = np.loadtxt(input_vtx2vox).astype(int)
mean_white = np.median(vol[vox_coords_all[:,0],vox_coords_all[:,1],vox_coords_all[:,2]])
    
#%%
   
if interpolation == "linear":
    
    res = []
    for i in range(len(line_coords)):
        
        y = vox_coords[0]
        x = line_coords[i]
        x_s = [np.floor(x), np.ceil(x)]
        y_s = [np.floor(y), np.ceil(y)]
        f_s = [[vol[int(y_s[0]),int(x_s[0]),int(np.round(vox_coords[2]))],vol[int(y_s[1]),int(x_s[0]),int(np.round(vox_coords[2]))]],
                    [vol[int(y_s[0]),int(x_s[1]),int(np.round(vox_coords[2]))],vol[int(y_s[1]),int(x_s[1]),int(np.round(vox_coords[2]))]]]
        
        res.append(linear_interpolation2d(x, y, x_s, y_s, f_s))

elif interpolation == "nearest":
    
    res = []
    res_y = []
    for i in range(len(line_coords)):
        
        x = int(vox_coords[0])
        y = int(line_coords[i])
        f_s = vol[x,y,np.round(vox_coords[2]).astype(int)]
        
        res.append(f_s)
        res_y.append(y)

else:
    print("choose a valid interpolation method!")


plt.plot(line_coords,res)

#%%


from scipy.interpolate import spline

test = [0]
for i in range(len(res) - 1):
    test.append(res[i+1] - res[i])



line_coords_new = np.linspace(line_coords.min(),line_coords.max(),1000) #300 represents number of points to make between T.min and T.max
test_smooth = spline(line_coords,test,line_coords_new)

plt.plot(line_coords_new,test_smooth)

#%%

# decide if in gm or wm
if vol[int(vox_coords[0]),int(vox_coords[1]),int(vox_coords[2])] > mean_white: 
    test_loc = line_coords_new[int(np.where(test_smooth==np.min(test_smooth))[0][0])] # in gm
else:
    test_loc = line_coords_new[int(np.where(test_smooth==np.max(test_smooth))[0][0])] # in wm

print(test_loc)