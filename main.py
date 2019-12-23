import os
import random
import numpy as np
import nibabel as nb
import matplotlib.pyplot as plt
from nibabel.freesurfer.io import write_geometry
from lib_gbb.utils import get_gradient
from nibabel.freesurfer.io import read_geometry
from lib.surface.vox2ras import vox2ras
from lib_gbb.normal import get_normal_direction
from lib_gbb.utils.get_adjm import get_adjm
from lib_gbb.utils.get_shift import get_shift
from lib_gbb.utils.cost_BBR import cost_BBR
from lib_gbb.neighbor.nn_3d import nn_3d
from lib_gbb.utils.update_mesh import update_mesh
from lib_gbb.utils.write_shift import write_shift

input_surf = "/home/daniel/projects/GBB/test_data/lh.layer10_def"
input_ref = "/home/daniel/projects/GBB/test_data/mean_data.nii"
input_vein = "/home/daniel/projects/GBB/test_data/vein.nii"
path_output = "/home/daniel/Schreibtisch/parameters1"
name_output = "lh.layer10_refined"

# parameters
t2s = True # underlying image contrast
line_dir = 2 # line direction in ras convention
line_length = 3 # line length in one direction in mm
cost_thres = 1e-8 # cost function threshold
coarse_iterations = 20000
medium_iterations = 20000
r_size = [20, 10, 5] # neighborhood
l_rate = [0.25, 0.15, 0.1] # learning rate
max_iter = 100000

# gradient preparation
sigma = 1
kernel_size = 3

# output
show_line = False
show_cost = True
write_intermediate = True
write_step = 10000

""" do not edit below """

# make output folder
if not os.path.exists(path_output):
    os.makedirs(path_output)

path_temp = os.path.join(path_output,"temp")
if not os.path.exists(path_temp):
    os.makedirs(path_temp)

# get transformation
print("get transformations")
vox2ras_tkr, ras2vox_tkr = vox2ras(input_ref)

print("load data")
vtx_old, fac_old = read_geometry(input_surf)
vol_array = nb.load(input_ref).get_fdata()
grad_array = get_gradient(input_ref, ras2vox_tkr, line_dir, sigma, kernel_size, write_output=None)
vein_array = nb.load(input_vein).get_fdata() 

# initialize some variables
cost_array = []
n_coords = np.arange(0,len(vtx_old),1)
vtx_new = vtx_old.copy()

# get array size
vol_max = np.shape(vol_array)

print("compute adjacency matrix")
adjm = get_adjm(input_surf)

# get normals
print("compute normals")
n = get_normal_direction(vtx_old, fac_old, line_dir)

print("do surface refinement")
i = 0
counter = 0
step = 0
while i < max_iter:
    
    if i == 0:
        print("start coarse iterations")
    elif i == coarse_iterations:
        print("start medium iterations")
        step = 1
    elif i == coarse_iterations + medium_iterations:
        print("start fine iterations")
        step = 2
        

    # get current vertex point
    n_vertex = random.choice(n_coords)
    
    # get cost function
    cost_array.append(cost_BBR(vtx_new, fac_old, vol_array, n, ras2vox_tkr, vol_max, 2, line_dir))
    
    # get shift
    vtx_shift = get_shift(vtx_new, fac_old, n, n_vertex, grad_array, vein_array, vox2ras_tkr, 
                          ras2vox_tkr, vol_max, line_length, line_dir, t2s, show_line)

    # update mesh
    if len(vtx_shift):
        nn_ind = nn_3d(vtx_new[n_vertex], vtx_new, r_size[step])
        vtx_new = update_mesh(vtx_new,vtx_shift,n_vertex,nn_ind,l_rate[step])
    else:
        counter += 1
    
    # exist criterion
    if i > 0 and not np.mod(i,1000):
        cost_old = np.mean(cost_array[-1050:-1000])
        cost_new = np.mean(cost_array[-50:])
        if np.abs(cost_new - cost_old) < cost_thres:
            break

    if show_cost and i > 0 and not np.mod(i,write_step):
        plt.plot(cost_array)
        plt.xlabel("iteration")
        plt.ylabel("Cost function")
        plt.pause(0.01)

    # write intermediate surfaces
    if not np.mod(i,write_step):
        write_geometry(os.path.join(path_temp,"temp_"+str(i)), vtx_new, fac_old)
    
    i += 1

# print some information
print("Number of iterations: "+str(i))
print("Number of skipped iterations "+str(counter))

# write output surface and vertex shifts
write_geometry(os.path.join(path_output,name_output), vtx_new, fac_old)
write_shift(vtx_new, vtx_old, path_output, name_output)