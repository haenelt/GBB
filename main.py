"""
GBB

This script executes the gradient-based boundary (GBB) surface reginement.

created by Daniel Haenelt
Date created: 26-12-2019
Last modified: 26-12-2019
"""
import os
import random
import numpy as np
import nibabel as nb
import matplotlib.pyplot as plt
from nibabel.freesurfer.io import write_geometry
from nibabel.freesurfer.io import read_geometry
from lib_gbb.utils import get_gradient
from lib.surface.vox2ras import vox2ras
from lib_gbb.normal import get_normal_direction
from lib_gbb.utils.get_shift import get_shift
from lib_gbb.utils.cost_BBR import cost_BBR
from lib_gbb.neighbor.nn_3d import nn_3d
from lib_gbb.utils.update_mesh import update_mesh
from lib_gbb.utils.write_shift import write_shift

input_surf = "/home/daniel/projects/GBB/test_data/lh.layer10_def"
input_ref = "/home/daniel/projects/GBB/test_data/mean_data.nii"
input_vein = "/home/daniel/projects/GBB/test_data/vein.nii"
path_output = "/home/daniel/Schreibtisch/parameters13"
name_output = "lh.layer10_refined"

# parameters
t2s = True # underlying image contrast
line_dir = 2 # line axis in ras convention
line_length = 3 # line length in one direction in mm
r_size = [5, 2.5, 1] # neighborhood radius in mm
l_rate = [0.1, 0.1, 0.1] # learning rate
max_iterations = [50000, 50000, 50000] # maximum iterations
cost_threshold = [1e-8,1e-8,1e-8] # cost function threshold

# gradient preparation
sigma = 1
kernel_size = 3

# output
show_line = False
show_cost = True
write_gradient = True
write_intermediate = True
write_step = 1000

""" do not edit below """

# make output folder
if not os.path.exists(path_output):
    os.makedirs(path_output)

path_temp = os.path.join(path_output,"temp")
if not os.path.exists(path_temp):
    os.makedirs(path_temp)

# get transformation
vox2ras_tkr, ras2vox_tkr = vox2ras(input_ref)

# load data
vtx_old, fac_old = read_geometry(input_surf)
vol_array = nb.load(input_ref).get_fdata()
grad_array = get_gradient(input_ref, ras2vox_tkr, line_dir, sigma, kernel_size, write_gradient)
vein_array = nb.load(input_vein).get_fdata() 

# get gradient to output folder
if write_gradient:
    os.rename(os.path.join(os.path.dirname(input_ref),"gradient.nii"),
              os.path.join(path_temp,"gradient.nii"))

# initialize some variables
cost_array = []
n_coords = np.arange(0,len(vtx_old),1)
vtx_new = vtx_old.copy()

# get array size
vol_max = np.shape(vol_array)

# get normals
n, vtx_norm = get_normal_direction(vtx_old, fac_old, line_dir)

# create textfile
file = open(os.path.join(path_output,name_output+"_info.txt"),"w")

# do surface refinement
i = 0
counter = 0
step = 0
c_cost_count = 0
c_cost_size = 0
c_cost_check = 0
n_cost_count = 50
n_cost_size = 1000
n_cost_check = 5
max_iter = np.sum(max_iterations)
while i < max_iter:
    
    if i == 0:
        print("start coarse registration at iteration: "+str(i))
        file.write("start coarse registration at iteration: "+str(i)) 
    elif i == max_iterations[0] or c_cost_check == n_cost_check:
        print("start medium registration at iteration: "+str(i))
        file.write("start medium registration at iteration: "+str(i))
        step = 1
    elif i == max_iterations[0] + max_iterations[1] or c_cost_check == 2*n_cost_check:
        print("start fine registration at iteration: "+str(i))
        file.write("start fine registration at iteration: "+str(i))
        step = 2
       
    if not np.mod(i,n_cost_size):
        c_cost_count = 0
        
    # get current vertex point
    n_vertex = random.choice(n_coords)
    
    # get cost function
    if c_cost_count < n_cost_count:
        c_cost_size += cost_BBR(vtx_new, fac_old, vtx_norm, vol_array, ras2vox_tkr, vol_max, t2s)
    
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
    if c_cost_count  == n_cost_count:
        cost_array.append(c_cost_size / n_cost_count)
        c_cost_size = 0
        if len(cost_array) > 1:
            cost_diff = np.abs(cost_array[-1] - cost_array[-2])
            if cost_diff < cost_threshold[step]:
                c_cost_check += 1

        if show_cost:
            plt.plot(cost_array)
            plt.xlabel("iteration")
            plt.ylabel("Cost function")
            plt.pause(0.01)
    
    if c_cost_check == 3*n_cost_check and step == 2:
        print("registration converged!")
        break

    # write intermediate surfaces
    if not np.mod(i,write_step):
        write_geometry(os.path.join(path_temp,"temp_"+str(i)), vtx_new, fac_old)
    
    i += 1
    c_cost_count += 1

# close textfile
file.close()

# print some information
print("Final number of iterations: "+str(i))
print("Final number of skipped iterations: "+str(counter))
file.write("Final number of iterations: "+str(i))
file.write("Final number of skipped iterations: "+str(counter))

# write output surface and vertex shifts
write_geometry(os.path.join(path_output,name_output), vtx_new, fac_old)
write_shift(vtx_new, vtx_old, path_output, name_output)

# save cost array
np.save(os.path.join(path_output,name_output+"_cost"), cost_array)