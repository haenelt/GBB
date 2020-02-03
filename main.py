"""
GBB

This script executes the gradient-based boundary (GBB) surface reginement.

created by Daniel Haenelt
Date created: 26-12-2019
Last modified: 31-01-2020
"""
import os
import shutil
import numpy as np
import nibabel as nb
import matplotlib.pyplot as plt
from nibabel.freesurfer.io import write_geometry
from nibabel.freesurfer.io import read_geometry
from lib.surface.vox2ras import vox2ras
from lib_gbb.normal.get_normal_direction import get_normal_direction
from lib_gbb.neighbor.nn_3d import nn_3d
from lib_gbb.utils.get_gradient import get_gradient
from lib_gbb.utils.get_shift import get_shift
from lib_gbb.utils.cost_BBR import cost_BBR
from lib_gbb.utils.update_mesh import update_mesh
from lib_gbb.utils.write_shift import write_shift
from lib_gbb.utils.deformation_field import deformation_field
from lib_gbb.utils.apply_shift import apply_shift

# input files
input_white = "/home/raid2/haenelt/projects/GBB/test_data/surf/rigid_fmap/lh.layer10_def2_smooth_def_smooth"
input_pial = "/home/raid2/haenelt/projects/GBB/test_data/surf/rigid_fmap/lh.layer0_def2_smooth_def_smooth"
input_ref = "/home/raid2/haenelt/projects/GBB/test_data/data/mean_epi_enhanced.nii"
input_vein = "/home/raid2/haenelt/projects/GBB/test_data/data/vein.nii"
path_output = "/home/raid2/haenelt/gbb_test/refinement/epi_enhanced/rigid_fmap"

# parameters
t2s = True # underlying image contrast (boolean)
line_dir = 2 # line axis in ras convention (0,1,2)
line_length = 3 # line length in one direction in mm
r_size = [5, 2.5, 1] # neighborhood radius in mm
l_rate = [0.1, 0.1, 0.1] # learning rate
max_iterations = [100000, 250000, 500000] # maximum iterations
cost_threshold = [1e-3, 5e-4, 1e-4] # cost function threshold
cleanup = False

# gradient preparation
g_sigma = 1 # gaussian filter
g_kernel = 3 # kernel size used by gradient calculation

# deformation field
o_sigma = 1 # gaussian filter

# output
show_cost = True # show temporary cost function
write_gradient = True # write gradient image
write_step = 10000 # step size to write intermediate surfaces (if set > 0)

""" do not edit below """

# basename for output
name_output = os.path.basename(input_white)

# make output folder
if not os.path.exists(path_output):
    os.makedirs(path_output)

tmp = np.random.randint(0, 10, 5)
tmp_string = ''.join(str(i) for i in tmp)
path_temp = os.path.join(path_output,"tmp_"+tmp_string)
if not os.path.exists(path_temp):
    os.makedirs(path_temp)

# get transformation
vox2ras_tkr, ras2vox_tkr = vox2ras(input_ref)

# load data
vtx_old, fac_old = read_geometry(input_white)
vol_array = nb.load(input_ref).get_fdata()
grad_array = get_gradient(input_ref, ras2vox_tkr, line_dir, g_sigma, g_kernel, write_gradient)
vein_array = nb.load(input_vein).get_fdata() 

# get gradient to output folder
if write_gradient:
    os.rename(os.path.join(os.path.dirname(input_ref),"gradient.nii"),
              os.path.join(path_temp,"gradient.nii"))

# get normals
n, vtx_norm = get_normal_direction(vtx_old, fac_old, line_dir)

# create textfile
file = open(os.path.join(path_output,name_output+"_info.txt"),"w")

print("start registration step 0 at iteration 0")
file.write("start registration step 0 at iteration 0\n")

# initialize some variables
vtx_new = vtx_old.copy()
n_coords = len(vtx_old)
max_iter = np.sum(max_iterations)
n_steps = len(r_size)
vol_max = np.shape(vol_array)

# cost function parameters
n_cost_count = 50 # number of time points for cost array average
n_cost_size = 1000 # step size between cost array points
n_cost_check = 5 # number successive time points below threshold for convergent registration

# do surface refinement
i = 0
counter = 0
step = 0
cost_array = []
c_steps = 0
c_cost_count = 0
c_cost_size = 0
c_cost_check = 0
while i < max_iter:
    
    if i == max_iterations[step] + c_steps or c_cost_check == n_cost_check:
        step += 1
        c_cost_check = 0
        c_steps += i
        print("start registration step "+str(step)+" at iteration "+str(i))
        file.write("start registration step "+str(step)+" at iteration "+str(i)+"\n")
       
    if not np.mod(i,n_cost_size):
        c_cost_count = 0
        
    # get current vertex point
    n_vertex = np.random.randint(n_coords)
    
    # get cost function
    if c_cost_count < n_cost_count:
        c_cost_size += cost_BBR(vtx_new, fac_old, vtx_norm, vol_array, ras2vox_tkr, vol_max, t2s)
    
    # get shift
    vtx_shift = get_shift(vtx_new, fac_old, n, n_vertex, grad_array, vein_array, vox2ras_tkr, 
                          ras2vox_tkr, vol_max, line_length, line_dir, t2s, False)

    # update mesh
    if len(vtx_shift):
        nn_ind = nn_3d(vtx_new[n_vertex], vtx_new, r_size[step])
        vtx_new = update_mesh(vtx_new,vtx_shift,n_vertex,nn_ind,l_rate[step])
    else:
        counter += 1
    
    # update exit criterion
    if c_cost_count  == n_cost_count:
        cost_array.append(c_cost_size / n_cost_count)
        c_cost_size = 0
        if len(cost_array) > 1:
            cost_diff = np.abs(cost_array[-1] - cost_array[-2])
            if c_cost_check > 0:
                cost_diff2 = np.abs(cost_array[-2] - cost_array[-3])
                if cost_diff < cost_threshold[step] and cost_diff2 < cost_threshold[step]:
                    c_cost_check += 1
                else:
                    c_cost_check = 0
            else:
                if cost_diff < cost_threshold[step]:
                    c_cost_check = 1

        if show_cost:
            plt.plot(cost_array)
            plt.xlabel("iteration")
            plt.ylabel("Cost function")
            plt.pause(0.01)
    
    # check exit criterion
    if c_cost_check == n_cost_check and step == n_steps - 1:
        print("Registration converged!")
        file.write("Registration converged!\n")
        break
    elif i == max_iterations[step] and step == n_steps - 1:
        print("Registration did not converge!")
        file.write("Registration did not converge!\n")
        break

    # write intermediate surfaces
    if write_step > 0 and not np.mod(i,write_step):
        write_geometry(os.path.join(path_temp,"temp_"+str(i)), vtx_new, fac_old)
    
    i += 1
    c_cost_count += 1

# print some information
print("Final number of iterations: "+str(i))
print("Final number of skipped iterations: "+str(counter))
file.write("Final number of iterations: "+str(i)+"\n")
file.write("Final number of skipped iterations: "+str(counter)+"\n")

# close textfile
file.close()

# get deformation field
print("write deformation field")
deformation_array = deformation_field(vtx_old, vtx_new, input_ref, line_dir, 
                                      vox2ras_tkr, ras2vox_tkr, 
                                      o_sigma, 
                                      path_output=path_output, 
                                      name_output=name_output+"_deformation", 
                                      write_output=True)

# get shift
write_shift(vtx_new, vtx_old, line_dir, 
            path_output=path_output, 
            name_output=name_output+"_shift")

# apply deformation to white and pial surface
print("apply deformation")
vtx, fac = read_geometry(input_white)
_ = apply_shift(vtx, fac, vox2ras_tkr, ras2vox_tkr, line_dir, deformation_array, 
                path_output=path_output, 
                name_output=os.path.basename(input_white)+"_refined", 
                write_output=True)

vtx, fac = read_geometry(input_pial)
_ = apply_shift(vtx, fac, vox2ras_tkr, ras2vox_tkr, line_dir, deformation_array, 
                path_output=path_output, 
                name_output=os.path.basename(input_pial)+"_refined", 
                write_output=True)

# save cost array
np.save(os.path.join(path_output,name_output+"_cost"), cost_array)

# remove intermediate files
if cleanup:
    shutil.rmtree(path_temp, ignore_errors=True)
