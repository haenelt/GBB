"""
GBB

This script executes the gradient-based boundary (GBB) surface refinement.

created by Daniel Haenelt
Date created: 26-12-2019
Last modified: 25-02-2020
"""
import os
import shutil
import numpy as np
import nibabel as nb
from nibabel.freesurfer.io import write_geometry
from nibabel.freesurfer.io import read_geometry
from lib.surface.vox2ras import vox2ras
from lib_gbb.normal.get_normal_direction import get_normal_direction
from lib_gbb.neighbor.nn_3d import nn_3d
from lib_gbb.plot.cost_plot import cost_plot
from lib_gbb.plot.slope_plot import slope_plot
from lib_gbb.utils.devein_mesh import devein_mesh
from lib_gbb.utils.get_gradient import get_gradient
from lib_gbb.utils.get_shift import get_shift
from lib_gbb.utils.cost_BBR import cost_BBR
from lib_gbb.utils.update_mesh import update_mesh
from lib_gbb.utils.write_shift import write_shift
from lib_gbb.utils.deformation_field import deformation_field
from lib_gbb.utils.apply_shift import apply_shift
from lib_gbb.utils.check_exit import check_exit

# input files
input_white = "/home/daniel/projects/GBB/test_data/lh.layer10_def2_smooth"
input_pial = "/home/daniel/projects/GBB/test_data/lh.layer0_def2_smooth"
input_ref = "/home/daniel/projects/GBB/test_data/mean_epi_enhanced.nii"
input_vein = "/home/daniel/projects/GBB/test_data/vein.nii"
path_output = "/home/daniel/Schreibtisch/test"

# parameters
t2s = True # underlying image contrast (boolean)
line_dir = 2 # line axis in ras convention (0,1,2)
line_length = 3 # line length in one direction in mm
r_size = [5, 2.5, 1] # neighborhood radius in mm
l_rate = [0.1, 0.25, 0.5] # learning rate
max_iterations = [100000, 250000, 500000] # maximum iterations
cost_threshold = [1e-3, 1e-5, 1e-6] # cost function threshold

# gradient preparation
g_sigma = 1 # gaussian filter
g_kernel = 3 # kernel size used by gradient calculation

# deveining parameters
deveining = True # start with deveining the surface mesh (boolean)
n_neighbor_deveining = 20 # number of neighbors in surface relaxation
n_smooth_deveining = 0 # final smoothing
max_iterations_deveining = 1000 # maximum iterations

# cost parameters
cost_step = 1000 # step size between cost array points
cost_fit_min = 5 # sample size for linear fit

# deformation field
o_sigma = 1 # gaussian filter

# output
show_cost = True # show temporary cost function
show_slope = False # show temporary slope function
write_gradient = True # write gradient image
write_step = 10000 # step size to write intermediate surfaces (if set > 0)
cleanup = False # remove intermediate files

""" do not edit below """

# basename for output
name_white = os.path.basename(input_white)

if input_pial:
    name_pial = os.path.basename(input_pial)

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

# load volume data
vol_array = nb.load(input_ref).get_fdata()
grad_array = get_gradient(input_ref, ras2vox_tkr, line_dir, g_sigma, g_kernel, write_gradient,
                          path_output)
vein_array = nb.load(input_vein).get_fdata() 

# create textfile
file = open(os.path.join(path_output,name_white+"_info.txt"),"w")

# run deveining and load surface
if deveining:
    file.write("apply deveining\n")
    vtx_devein = devein_mesh(input_white, 
                             input_ref, 
                             input_vein, 
                             os.path.join(path_temp,name_white+"_devein"), 
                             n_neighbor_deveining, 
                             line_dir, 
                             n_smooth_deveining, 
                             max_iterations_deveining)
    
    # get deformation field
    print("write deformation field from deveining")
    vtx_old, _ = read_geometry(input_white)
    deformation_array = deformation_field(vtx_old, vtx_devein, input_ref, line_dir, 
                                          vox2ras_tkr, ras2vox_tkr, 
                                          o_sigma, 
                                          path_output=path_temp, 
                                          name_output=name_white+"_devein_deformation", 
                                          write_output=True)
    
    # get shift
    write_shift(vtx_devein, vtx_old, line_dir, 
                path_output=path_temp, 
                name_output=name_white+"_devein_shift")
      
    # apply deformation to white and pial surface
    print("apply deformation")
    vtx, fac = read_geometry(input_white)
    _ = apply_shift(vtx, fac, vox2ras_tkr, ras2vox_tkr, line_dir, deformation_array, 
                    path_output=path_temp, 
                    name_output=name_white+"_devein_refined", 
                    write_output=True)
    
    if input_pial:
        vtx, fac = read_geometry(input_pial)
        _ = apply_shift(vtx, fac, vox2ras_tkr, ras2vox_tkr, line_dir, deformation_array, 
                        path_output=path_temp, 
                        name_output=name_pial+"_devein_refined", 
                        write_output=True)
    
    # load deveined surface data
    vtx_old, fac_old = read_geometry(os.path.join(path_temp,name_white+"_devein_refined"))
else:
    # load surface data
    vtx_old, fac_old = read_geometry(input_white)

# get normals
n, vtx_norm = get_normal_direction(vtx_old, fac_old, line_dir)

print("start registration step 0 at iteration 0")
file.write("start registration step 0 at iteration 0\n")

# initialize some variables
vtx_new = vtx_old.copy()
n_coords = len(vtx_old)
n_steps = len(r_size)
vol_max = np.shape(vol_array)

# do surface refinement
i = 0
j = 0
p = 0
q = 0
counter = 0
step = 0
cost_array = []
m_array = []
n_array = [] 
while True:
            
    # get current vertex point
    n_vertex = np.random.randint(n_coords)
    
    # get shift
    vtx_shift = get_shift(vtx_new, fac_old, n, n_vertex, grad_array, vein_array, vox2ras_tkr, 
                          ras2vox_tkr, vol_max, line_length, line_dir, t2s, False)

    # update mesh
    if len(vtx_shift):
        nn_ind, _ = nn_3d(vtx_new[n_vertex], vtx_new, r_size[step])
        vtx_new = update_mesh(vtx_new,vtx_shift,n_vertex,nn_ind,l_rate[step])
    else:
        counter += 1
        continue
    
    # get cost function
    if p >= cost_step:
        J = cost_BBR(vtx_new, fac_old, vtx_norm, vol_array, ras2vox_tkr, vol_max, t2s)
        cost_array = np.append(cost_array, J)
        q += 1
        p = 0
        if len(cost_array) >= cost_fit_min:
            
            m_fit, n_fit, exit_crit = check_exit(np.arange(q-cost_fit_min,q), 
                                                 cost_array[-cost_fit_min:], 
                                                 cost_threshold=cost_threshold[step])
            m_array = np.append(m_array, m_fit)
            n_array = np.append(n_array, n_fit)
            
            if show_cost:
                set_title = "Exit criterion: "+str(exit_crit)+", Step: "+str(step)
                cost_plot(q, cost_array, m_fit, n_fit, set_title, save_plot=False, 
                          path_output=False, name_output=False)
                
            if show_slope:
                set_title = "Exit criterion: "+str(exit_crit)+", Step: "+str(step)
                slope_plot(q, m_array, set_title, save_plot=False, path_output=False, 
                           name_output=False)
    else:
        p += 1
        exit_crit = False

    # check exit
    if exit_crit and step < len(max_iterations) - 1 or j == max_iterations[step]:
        step += 1
        j = 0
        print("start registration step "+str(step)+" at iteration "+str(i))
        file.write("start registration step "+str(step)+" at iteration "+str(i)+"\n")
    elif exit_crit and step == len(max_iterations) - 1:
        print("Registration converged!")
        file.write("Registration converged!\n")        
        break
    elif step == len(max_iterations) - 1 and j == max_iterations[-1]:
        print("Registration did not converge!")
        file.write("Registration did not converge!\n")
        break

    # write intermediate surfaces
    if write_step and not np.mod(i,write_step):
        write_geometry(os.path.join(path_temp,"temp_"+str(i)), vtx_new, fac_old)
    
    i += 1
    j += 1

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
                                      name_output=name_white+"_deformation", 
                                      write_output=True)

# get shift
write_shift(vtx_new, vtx_old, line_dir, 
            path_output=path_output, 
            name_output=name_white+"_shift")

# apply deformation to white and pial surface
print("apply deformation")
vtx, fac = read_geometry(input_white)
_ = apply_shift(vtx, fac, vox2ras_tkr, ras2vox_tkr, line_dir, deformation_array, 
                path_output=path_output, 
                name_output=name_white+"_refined", 
                write_output=True)

if input_pial:
    vtx, fac = read_geometry(input_pial)
    _ = apply_shift(vtx, fac, vox2ras_tkr, ras2vox_tkr, line_dir, deformation_array, 
                    path_output=path_output, 
                    name_output=name_pial+"_refined", 
                    write_output=True)

# save cost array and slope and y-axis intercept arrays of linear fits
np.savez(os.path.join(path_output,name_white+"_cost"), J=cost_array, m=m_array, n=n_array)

# remove intermediate files
if cleanup:
    shutil.rmtree(path_temp, ignore_errors=True)