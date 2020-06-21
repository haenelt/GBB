"""
GBB

This script executes the gradient-based boundary (GBB) surface refinement.

created by Daniel Haenelt
Date created: 26-12-2019
Last modified: 18-05-2020
"""
import os
import numpy as np
from lib_gbb.io.load_data import load_data
from lib_gbb.io.write_readme import write_readme
from lib_gbb.io.write_shift import write_shift
from lib_gbb.methods.devein_mesh import devein_mesh
from lib_gbb.methods.anchor_mesh import anchor_mesh
from lib_gbb.methods.gbb_mesh import gbb_mesh
from lib_gbb.utils.deformation_field import deformation_field
from lib_gbb.utils.apply_deformation import apply_deformation

# input and output parameters
io_file = dict()
io_file["i_white"] = "/data/pt_01880/Experiment1_ODC/p1/anatomy/dense_deformed/rh.white_def2"
io_file["i_pial"] = "/data/pt_01880/Experiment1_ODC/p1/anatomy/dense_deformed/rh.pial_def2"
io_file["i_ref"] = "/data/pt_01880/Experiment1_ODC/p1/resting_state/mean_udata_enhanced.nii"
io_file["i_vein"] = "/data/pt_01880/Experiment1_ODC/p1/resting_state/vein/native/vein.nii"
io_file["i_ignore"] = None
io_file["i_anchor"] = None
io_file["o_output"] = "/data/pt_01880/Experiment1_ODC/p1/anatomy/dense_refined/gbb"
io_file["o_sigma"] = 2 # gaussian filter for deformation field

# deveining parameters
devein_params = dict()
devein_params["run"] = False
devein_params["n_neighbor"] = 10 # number of neighbors in surface relaxation
devein_params["n_smooth"] = 0 # final smoothing
devein_params["max_iter"] = 1000 # maximum iterations

# anchoring parameters
anchor_params = dict()
anchor_params["run"] = False
anchor_params["n_neighbor"] = 20 # number of neighbors
anchor_params["n_smooth"] = 0 # final smoothing

# registration parameter
reg_params = dict()
reg_params["run"] = True
reg_params["t2s"] = True # underlying image contrast (boolean)
reg_params["line_dir"] = 2 # line axis in ras convention (0,1,2,3)
reg_params["line_length"] = 3 # line length in one direction in mm
reg_params["r_size"] = [5, 2.5, 1] # neighborhood radius in mm
reg_params["l_rate"] = [0.1, 0.1, 0.1] # learning rate
reg_params["max_iter"] = [250000, 500000, 1000000] # maximum iterations
reg_params["cost_threshold"] = [1e-4, 1e-5, 1e-6] # cost function threshold
reg_params["gradient_sigma"] = 1 # gaussian blurring used by gradient calculation
reg_params["gradient_kernel"] = 3 # kernel size used by gradient calculation
reg_params["gradient_write"] = True # write gradient image
reg_params["overwrite_control"] = False # do not lock control points
reg_params["cost_step"] = 1000 # step size between cost array points
reg_params["cost_sample"] = 10 # sample size for linear fit
reg_params["show_cost"] = True # show temporary cost function
reg_params["show_slope"] = False # show temporary slope function
reg_params["intermediate_write"] = 10000 # step size to write intermediate surfaces (if set > 0)

""" do not edit below """

# make output folder
if not os.path.exists(io_file["o_output"]):
    os.makedirs(io_file["o_output"])

# load input
volume, T, surf, point, basename = load_data(io_file, reg_params)

# run deveining
niter_devein = 0
if devein_params["run"]:        
    surf["vtx_white"], niter_devein = devein_mesh(surf["vtx_white"], 
                                                  surf["fac_white"], 
                                                  surf["n_white"],
                                                  volume["vein"], 
                                                  volume["ignore"],  
                                                  T["adjm"], 
                                                  T["ras2vox"],
                                                  devein_params["n_neighbor"], 
                                                  3, # along all directions
                                                  devein_params["n_smooth"], 
                                                  devein_params["max_iter"])   

# run anchoring
ind_control = []
if anchor_params["run"]:
    surf["vtx_white"], ind_control = anchor_mesh(surf["vtx_white"], 
                                                 surf["fac_white"], 
                                                 T["adjm"], 
                                                 point["anchor"],
                                                 anchor_params["n_neighbor"], 
                                                 anchor_params["n_smooth"])

# run gbb
gbb_params = None
if reg_params["run"]:
    
    # do not consider control points in gbb
    if reg_params["overwrite_control"]:
        ind_control = []
    
    surf["vtx_white"], gbb_params = gbb_mesh(surf["vtx_white"],
                                             surf["fac_white"], 
                                             surf["n_white"],
                                             ind_control, 
                                             volume["ref"], 
                                             volume["gradient"], 
                                             volume["vein"], 
                                             volume["ignore"], 
                                             reg_params["t2s"], 
                                             T["vox2ras"], 
                                             T["ras2vox"], 
                                             reg_params["line_dir"], 
                                             reg_params["line_length"], 
                                             reg_params["r_size"], 
                                             reg_params["l_rate"], 
                                             reg_params["max_iter"], 
                                             reg_params["cost_threshold"], 
                                             reg_params["cost_step"], 
                                             reg_params["cost_sample"], 
                                             io_file["o_output"], 
                                             reg_params["show_cost"], 
                                             reg_params["show_slope"], 
                                             reg_params["intermediate_write"])
    
    # save cost array and slope and y-axis intercept arrays of linear fits
    np.savez(os.path.join(io_file["o_output"],basename["white"]+"_cost"), 
             J= gbb_params["cost_array"], m=gbb_params["m_array"], n=gbb_params["n_array"])

# write deformation field
print("write deformation field")
arr_deformation = deformation_field(surf["vtx_white_archive"],
                                    surf["vtx_white"],
                                    io_file["i_ref"],
                                    T["vox2ras"],
                                    T["ras2vox"],
                                    io_file["o_sigma"],
                                    io_file["o_output"],
                                    basename["white"]+"_deformation",
                                    True)

# get shift
print("write shift")
write_shift(surf["vtx_white_archive"], 
            surf["vtx_white"],
            reg_params["line_dir"], 
            io_file["o_output"], 
            basename["white"]+"_shift")
  
# apply deformation to white and pial surface
print("apply deformation")
_, _ = apply_deformation(surf["vtx_white_archive"],
                         surf["fac_white"],
                         arr_deformation,
                         T["vox2ras"],
                         T["ras2vox"],
                         io_file["o_output"],
                         basename["white"]+"_refined", 
                         True)

if surf["vtx_pial"] is not None:
    _, _ = apply_deformation(surf["vtx_pial"], 
                             surf["fac_pial"], 
                             arr_deformation,
                             T["vox2ras"],
                             T["ras2vox"],
                             io_file["o_output"], 
                             basename["pial"]+"_refined", 
                             True)

# write readme
write_readme(devein_params, anchor_params, reg_params, gbb_params, io_file["o_sigma"], niter_devein, 
             len(ind_control), io_file["o_output"], basename["white"])
