# -*- coding: utf-8 -*-


"""
List of parameters used in GBB

Parameters for deveining, anchoring and registration which are read by the main
module. When the main module is run, all parameters are stored in a json file
for later checks.

created by Daniel Haenelt
Date created: 03-10-2020
Last modified: 03-10-2020
"""

# io parameters
io_params = dict()
io_params["o_sigma"] = 2 # gaussian filter for deformation field

# deveining parameters
devein_params = dict()
devein_params["run"] = False # run deveining
devein_params["n_neighbor"] = 10 # number of neighbors in surface relaxation
devein_params["n_smooth"] = 0 # final smoothing
devein_params["max_iter"] = 1000 # maximum iterations

# anchoring parameters
anchor_params = dict()
anchor_params["run"] = False # run anchoring
anchor_params["n_neighbor"] = 10 # number of neighbors
anchor_params["n_smooth"] = 0 # final smoothing

# registration parameter
reg_params = dict()
reg_params["run"] = True # run gbb
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
