# -*- coding: utf-8 -*-


"""
List of parameters used in GBB

Parameters for deveining, anchoring, registration and the used cost function 
which are read by the main module. When the main module is run, all defined 
parameters are stored in a separate json file for later checks.

created by Daniel Haenelt
Date created: 03-10-2020
Last modified: 11-10-2020
"""

# io parameters
#-------------------------------------------------------------------------------

io_params = dict()
io_params["o_sigma"] = 2 # gaussian filter for deformation field

# deveining parameters
#-------------------------------------------------------------------------------

devein_params = dict()
devein_params["run"] = False # run deveining
devein_params["n_neighbor"] = 10 # number of neighbors in surface relaxation
devein_params["n_smooth"] = 0 # final smoothing
devein_params["max_iter"] = 1000 # maximum iterations

# anchoring parameters
#-------------------------------------------------------------------------------

anchor_params = dict()
anchor_params["run"] = False # run anchoring
anchor_params["n_neighbor"] = 10 # number of neighbors
anchor_params["n_smooth"] = 0 # final smoothing

# registration parameters
#-------------------------------------------------------------------------------

# NB: If the parameter line_dir is set to 3, the mesh is not shifted in
# direction of one of the 3 volume axes but in direction normal to the surface
# mesh.

# NB: Multiple neighborhood sizes (r_size) and learning rates (l_rate) can be 
# defined. During processing, the next neighborhood size and/or learning rate 
# is taken if the cost function is below the set threshold (cost_threshold) or 
# the maximum number of iterations (max_iter) is reached. To run properly, the 
# lists r_size, l_rate, max_iter and cost_threshold need to have the same number
# of elements.

gbb_params = dict()
gbb_params["run"] = True # run gbb
gbb_params["line_dir"] = 2 # line axis in ras convention (0,1,2,3)
gbb_params["line_length"] = 3 # line length in one direction in mm
gbb_params["r_size"] = [5, 2.5, 1] # neighborhood radius in mm
gbb_params["l_rate"] = [0.1, 0.1, 0.1] # learning rate
gbb_params["max_iter"] = [250000, 500000, 1000000] # maximum iterations
gbb_params["cost_threshold"] = [1e-4, 1e-5, 1e-6] # cost function threshold
gbb_params["gradient_sigma"] = 1 # gaussian blurring used by gradient calculation
gbb_params["gradient_kernel"] = 3 # kernel size used by gradient calculation
gbb_params["gradient_write"] = True # write gradient image
gbb_params["overwrite_control"] = False # do not lock control points
gbb_params["cost_step"] = 1000 # step size between cost array points
gbb_params["cost_sample"] = 10 # sample size for linear fit
gbb_params["show_cost"] = True # show temporary cost function
gbb_params["show_slope"] = False # show temporary slope function
gbb_params["intermediate_write"] = 10000 # step size to write intermediate surfaces (if set > 0)

# bbr parameters (cost function)
#-------------------------------------------------------------------------------

# NB: Cost function parameters are explained in the original BBR publication.
# Greve DN and Fischl B, Accurate and robust brain image alignment using 
# boundary-based registration, 48(1), 63--72 (2009).

bbr_params = dict()
bbr_params["t2s"] = True # underlying image contrast (boolean)
bbr_params["Q0"] = 0 # offset parameter in percent contrast measure
bbr_params["M"] = 0.5 # slope parameter in percent contrast measure
bbr_params["h"] = 1 # weight for each vertex in percent contrast measure
bbr_params["s"] = 1 # distance scaling for sampling normal to the surface
