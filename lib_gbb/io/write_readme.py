def write_readme(devein_params, anchor_params, reg_params, gbb_params, niter_d, niter_a, path_output, 
                 name_output):
    """
    This function writes a textfile summarizing input and output parameters. 
    Inputs.
        *devein_params: deveining parameters (dict).
        *anchor_params: anchoring parameters (dict).
        *reg_params: registration parameters (dict).
        *gbb_params: gbb performance parameters (dict).
        *niter_d: number of deveining iterations.
        *niter_a: number of anchoring iterations.
        *path_output: path where output is saved.
        *name_output: basename of output file.

    created by Daniel Haenelt
    Date created: 12-05-2020 
    Last modified: 13-05-2020  
    """
    import os
    
    # create textfile
    file = open(os.path.join(path_output,name_output+"_info.txt"),"w")
    
    # deveining
    file.write("apply deveining: "+str(devein_params["run"])+"\n")
    if devein_params["run"]:
        
        # deveining parameters
        file.write("Number of neighbors: "+str(devein_params["n_neighbor"])+"\n")
        file.write("Final smoothing iterations: "+str(devein_params["n_smooth"])+"\n")
        file.write("Maximum iterations: "+str(devein_params["max_iter"])+"\n")
        
        # number of deveining iterations
        file.write("Number of deveining iterations: "+str(niter_d)+"\n")
        
    # anchoring
    file.write("apply anchoring: "+str(anchor_params["run"]))
    if anchor_params["run"]:
        
        # anchoring parameters
        file.write("Number of neighbors: "+str(anchor_params["n_neighbor"])+"\n")
        file.write("Final smoothing iterations: "+str(anchor_params["n_smooth"])+"\n")
        
        # number of anchoring iterations
        file.write("Number of control points: "+str(niter_a)+"\n")
       
    # registration
    file.write("apply registering: "+str(reg_params["run"]))
    if reg_params["run"]:
        
        # registration parameters
        file.write("T2s contrast: "+str(reg_params["t2s"])+"\n")
        file.write("Shift direction: "+str(reg_params["line_dir"])+"\n")
        file.write("Shift length: "+str(reg_params["line_length"])+"\n")
        file.write("Neighborhood radius: "+str(reg_params["r_size"])+"\n")
        file.write("Learning rate: "+str(reg_params["l_rate"])+"\n")
        file.write("Maximum iterations: "+str(reg_params["max_iter"])+"\n")
        file.write("Cost (threshold): "+str(reg_params["cost_threshold"])+"\n")        
        file.write("Gradient (gaussian): "+str(reg_params["gradient_sigma"])+"\n")
        file.write("Gradient (kernel): "+str(reg_params["gradient_kernel"])+"\n")
        file.write("Overwrite control points: "+str(reg_params["overwrite_control"])+"\n")
        file.write("Cost (step size): "+str(reg_params["cost_step"])+"\n")
        file.write("Cost (sample size): "+str(reg_params["cost_sample"])+"\n")
        file.write("Deformation (gaussian): "+str(reg_params["deformation_sigma"])+"\n")
        
        # number of registration iterations
        if gbb_params is not None:
            for i in range(len(gbb_params["n_iter_step"])):
                file.write("Number of iterations at step "+str(i)+": "+str(gbb_params["n_iter_step"][i])+"\n")
        
            file.write("Registration converged: "+str(gbb_params["convergence"])+"\n")
            file.write("Final number of iterations: "+str(gbb_params["n_iter"])+"\n")
            file.write("Final number of skipped iterations: "+str(gbb_params["n_skip"])+"\n")
    
    # close textfile
    file.close()