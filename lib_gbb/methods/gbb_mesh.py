def gbb_mesh(vtx, fac, n_dir, normal, ind_control, arr_ref, arr_gradient, arr_vein, arr_ignore, t2s, 
             vox2ras_tkr, ras2vox_tkr, line_dir, line_length, r_size, l_rate, max_iterations, 
             cost_threshold, cost_step=1000, cost_sample=10, path_output = "", show_cost=True, 
             show_slope=False, write_temp=10000):
    """
    This function applies the core function of the gbb method. 
    Inputs.
        *vtx: array of vertex points.
        *fac: array of corresponding faces.
        *n_dir: array normal directions.
        *normal: array of normal vectors.
        *ind_control: array of vertex indices matching control points.
        *arr_ref: array of reference volume.
        *arr_gradient: array of gradient image.
        *arr_vein: array of vein mask.
        *arr_ignore: array of ignore mask.
        *t2s: image contrast (boolean).
        *vox2ras_tkr: transformation matrix from voxel space to ras.
        *ras2vox_tkr: transformation matrix from ras to voxel space.
        *line_dir: line axis in ras convention (0,1,2,3).
        *line_length: line length in one direction in mm.
        *r_size: neighborhood radius in mm (list).
        *l_rate: learning rate (list).
        *max_iterations: maximum iterations (list).
        *cost_threshold: cost function threshold (list).
        *cost_step: step size between cost array points.
        *cost_sample: sample size for linear fit.
        *path_output: path where intermediate folder is created.
        *show_cost: show temporary cost function (boolean).
        *show_slope: show temporary slope function (boolean).
        *write_temp: step size to write intermediate surfaces (if set > 0).
    Outputs:
        *vtx: shifted array of vertex points.
        *gbb: collection of performance parameters (dict).
    
    created by Daniel Haenelt
    Date created: 06-02-2020          
    Last modified: 13-05-2020
    """
    import os
    import numpy as np
    from nibabel.freesurfer.io import write_geometry
    from lib_gbb.neighbor.nn_3d import nn_3d
    from lib_gbb.utils.get_shift import get_shift
    from lib_gbb.utils.cost_BBR import cost_BBR
    from lib_gbb.utils.update_mesh import update_mesh
    from lib_gbb.utils.check_exit import check_exit
    from lib_gbb.utils.get_ignore import get_ignore
    from lib_gbb.plot.cost_plot import cost_plot
    from lib_gbb.plot.slope_plot import slope_plot

    # make intermediate folder
    if write_temp > 0 and path_output:
        tmp = np.random.randint(0, 10, 5)
        tmp_string = ''.join(str(i) for i in tmp)
        path_temp = os.path.join(path_output,"tmp_"+tmp_string)
        if not os.path.exists(path_temp):
            os.makedirs(path_temp)
    
    # get vertices to ignore
    if arr_ignore is not None:
        _, ind_ignore = get_ignore(vtx, arr_ignore, ras2vox_tkr, write_output=False, path_output=False)
    else:
        ind_ignore = []
    
    print("start mesh initialization (gbb)")    
    
    # run surface reginement
    i = 0
    j = 0
    p = 0
    q = 0
    counter = 0
    step = 0
    cost_array = []
    m_array = []
    n_array = [] 
    n_coords = len(vtx)
    n_iter_step = np.zeros(len(r_size))
    vol_max = np.shape(arr_ref)
    while True:
        
        # get current vertex point
        n_vertex = np.random.randint(n_coords)
        if n_vertex in ind_ignore:
            counter += 1
            continue
        
        # get shift
        vtx_shift = get_shift(vtx, fac, n_dir, n_vertex, 
                              arr_gradient, arr_vein, vox2ras_tkr, ras2vox_tkr, 
                              vol_max, line_length, 
                              line_dir, t2s, False)
           
        # update mesh
        if len(vtx_shift):
            nn_ind, _ = nn_3d(vtx[n_vertex], vtx, r_size[step])
            if np.any(np.in1d(nn_ind,ind_control)):
                counter += 1
                continue     
            vtx = update_mesh(vtx,vtx_shift,n_vertex,nn_ind,l_rate[step])
        else:
            counter += 1
            continue        
        
        # get cost function
        if p >= cost_step:
            J = cost_BBR(vtx, normal, arr_ref, ras2vox_tkr, vol_max, 
                         t2s)
            cost_array = np.append(cost_array, J)
            q += 1
            p = 0
            if len(cost_array) >= cost_sample:
                
                # check exit criterion
                m_fit, n_fit, exit_crit = check_exit(np.arange(q-cost_sample,q), 
                                                     cost_array[-cost_sample:], 
                                                     cost_threshold=cost_threshold[step])
                
                # save computed slope and y-axis intercept of liner fit
                m_array = np.append(m_array, m_fit)
                n_array = np.append(n_array, n_fit)
                
                # shot plots
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
            n_iter_step[step] = j
            step += 1
            j = 0
            print("start registration step "+str(step)+" at iteration "+str(i))
        elif exit_crit and step == len(max_iterations) - 1:
            n_iter_step[step] = j
            gbb_converged = True
            print("Registration converged!")   
            break
        elif step == len(max_iterations) - 1 and j == max_iterations[-1]:
            n_iter_step[step] = j
            gbb_converged = False
            print("Registration did not converge!")
            break
    
        # write intermediate surfaces
        if write_temp and path_output and not np.mod(i,write_temp):
            write_geometry(os.path.join(path_temp,"temp_"+str(i)), vtx, fac)
        
        i += 1
        j += 1
    
    # print some information
    print("Final number of iterations: "+str(i))
    print("Final number of skipped iterations: "+str(counter))
    
    # collect some descriptive variables
    gbb = dict()
    gbb["convergence"] = gbb_converged
    gbb["cost_array"] = cost_array
    gbb["m_array"] = m_array
    gbb["n_array"] = n_array
    gbb["n_iter"] = i
    gbb["n_skip"] = counter
    gbb["n_iter_step"] = n_iter_step
    
    return vtx, gbb