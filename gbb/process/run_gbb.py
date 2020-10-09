# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
from nibabel.freesurfer.io import write_geometry

# local inputs
from gbb.neighbor.nn_3d import nn_3d
from gbb.utils.get_shift import get_shift
from gbb.utils.cost_BBR import cost_BBR
from gbb.utils.update_mesh import update_mesh
from gbb.utils.check_exit import check_exit
from gbb.utils.get_ignore import get_ignore
from gbb.plot.cost_plot import cost_plot
from gbb.plot.slope_plot import slope_plot


def run_gbb(vtx, fac, vtx_n, ind_control, arr_ref, arr_gradient, arr_vein, 
            arr_ignore, t2s, vox2ras_tkr, ras2vox_tkr, line_dir, line_length, 
            r_size, l_rate, max_iterations, cost_threshold, cost_step=1000, 
            cost_sample=10, path_output = "", show_cost=True, show_slope=False, 
            write_temp=10000, Q0=0, M=0.5, h=1, s=1):    
    """
    This function applies the core function of the gbb method. 
    Inputs.
        *vtx (arr): array of vertex points.
        *fac (arr): array of corresponding faces.
        *vtx_n (arr): array of vertex normals.
        *ind_control (arr): array of vertex indices matching control points.
        *arr_ref (arr): array of reference volume.
        *arr_gradient (arr): array of gradient image.
        *arr_vein (arr): array of vein mask.
        *arr_ignore (arr): array of ignore mask.
        *t2s (bool): image contrast.
        *vox2ras_tkr (arr): transformation matrix from voxel space to ras.
        *ras2vox_tkr (arr): transformation matrix from ras to voxel space.
        *line_dir (int): line axis in ras convention (0,1,2,3).
        *line_length (float): line length in one direction in mm.
        *r_size (list): neighborhood radius in mm.
        *l_rate (list): learning rate.
        *max_iterations (list): maximum iterations.
        *cost_threshold (list): cost function threshold.
        *cost_step (int): step size between cost array points.
        *cost_sample (int): sample size for linear fit.
        *path_output (str): path where intermediate folder is created.
        *show_cost (bool): show temporary cost function.
        *show_slope (bool): show temporary slope function.
        *write_temp (int): step size to write intermediate surfaces (if set > 0).
        *Q0 (float): const function offset parameter.
        *M (float): cost function slope parameter.
        *h (float): cost function weight parameter.
        *s (float): cost function distance parameter.
    Outputs:
        *vtx (arr): shifted array of vertex points.
        *gbb (dict): collection of performance parameters.
    
    created by Daniel Haenelt
    Date created: 06-02-2020          
    Last modified: 08-10-2020
    """

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
    n_iter_step = []
    n_coords = len(vtx)
    while True:
        
        # get current vertex point
        n_vertex = np.random.randint(n_coords)
        if n_vertex in ind_ignore:
            counter += 1
            continue
        
        # get shift        
        vtx_shift = get_shift(vtx, fac, vtx_n, n_vertex, arr_gradient, 
                              vox2ras_tkr, ras2vox_tkr, arr_vein, line_dir, 
                              line_length, t2s, False)
        
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
            J = cost_BBR(vtx, vtx_n, arr_ref, ras2vox_tkr, Q0=Q0, M=M, h=h, s=s, t2s=t2s)
            cost_array.append(J)
            q += 1
            p = 0
            if len(cost_array) >= cost_sample:
                
                # check exit criterion
                m_fit, n_fit, exit_crit = check_exit(np.arange(q-cost_sample,q), 
                                                     cost_array[-cost_sample:], 
                                                     cost_threshold=cost_threshold[step])
                
                # save computed slope and y-axis intercept of liner fit
                m_array.append(m_fit)
                n_array.append(n_fit)
                
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
            n_iter_step.append(j)
            step += 1
            j = 0
            print("start registration step "+str(step)+" at iteration "+str(i))
        elif exit_crit and step == len(max_iterations) - 1:
            n_iter_step.append(j)
            gbb_converged = True
            print("registration converged!")   
            break
        elif step == len(max_iterations) - 1 and j == max_iterations[-1]:
            n_iter_step.append(j)
            gbb_converged = False
            print("registration did not converge!")
            break
    
        # write intermediate surfaces
        if write_temp and path_output and not np.mod(i,write_temp):
            write_geometry(os.path.join(path_temp,"temp_"+str(i)), vtx, fac)
        
        i += 1
        j += 1
    
    # print some information
    print("final number of iterations: "+str(i))
    print("final number of skipped iterations: "+str(counter))
    
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
