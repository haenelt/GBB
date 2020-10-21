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
from gbb.plot.plot_cost import plot_cost
from gbb.plot.plot_slope import plot_slope


def run_gbb(vtx, fac, vtx_n, ind_control, arr_ref, arr_gradient, arr_vein, 
            arr_ignore, t2s, vox2ras_tkr, ras2vox_tkr, line_dir, line_length, 
            r_size, l_rate, max_iterations, cost_threshold, cost_step=1000, 
            cost_sample=10, path_output = "", show_cost=True, show_slope=False, 
            write_temp=10000, Q0=0, M=0.5, h=1, s=1):
    """ Run GBB

    This function applies the core function of the gbb method.     

    Parameters
    ----------
    vtx : ndarray
        Array of vertex points.
    fac : ndarray
        Array of corresponding faces.
    vtx_n : ndarray
        Array of vertex normals.
    ind_control : ndarray
        Array of vertex indices matching control points.
    arr_ref : ndarray
        Array of reference volume.
    arr_gradient : ndarray
        Array of gradient image.
    arr_vein : ndarray
        Array of vein mask.
    arr_ignore : ndarray
        Array of ignore mask.
    t2s : bool
        Image contrast.
    vox2ras_tkr : ndarray
        Transformation matrix from voxel space to ras.
    ras2vox_tkr : ndarray
        Transformation matrix from ras to voxel space.
    line_dir : int
        Line axis in ras convention (0,1,2,3).
    line_length : float
        Line length in one direction in mm.
    r_size : list
        Neighborhood radius in mm.
    l_rate : list
        Learning rate.
    max_iterations : list
        Maximum iterations.
    cost_threshold : list
        Cost function threshold.
    cost_step : int, optional
        Step size between cost array points. The default is 1000.
    cost_sample : int, optional
        Sample size for linear fit. The default is 10.
    path_output : str, optional
        Path where intermediate folder is created. The default is "".
    show_cost : bool, optional
        Show temporary cost function. The default is True.
    show_slope : bool, optional
        Show temporary slope function. The default is False.
    write_temp : int, optional
        Step size to write intermediate surfaces (if set > 0). The default is 
        10000.
    Q0 : float, optional
        Const function offset parameter. The default is 0.
    M : float, optional
        Cost function slope parameter. The default is 0.5.
    h : float, optional
        Cost function weight parameter. The default is 1.
    s : float, optional
        Cost function distance parameter. The default is 1.

    Returns
    -------
    vtx : ndarray
        Shifted array of vertex points.
    gbb : dict
        Collection of performance parameters.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 06-02-2020          
    Last modified: 21-10-2020

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
                    plot_cost(q, cost_array, m_fit, n_fit, set_title, save_plot=False, 
                              path_output=False, name_output=False)
                    
                if show_slope:
                    set_title = "Exit criterion: "+str(exit_crit)+", Step: "+str(step)
                    plot_slope(q, m_array, set_title, save_plot=False, path_output=False, 
                               name_output=False)
        else:
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
        p += 1
    
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
