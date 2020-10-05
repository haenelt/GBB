# -*- coding: utf-8 -*-

# external inputs
import numpy as np


def check_exit(arr_x, arr_y, cost_threshold=None):
    """
    This function computes a linear fit from x- and y-coordinates of given 
    sample points. It is checked if the slope of the fit is below a set 
    threshold value.
    Inputs.
        *arr_x (arr): input vector of x-coordinates.
        *arr_y (arr): input vector of y-coordinates.
        *cost_threshold (float): threshold value.
    Outputs:
        *m (float): slope of the linear fit.
        *n (float): y-intercept of the linear fit.
        *exit_criterion (bool): if m is below threshold value.
    
    created by Daniel Haenelt
    Date created: 09-02-2020
    Last modified: 05-10-2020
    """
    
    # make linear fit
    m, n = np.polyfit(arr_x, arr_y, deg=1, rcond=None, full=False, w=None, cov=False)
        
    # check threshold
    if cost_threshold and np.abs(m) < cost_threshold:
        exit_criterion = True
    else:
        exit_criterion = False
        
    return m, n, exit_criterion
