# -*- coding: utf-8 -*-

# external inputs
import numpy as np


def check_exit(arr_x, arr_y, cost_threshold=None):
    """ Check exit

    This function computes a linear fit from x- and y-coordinates of given 
    sample points. It is checked if the slope of the fit is below a set 
    threshold value.    

    Parameters
    ----------
    arr_x : ndarray
        Input vector of x-coordinates.
    arr_y : ndarray
        Input vector of y-coordinates.
    cost_threshold : float, optional
        Threshold value. The default is None.

    Returns
    -------
    m : float
        Slope of the linear fit.
    n : float
        y-intercept of the linear fit.
    exit_criterion : bool
        True if m is below threshold value.

    Notes
    -------
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
