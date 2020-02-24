def check_exit(x_array, y_array, cost_threshold=None):
    """
    This function computes a linear fit from x- and y-coordinates of given sample points. It is 
    checked if the slope of the fit is below a set threshold value.
    Inputs.
        *x_array: input vector of x-coordinates.
        *y_array: input vector of y-coordinates.
        *cost_threshold: threshold value.
    Outputs:
        *m: slope of the linear fit.
        *n: y-intercept of the linear fit.
        *exit_criterion: boolean which indicate if m is below threshold value.
    
    created by Daniel Haenelt
    Date created: 09-02-2020
    Last modified: 23-02-2020
    """
    import numpy as np
    
    # make linear fit
    m, n = np.polyfit(x_array, y_array, deg=1, rcond=None, full=False, w=None, cov=False)
        
    # check threshold
    if cost_threshold and np.abs(m) < cost_threshold:
        exit_criterion = True
    else:
        exit_criterion = False
        
    return m, n, exit_criterion