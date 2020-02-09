def check_exit(y_array, cost_threshold=None):
    """
    This function computes a linear fit from y-coordinates of given sample points. x-coordinates are
    assumed to be equidistantly spaced integer values. It is checked if the slope of the fit is
    below a set threshold value.
    Inputs.
        *y_array: input vector of y-coordinates.
        *cost_threshold: threshold value.
    Outputs:
        *m: slope of the lienear fit.
        *exit_criterion: boolean which indicate if m is below threshold value.
    
    created by Daniel Haenelt
    Date created: 09-02-2020
    Last modified: 09-02-2020
    """
    import numpy as np

    # define x-array
    x_array = np.arange(len(y_array))
    
    # make linear fit
    m, _ = np.polyfit(x_array, y_array, deg=1, rcond=None, full=False, w=None, cov=False)
    
    # make abs
    m = np.abs(m)
    
    # check threshold
    if cost_threshold and m < cost_threshold:
        exit_criterion = True
    else:
        exit_criterion = False
        
    return m, exit_criterion