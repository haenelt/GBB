def linear_interpolation2d(x, y, x_s, y_s, f_s):
    """
    This function computes a bilinear interpolation at the point [x,y] from four points [x_1,y_1], 
    [x_1, y_2], [x_2, y_1] and [x_2, y_2] on a two-dimensional plane.
    Inputs:
        *x: coordinate for which the interpolation is calculated.
        *y: coordinate for which the interpolation is calculated.
        *x_s: list [x_1, x_2] of grid coordinates.
        *y_s: list [y_1, y_2] of grid coordinates.
        *f_s: list [[f_11, f_12],[f_21, f_22]] of grid values.
    Outputs:
        *f: interpolated value for [x, y].
        
    created by Daniel Haenelt
    Date created: 29-10-2019           
    Last modified: 29-10-2019
    """
    import numpy as np
    
    # convert to numpy array
    x_s = np.array(x_s)
    y_s = np.array(y_s)
    f_s = np.array(f_s)
    
    # get distance between grid points
    Delta_x = x_s[1] - x_s[0]
    Delta_y = y_s[1] - y_s[0]
    
    # compute interpolation
    A = 1 / (Delta_x*Delta_y)
    B = np.array([x_s[1] - x, x - x_s[0]])
    C = [y_s[1] - y, y - y_s[0]]
    
    f = A * np.dot(B,np.dot(f_s,C))
    
    return f