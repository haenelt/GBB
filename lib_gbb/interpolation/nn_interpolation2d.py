def nn_interpolation2d(x, y, f_array):
    """
    This function computes a nearest neighbour interpolation at the point [x,y].
    Inputs:
        *x: coordinate for which the interpolation is calculated.
        *y: coordinate for which the interpolation is calculated.
        *f_array: 2D array with input values.
    Outputs:
        *f: interpolated value for [x, y].
        
    created by Daniel Haenelt
    Date created: 30-10-2019           
    Last modified: 30-10-2019
    """
    import numpy as np
    
    # get nearest neighbour grid point
    x_s = np.round(x).astype(int)
    y_s = np.round(y).astype(int) 
    
    f = f_array[y_s,x_s]
    
    return f