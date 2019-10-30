def linear_interpolation2d(x, y, f_array):
    """
    This function computes a bilinear interpolation at the point [x,y] from the nearest neighbour 
    points [x_s,y_s].
    Inputs:
        *x: coordinate for which the interpolation is calculated.
        *y: coordinate for which the interpolation is calculated.
        *f_array: 2D array with input values.
    Outputs:
        *f: interpolated value for [x, y].
        
    created by Daniel Haenelt
    Date created: 29-10-2019           
    Last modified: 30-10-2019
    """
    import numpy as np
    
    # get nearest neighbour grid points
    x_s = np.array([np.floor(x).astype(int), np.ceil(x).astype(int)])
    y_s = np.array([np.floor(y).astype(int), np.ceil(y).astype(int)]) 
    
    # get grid point array values
    f_s = [[f_array[y_s[0],x_s[0]],f_array[y_s[1],x_s[0]]],
           [f_array[y_s[0],x_s[1]],f_array[y_s[1],x_s[1]]]]
       
    # get distance between grid points
    Delta_x = x_s[1] - x_s[0]
    Delta_y = y_s[1] - y_s[0]
    
    # compute interpolation
    A = 1 / (Delta_x*Delta_y)
    B = np.array([x_s[1] - x, x - x_s[0]])
    C = [y_s[1] - y, y - y_s[0]]
    
    f = A * np.dot(B,np.dot(f_s,C))
    
    return f