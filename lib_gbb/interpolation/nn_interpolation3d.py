def nn_interpolation3d(x, y, z, arr_c):
    """
    This function computes the nearest neighbor interpolation of an 3D array at points [x,y,z].
    Coordinates are in voxel space.
    Inputs:
        *x: list of x coordinates.
        *y: list of y coordinates.
        *z: list of z coordinates.
        *arr_c: 3D array with input values.
    Outputs:
        *c: list of interpolated values for [x,y,z].
        
    created by Daniel Haenelt
    Date created: 30-10-2019
    Last modified: 20-12-2019
    """
    import numpy as np

    # get nearest neighbour grid points
    x0 = np.round(x).astype(int)
    y0 = np.round(y).astype(int) 
    z0 = np.round(z).astype(int)
    
    c = arr_c[x0,y0,z0]
    
    return c