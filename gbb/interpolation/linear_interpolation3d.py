# -*- coding: utf-8 -*-

# external inputs
import numpy as np


def linear_interpolation3d(x, y, z, arr_c):
    """
    This function computes the trilinear interpolation for a list of 3D points 
    with coordinates [x,y,z]. Coordinates are in voxel space.
    Inputs:
        *x (list): x-coordinates.
        *y (list): y-coordinates.
        *z (list): z-coordinates.
        *arr_c (arr): 3D array with input values.
    Outputs:
        *c (list): interpolated values for [x,y,z].
        
    created by Daniel Haenelt
    Date created: 29-10-2019
    Last modified: 05-10-2020
    """

    # corner points
    x0 = np.floor(x).astype(int)
    x1 = np.ceil(x).astype(int)
    y0 = np.floor(y).astype(int)
    y1 = np.ceil(y).astype(int)
    z0 = np.floor(z).astype(int)
    z1 = np.ceil(z).astype(int)
    
    # distances to corner points
    xd = (x - x0) / (x1 - x0)
    yd = (y - y0) / (y1 - y0)
    zd = (z - z0) / (z1 - z0)
    
    # corner values
    c000 = arr_c[x0,y0,z0]
    c001 = arr_c[x0,y0,z1]
    c010 = arr_c[x0,y1,z0]
    c011 = arr_c[x0,y1,z1]
    c100 = arr_c[x1,y0,z0]
    c101 = arr_c[x1,y0,z1]
    c110 = arr_c[x1,y1,z0]
    c111 = arr_c[x1,y1,z1]
    
    # interpolation along x-axis
    c00 = c000 * (1 - xd) + c100 * xd
    c01 = c001 * (1 - xd) + c101 * xd
    c10 = c010 * (1 - xd) + c110 * xd
    c11 = c011 * (1 - xd) + c111 * xd
    
    # interpolation along y-axis
    c0 = c00 * (1 - yd) + c10 * yd
    c1 = c01 * (1 - yd) + c11 * yd
    
    # interpolation along z-axis
    c = c0 * (1 - zd) + c1 * zd
    
    return c
