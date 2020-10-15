# -*- coding: utf-8 -*-

# external inputs
import numpy as np


def nn_interpolation3d(x, y, z, arr_c):
    """ NN interpolation 3D

    This function computes the nearest neighbor interpolation for a list of 3D 
    points with coordinates [x,y,z]. Coordinates are in voxel space.    

    Parameters
    ----------
    x : list
        x-coordinates.
    y : list
        y-coordinates.
    z : list
        z-coordinates.
    arr_c : list
        3D array with input values.

    Returns
    -------
    c : list
        Interpolated values for [x,y,z].

    Notes
    -------
    created by Daniel Haenelt
    Date created: 30-10-2019
    Last modified: 05-10-2020

    """
    
    # get nearest neighbour grid points
    x0 = np.round(x).astype(int)
    y0 = np.round(y).astype(int) 
    z0 = np.round(z).astype(int)
    
    c = arr_c[x0,y0,z0]
    
    return c
