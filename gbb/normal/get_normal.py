# -*- coding: utf-8 -*-

# external inputs
import numpy as np


def get_normal(vtx, fac):
    """
    This function computes the surfaces normals per vertex from an input surface 
    mesh. The code is taken from https://sites.google.com/site/dlampetest/
    python/calculating-normals-of-a-triangle-mesh-using-numpy and adapted to my 
    own purposes.
    Inputs:
        *vtx (arr): vertex array.
        *fac (arr): face array.
    Outputs:
        *norm (arr): normal per vertex.
        
    created by Daniel Haenelt
    Date created: 13-12-2019
    Last modified: 05-10-2020
    """

    def normalize_v3(arr):
        """ Normalize a numpy array of 3 component vectors shape=(n,3) """
        lens = np.sqrt( arr[:,0]**2 + arr[:,1]**2 + arr[:,2]**2 )
        res = np.zeros_like(arr)
        res[:,0] = arr[:,0] / lens
        res[:,1] = arr[:,1] / lens
        res[:,2] = arr[:,2] / lens
        return res

    def f2v(length, fac, nt):
        """ get average vertex-wise normal by adding up all face-wise normals
        around vertex """
        nv = np.zeros((length,3))
        for i in range(len(fac)):       
            nv[fac[i,0],:] += nt[i,:]
            nv[fac[i,1],:] += nt[i,:]
            nv[fac[i,2],:] += nt[i,:]
        
        return nv

    # initialize array with same type and shape as vtx
    norm = np.zeros( vtx.shape, dtype=vtx.dtype )

    # indexed view into the vertex array
    tris = vtx[fac]   
    
    # calculate the normal for all triangles by taking the cross product of 
    # the vectors v1-v0 and v2-v0 in each triangle
    n = np.cross( tris[::,1 ] - tris[::,0]  , tris[::,2 ] - tris[::,0] )
    
    # calculate vertex-wise normals and normalize
    norm = f2v(len(vtx), fac, n)
    norm = normalize_v3(norm)
    
    return norm
