def get_normal(vtx, fac):
    """
    This function computes the surfaces normals per vertex from an input surface mesh. The code is 
    taken from https://sites.google.com/site/dlampetest/python/calculating-normals-of-a-triangle-
    mesh-using-numpy
    Inputs:
        *vtx: vertex array.
        *fac: face array.
    Outputs:
        *norm: normal per vertex.
        
    created by Daniel Haenelt
    Date created: 13-12-2019
    Last modified: 18-12-2019
    """
    import numpy as np

    # normalize vectors
    def normalize_v3(arr):
        """ Normalize a numpy array of 3 component vectors shape=(n,3) """
        lens = np.sqrt( arr[:,0]**2 + arr[:,1]**2 + arr[:,2]**2 )
        arr[:,0] /= lens
        arr[:,1] /= lens
        arr[:,2] /= lens
        return arr

    # initialize array with same type and shape as vtx
    norm = np.zeros( vtx.shape, dtype=vtx.dtype )

    # indexed view into the vertex array
    tris = vtx[fac]    
    
    # calculate the normal for all triangles by taking the cross product of the vectors v1-v0 and
    # v2-v0 in each triangle
    n = np.cross( tris[::,1 ] - tris[::,0]  , tris[::,2 ] - tris[::,0] )
    
    # n is array of normal per triangle. Now we normalize all normals
    normalize_v3(n)
    
    # now we have a normalized array of normals, one per triangle. We add each vertex in that
    # triangle to get one normal per vertex.
    norm[ fac[:,0] ] += n
    norm[ fac[:,1] ] += n
    norm[ fac[:,2] ] += n
    normalize_v3(norm)
    
    return norm
