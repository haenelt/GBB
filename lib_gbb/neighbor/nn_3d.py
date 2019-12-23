def nn_3d(vtx0, vtx, r_size):
    """
    This function all nearest neighbors found within a sphere with a radius defined in ras
    coordinates. Note that the defined neighborhood does not have to be fully connected in this
    case. Vertex coordinates are in ras space.
    Inputs:
        *vtx0: vertex point.
        *vtx: array of vertices.
        *r_size: radius of sphere in ras coordinates.
    Outputs:
        *nn: array of neighbor indices.
        
    created by Daniel Haenelt
    Date created: 22-12-2019     
    Last modified: 22-12-2019
    """
    import numpy as np

    rx = ( vtx[:,0] - vtx0[0] ) ** 2
    ry = ( vtx[:,1] - vtx0[1] ) ** 2
    rz = ( vtx[:,2] - vtx0[2] ) ** 2
    
    r = np.sqrt( rx + ry + rz )
    
    nn = np.where(r < r_size)[0]
    
    return nn