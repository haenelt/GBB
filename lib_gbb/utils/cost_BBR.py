def cost_BBR(vtx, fac, vol_array, n, ras2vox, delta_size=2, delta_dir=2, Q0=0, M=0.5, h=1, 
             t2s=True):
    """
    This function the cost function defined in the original BBR paper (Fischl et al., 2009). 
    Vertices should be in freesurfer ras_tkr space. First, vertex coordinates on both sides of the
    surface mesh along one axis are computed. GM and WM coordinates are the sorted based on surface
    normals. Based on a tranformation matrix, vertex coordinates are transformed into the voxel
    space of the input array. Vertices are excluded with normals perpendicular to the offset 
    direction and with coordinates exceeding the array limits. GM and WM values are are sampled 
    using linear interpolation. A percent contrast measure and the cost function value are then 
    calculated.
    Inputs:
        *vtx: array of vertices.
        *fac: array of faces.
        *vol_array: 3D array of image volume.
        *n: normal direction.
        *ras2vox: transformation matrix to voxel space.
        *delta_size: offset from surface in mm.
        *delta_dir: axis for offset in the ras coordinate system (0,1,2). 
        *Q0: offset parameter in percent contrast measure.
        *M: slope parameter in percent contrast measure.
        *h: weight for each vertex in percent contrast measure.
        *t2s: t2star weighing, i.e., gm > wm intensity (boolean).
    Outputs:
        *J: cost function value.
        
    created by Daniel Haenelt
    Date created: 21-12-2019
    Last modified: 21-12-2019
    """
    import numpy as np
    from nibabel.affines import apply_affine
    from lib_gbb.interpolation import linear_interpolation3d

    # get maximum vertex coordinates in voxel space
    vol_xmax = np.shape(vol_array)[0] - 1
    vol_ymax = np.shape(vol_array)[1] - 1
    vol_zmax = np.shape(vol_array)[2] - 1
    
    # get shifted vertices
    vtx_up = vtx.copy()
    vtx_down = vtx.copy()
    vtx_up[:,delta_dir] += delta_size
    vtx_down[:,delta_dir] -= delta_size
    
    # sort offset in two groups according to normal direction
    gm_pts = np.asarray([vtx_up[i,:] if n[i] < 0 else vtx_down[i,:] for i in range(len(vtx))])
    wm_pts = np.asarray([vtx_up[i,:] if n[i] > 0 else vtx_down[i,:] for i in range(len(vtx))])
    
    # remove vertices with normals below threshold (normals perpendicular to offset direction)
    gm_pts = gm_pts[n != 0, :]
    wm_pts = wm_pts[n != 0, :]
    
    # ras2vox transformation
    gm_pts = apply_affine(ras2vox, gm_pts)
    wm_pts = apply_affine(ras2vox, wm_pts)
    
    # get location of outlier coordinates
    outlier = np.zeros(len(gm_pts), dtype=np.int8())
    outlier[gm_pts[:,0] > vol_xmax] = 1
    outlier[gm_pts[:,1] > vol_ymax] = 1
    outlier[gm_pts[:,2] > vol_zmax] = 1
    outlier[wm_pts[:,0] > vol_xmax] = 1
    outlier[wm_pts[:,1] > vol_ymax] = 1
    outlier[wm_pts[:,2] > vol_zmax] = 1
    outlier[gm_pts[:,0] < 0] = 1
    outlier[gm_pts[:,1] < 0] = 1
    outlier[gm_pts[:,2] < 0] = 1
    outlier[wm_pts[:,0] < 0] = 1
    outlier[wm_pts[:,1] < 0] = 1
    outlier[wm_pts[:,2] < 0] = 1
    
    # get values in GM and WM
    gm_val = linear_interpolation3d(gm_pts[:,0], gm_pts[:,1], gm_pts[:,2], vol_array)
    wm_val = linear_interpolation3d(wm_pts[:,0], wm_pts[:,1], wm_pts[:,2], vol_array)
    
    # get rid of outliers
    gm_val = gm_val[outlier == 0]
    wm_val = wm_val[outlier == 0]
    
    # percent contrast measure
    if t2s == True:
        Q = 100 * ( wm_val - gm_val ) / ( 0.5 * ( gm_val + wm_val ) )
    else:
        Q = 100 * ( gm_val - wm_val ) / ( 0.5 * ( gm_val + wm_val ) )
    
    # cost value
    J = 1 / len(Q) * np.sum( h * ( 1 + np.tanh(M * ( Q - Q0 )) ) )
    
    return J