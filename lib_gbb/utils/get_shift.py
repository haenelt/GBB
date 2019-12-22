def get_shift(vtx, fac, n, ind, grad_array, vox2ras_tkr, ras2vox_tkr, vol_max, mm_line=3, 
              n_line=1000, delta_dir=2, t2s=True, show_plot=True):
    """
    This function computes the vertex shift in one direction towards the highest GM/WM gradient. 
    First, the vertex end points of the line are computed and it is checked which side points
    towards WM. A line between both points is defined and values of the second order gradient are
    interpolated onto line points. Starting from the WM end, the first zero crossing is found (it
    exists). If a t2s contrast is considered, the intensity values are goind from dark to bright for
    WM -> GM. Therefore, we expected a transition from positive to negative in this case for the
    around the zero crossing in the second order gradient.
    Inputs:
        *vtx: array of vertex points.
        *fac: array of corresponding faces.
        *n: array of corresponding directions along one axis.
        *ind: current vertex index.
        *grad_array: 3D array of second order gradient values along one axis. 
        *vox2ras: voxel to ras transformation matrix.
        *ras2vox: ras to voxel transformation matrix.
        *vol_max: array of maximum voxel coordinates in x-, y-, and z-direction. 
        *mm_line: length of vertex shift in one direction.
        *n_line: number of line poits.
        *delta_dir: line direction.
        *t2s: wm darker than gm (boolean).
        *show_plot: show line plot in command window.
    Outputs:
        *shift_curr: shift of vertex in ras coordinates.
        
    created by Daniel Haenelt
    Date created: 21-12-2019    
    Last modified: 22-12-2019
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from nibabel.affines import apply_affine
    from lib_gbb.interpolation import linear_interpolation3d

    # initialize shift
    shift_curr = None
    
    # get current vertex and normal
    vtx_curr = vtx[ind,:]
    n_curr = n[ind]
    
    # get line a -> b (WM -> GM)
    pt_start = vtx_curr.copy()
    pt_end = vtx_curr.copy()
    if n_curr == 1:
        pt_start[delta_dir] += mm_line
        pt_end[delta_dir] -= mm_line
    elif n_curr == -1:
        pt_start[delta_dir] -= mm_line
        pt_end[delta_dir] += mm_line
    else:
        return shift_curr
        
    line_curr = np.zeros((n_line,3), dtype=np.float)
    line_curr[:,delta_dir] = np.linspace(0,1,n_line)    
    line_curr = (pt_end-pt_start) * line_curr + pt_start
    line_curr = apply_affine(ras2vox_tkr, line_curr)
    
    # remove coordinates exceeding the array limits
    outlier = np.zeros(len(line_curr), dtype=np.int8)
    outlier[line_curr[:,0] < 0] = 1
    outlier[line_curr[:,1] < 0] = 1
    outlier[line_curr[:,2] < 0] = 1
    outlier[line_curr[:,0] > vol_max[0]] = 1
    outlier[line_curr[:,1] > vol_max[1]] = 1
    outlier[line_curr[:,2] > vol_max[2]] = 1
    
    line_curr = line_curr[outlier == 0,:]
    
    # get gradient values along line
    grad_curr = linear_interpolation3d(line_curr[:,0],line_curr[:,1],line_curr[:,2],grad_array)
    
    # show line plot
    if show_plot:
        plt.plot(grad_curr)
        plt.xlabel("WM -> GM")
        plt.ylabel("Second order gradient")
    
    # get point of zero crossing (move from wm to gm)
    i = 0
    while i < n_line - 1:
        
        if grad_curr[i] > 0 and grad_curr[i+1] < 0 and t2s:
            zero_curr = apply_affine(vox2ras_tkr, line_curr[i,:])
            shift_curr = vtx_curr - zero_curr
            break
        elif grad_curr[i] < 0 and grad_curr[i+1] > 0 and not t2s:
            zero_curr = apply_affine(vox2ras_tkr, line_curr[i,:])
            shift_curr = vtx_curr - zero_curr
            break
        
        i += 1
    
    return shift_curr