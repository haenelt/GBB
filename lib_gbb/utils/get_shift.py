def get_shift(vtx, fac, n, ind, arr_grad, arr_vein, vox2ras_tkr, ras2vox_tkr, vol_max, 
              line_length=3, line_dir=2, t2s=True, show_plot=True):
    """
    This function computes the vertex shift in one direction towards the highest GM/WM gradient. 
    First, vertex end points of the line are computed and it is checked which side points towards 
    WM. A line between both points is calculated and values of the second order gradient are 
    sampled onto the line. Starting from the current vertex point, the nearest zero crossing
    is found (if it exists). If a t2s contrast is considered, the intensity values are going from 
    dark to bright for WM -> GM. Therefore, we expect a transition from positive to negative in 
    this case for the zero crossing in the second order gradient. The found shift is only considered
    if no vein is found in shift direction along the line. Vertex coordinates are in ras space.
    Inputs:
        *vtx: array of vertex points.
        *fac: array of corresponding faces.
        *n: array of corresponding directions along one axis.
        *ind: current vertex index.
        *arr_grad: 3D array of second order gradient values along one axis.
        *arr_vein: 3D array with masked veins.
        *vox2ras: voxel to ras transformation matrix.
        *ras2vox: ras to voxel transformation matrix.
        *vol_max: array of maximum voxel coordinates in x-, y-, and z-direction.
        *line_length: length of vertex shift in one direction in mm.
        *line_dir: line direction in ras conventions.
        *t2s: wm darker than gm (boolean).
        *show_plot: show line plot in command window.
    Outputs:
        *shift_curr: shift of vertex in ras coordinates.
        
    created by Daniel Haenelt
    Date created: 21-12-2019
    Last modified: 13-05-2020
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from nibabel.affines import apply_affine
    from lib_gbb.interpolation import linear_interpolation3d

    # initialize
    shift_curr = []
    zero_found = []
    
    # fix parameters
    n_line = 1000 # number of line points
    
    # get current vertex and normal
    vtx_curr = vtx[ind,:]
    n_curr = n[ind]
    
    # get line a -> b (WM -> GM)
    pt_start = vtx_curr.copy()
    pt_end = vtx_curr.copy()
    if n_curr == 1:
        pt_start[line_dir] += line_length
        pt_end[line_dir] -= line_length
    elif n_curr == -1:
        pt_start[line_dir] -= line_length
        pt_end[line_dir] += line_length
    else:
        return []
       
    line_curr = np.zeros((n_line,3), dtype=np.float)
    line_curr[:,line_dir] = np.linspace(0,1,n_line)    
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
    n_line = len(line_curr) # update line length
    
    # get gradient and vein data along line
    grad_curr = linear_interpolation3d(line_curr[:,0],line_curr[:,1],line_curr[:,2],arr_grad)
    vein_curr = linear_interpolation3d(line_curr[:,0],line_curr[:,1],line_curr[:,2],arr_vein)
    
    # show line plot
    if show_plot:
        plt.plot(grad_curr)
        plt.xlabel("WM -> GM")
        plt.ylabel("Second order gradient")
    
    # mid-point of line
    i = np.floor(n_line / 2).astype(int)
    
    # get point of zero crossing (closest to mid-point) if grad_curr exists, i.e., if it does not 
    # contain only outliers
    j = 0
    switch = 0
    if len(grad_curr):
    
        # look for veins in line
        vein_up = np.sum(vein_curr[i:])
        vein_down = np.sum(vein_curr[:i])

        # start search
        while i > 0 and i < n_line - 2:
            
            # if vein up -> i nur in negativer richtung
            # if vein down -> i nur in positiver richtung
            # if not -> look in both directions (switch)
            # if both -> return empty array
            if vein_up and not vein_down:
                i -= 1
            elif not vein_up and vein_down:
                i += 1
            elif not vein_up and not vein_down:
                j += 1
                if switch == 0:
                    switch = 1
                    i += j
                else:
                    switch = 0
                    i -= j
            else:
                break
            
            if grad_curr[i] < 0 and grad_curr[i+1] > 0 and t2s:
                zero_found = 1
                break
            elif grad_curr[i] > 0 and grad_curr[i+1] < 0 and not t2s:
                zero_found = 1
                break
        
    # only consider shift if not positioned within a vein
    if zero_found:
        zero_curr = apply_affine(vox2ras_tkr, line_curr[i,:])
        shift_curr = vtx_curr - zero_curr
    
    return shift_curr