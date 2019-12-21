def get_shift(vtx, fac, n, ind, grad_array, vox2ras, ras2vox, vol_xmax, vol_ymax, vol_zmax, 
              mm_line=3, n_line=1000, delta_dir=2, t2s=True, show_plot=True):
    """
    This function blas...
    for t2s: should be positive gradient (dark to light)
    Inputs:
        *vtx:
        *fac:
        *n:
        *ind:
        *grad_array:
        *vox2ras:
        *ras2vox:
        *vol_xmax:
        *vol_ymax:
        *vol_zmax:
        *mm_line:
        *n_line:
        *delta_dir:
        *t2s:
        *show_plot:
    Outputs:
        *shift_curr: bla
        
    created by Daniel Haenelt
    Date created: 21-12-2019    
    Last modified: 21-12-2019
    """
    import numpy as np
    from lib_gbb.interpolation import linear_interpolation3d
    from lib.surface.vox2ras import vox2ras
    from nibabel.affines import apply_affine
    import matplotlib.pyplot as plt

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
    line_curr = apply_affine(ras2vox, line_curr)
    
    # remove coordinates exceeding the array limits
    line_curr[line_curr[:,0] < 0,0] = np.nan
    line_curr[line_curr[:,1] < 0,0] = np.nan
    line_curr[line_curr[:,2] < 0,0] = np.nan
    line_curr[line_curr[:,0] > vol_xmax,0] = np.nan
    line_curr[line_curr[:,1] > vol_ymax,0] = np.nan
    line_curr[line_curr[:,2] > vol_zmax,0] = np.nan
    
    # get gradient values along line
    grad_curr = linear_interpolation3d(line_curr[:,0],line_curr[:,1],line_curr[:,2],grad_array)
    
    # shot line plot
    if show_plot:
        plt.plot(grad_curr)
        plt.xlabel("WM -> GM")
        plt.ylabel("Second order gradient")
    
    # get point of zero crossing
    i = 0
    while i <= n_line - 2:        
        
        if grad_curr[i] < 0 and grad_curr[i+1] > 0 and t2s:
            zero_curr = apply_affine(vox2ras, line_curr[i,:])
            shift_curr = vtx_curr - zero_curr
            break
        elif grad_curr[i] > 0 and grad_curr[i+1] < 0 and not t2s:
            zero_curr = apply_affine(vox2ras, line_curr[i,:])
            shift_curr = vtx_curr - zero_curr
            break
        
        i += 1
    
    return shift_curr