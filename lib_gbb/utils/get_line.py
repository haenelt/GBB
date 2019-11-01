def get_line(vox_coords, input_array, line_dir, line_size, step_size, interpolation="linear", show_plot=False):
    """
    This function computes an interpolated line in a specific direction within an input array. The 
    line final line is interpolated using spline interpolation.
    Inputs:
        *vox_coords: [x,y,z] coordinates of starting point in input array.
        *input_array: 3d nifti array from which line values are computed.
        *line_dir: direction of line either horizontally (0) or vertically (1), not through slices.
        *line_size: line length 2*line_size in mm.
        *step_size: distance between sampled line points in mm.
        *interpolation: interpolation method for line computation (either linear or nearest).
        *show_plot: show plot of final line values.
    Outputs:
        *line_coords_new: array of line coordinates.
        *line_smooth: array of line values.
        
    created by Daniel Haenelt
    Date created: 31-10-2019 
    Last modified: 31-10-2019
    """
    import numpy as np
    from scipy.interpolate import splev, splrep
    from lib_gbb.interpolation import linear_interpolation2d, nn_interpolation2d
    import matplotlib.pyplot as plt
    
    # get vox coordinates of line
    line_coords = np.arange(vox_coords[line_dir] - line_size, vox_coords[line_dir] + line_size, step_size)
    
    # exclude coordinates at volume edges
    line_coords = line_coords[line_coords > 0]
    line_coords = line_coords[line_coords < np.shape(input_array)[line_dir] -1] 

    line = []
    z = np.round(vox_coords[2]).astype(int) # slice
    for i in range(len(line_coords)):
    
        # set coordinates for chosen line direction
        if line_dir == 0:
            x = line_coords[i]
            y = vox_coords[1]
        elif line_dir == 1:
            x = vox_coords[0]
            y = line_coords[i]
        else:
            print("choose a valid line direction!")
            
        # get interpolated line point
        if interpolation == "linear":
            line.append(linear_interpolation2d(x, y, input_array[:,:,z]))
        elif interpolation == "nearest":
            line.append(nn_interpolation2d(x, y, input_array[:,:,z]))
        else:
            print("choose a valid interpolation method!")      

    # spline interpolation of line data
    line_coords_new = np.linspace(line_coords.min(), line_coords.max(), 1000)
    line_spl = splrep(line_coords, line)
    line_smooth = splev(line_coords_new, line_spl)
    
    if show_plot:
        plt.plot(line_coords, line, 'o', line_coords_new, line_smooth)
        plt.show()
    
    return line_coords_new, line_smooth