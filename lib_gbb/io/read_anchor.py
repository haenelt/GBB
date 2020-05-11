def read_anchor(input_anchor):
    """
    This function reads a freesurfer point set file and returns the points in a numpy array.
    Inputs:
        *input_anchor: filename of the point set textfile.
    Outputs:
        *data: numpy array with point coordinates.
        
    created by Daniel Haenelt
    Date created: 11-05-2020
    Last modified: 11-05-2020
    """
    import numpy as np

    with open(input_anchor, "r") as f:
        x = f.readlines()
        
        # skip last 3 lines
        n_points = len(x) - 3
        
        # get anchor coordinates into numpy array
        data = np.zeros((n_points,3))
        for i in range(n_points):
            line = x[i].split()
        
            data[i,0] = line[0]
            data[i,1] = line[1]
            data[i,2] = line[2]
    
    return data