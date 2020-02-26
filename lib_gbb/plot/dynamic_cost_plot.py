def dynamic_cost_plot(file_in, path_output, name_output):
    """
    This function creates an animation to illustrate the progression of the cost function with 
    increasing number of iterations.
    Inputs:
        *file_in: filename of the saved compressed numpy variables.
        *path_output: path where output is written.
        *name_output: basename of saved animation.
        
    created by Daniel Haenelt
    Date created: 24-02-2020     
    Last modified: 24-02-2020
    """
    import os
    import glob
    import numpy as np
    from lib_gbb.plot.cost_plot import cost_plot
    from lib.img.get_gif import get_gif

    # load numpy variables
    data = np.load(file_in)
    
    # get cost array, slop and y-axis intercept
    J = data["J"]
    m = data["m"]
    n = data["n"]
    
    # save single frames
    i = len(J) - len(m)
    while i < len(J):
        cost_plot(i, 
                  J[:i], 
                  m[i - (len(J) - len(m))], 
                  n[i - (len(J) - len(m))], 
                  set_title=None,
                  save_plot=True, 
                  path_output=path_output,
                  name_output="%04d" % i)
        
        i += 1
       
    # get all single frames
    img_file = glob.glob(os.path.join(path_output,"*.png"))
    img_file = sorted(img_file)
    
    # write gif
    get_gif(img_file, path_output, name_output, 1, 0.25)