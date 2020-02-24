def slope_plot(x_max, m_array, set_title=False, save_plot=False, path_output=False, 
               name_output=False):
    """
    This function plots the slopes from the computed linear fits.
    Inputs:
        *x_max: maximum x-coordinate.
        *m_array: slope array.
        *set_title: optional plot title.
        *save_plot: write out image file.
        *path_output: path where output is written.
        *name_output: basename of saved plot.

    created by Daniel Haenelt
    Date created: 24-02-2020         
    Last modified: 24-02-2020
    """
    import os
    import numpy as np
    import matplotlib.pyplot as plt

    # make output folder
    if not os.path.exists(path_output):
        os.makedirs(path_output)
    
    # show plot
    plt.plot(np.arange(x_max), m_array, label="m: "+str("%.10e" % m_array[-1]))
    plt.legend(loc=1)
    plt.xlabel("iteration")
    plt.ylabel("Cost function")

    # make title
    if set_title:
        plt.title(set_title)   
    
    # save figure
    if save_plot:
        plt.savefig(os.path.join(path_output,name_output+".png"))

    plt.pause(0.01)