# -*- coding: utf-8 -*-

# python standard library inputs
import os
import sys

# external inputs
import numpy as np
import matplotlib.pyplot as plt


def plot_cost(x_max, cost_array, m, n, set_title="", save_plot=False, 
              path_output="", name_output=""):
    """ Plot cost

    This function plots the cost function array with a corresponding linear fit.    

    Parameters
    ----------
    x_max : float
        Maximum x-coordinate.
    cost_array : ndarray
        Array with cost function values J.
    m : float
        Slope of linear fit.
    n : float
        y-axis intercept of linear fit.
    set_title : str, optional
        Plot title. The default is "".
    save_plot : bool, optional
        Write out image file. The default is False.
    path_output : str, optional
        Path where output is written. The default is "".
    name_output : str, optional
        Basename of saved plot. The default is "".

    Returns
    -------
    None.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 24-02-2020         
    Last modified: 15-10-2020
    
    """
    
    # make output folder    
    try:
        if not os.path.exists(path_output):
            os.makedirs(path_output)
    except TypeError:
        sys.exit("error: Output directory not defined!")
    
    # compute line
    line_fit = np.arange(x_max) * m + n
    
    # show plot
    plt.clf()
    plt.plot(np.arange(x_max), cost_array, label="J")
    plt.plot(np.arange(x_max), line_fit, label="fit (m: "+str("%.10e" % m)+")")
    plt.legend(loc=1)
    plt.xlabel("iteration")
    plt.ylabel("cost function")
    
    # make title
    if set_title:
        plt.title(set_title)   
    
    # save figure
    if save_plot:
        plt.savefig(os.path.join(path_output,name_output+".png"))

    plt.pause(0.01)
    