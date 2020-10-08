# -*- coding: utf-8 -*-

# python standard library inputs
import os
import sys

# external inputs
import numpy as np
import matplotlib.pyplot as plt


def cost_plot(x_max, cost_array, m, n, set_title="", save_plot=False, 
              path_output="", name_output=""):
    """
    This function plots the cost function array with a corresponding linear fit.
    Inputs:
        *x_max (float): maximum x-coordinate.
        *cost_array (arr): array with cost function values J.
        *m (float): slope of linear fit.
        *n (float): y-axis intercept of linear fit.
        *set_title (str): plot title.
        *save_plot (bool): write out image file.
        *path_output (str): path where output is written.
        *name_output (str): basename of saved plot.

    created by Daniel Haenelt
    Date created: 24-02-2020         
    Last modified: 08-10-2020
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
    