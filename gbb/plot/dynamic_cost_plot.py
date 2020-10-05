# -*- coding: utf-8 -*-

# python standard library inputs
import os
import glob

# external inputs
import numpy as np

# local inputs
from gbb.plot.cost_plot import cost_plot
from gbb.plot.get_gif import get_gif


def dynamic_cost_plot(file_in, path_output, name_output):
    """
    This function creates an animation to illustrate the progression of the cost 
    function with increasing number of iterations.
    Inputs:
        *file_in (str): filename of the saved compressed numpy variables.
        *path_output (str): path where output is written.
        *name_output (str): basename of saved animation.
        
    created by Daniel Haenelt
    Date created: 24-02-2020     
    Last modified: 05-10-2020
    """

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
    