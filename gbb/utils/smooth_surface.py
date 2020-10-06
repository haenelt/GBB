# -*- coding: utf-8 -*-

# python standard library inputs
import os
import sys
import subprocess

# local inputs
from gbb.io.get_filename import get_filename


def smooth_surface(file_in, file_out, n_iter):
    """
    This function smoothes a surface mesh using freesurfer.
    Inputs:
        *file_in (str): filename of input surface.
        *file_out (str): filename of output surface.
        *n_iter (int): number of smoothing iterations.
        
    created by Daniel Haenelt
    Date created: 13-07-2019
    Last modified: 05-10-2020
    """
       
    # make output folder
    path_output, _, _ = get_filename(file_out)
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # smooth surface
    try:
        subprocess.run(['mris_smooth', 
                        '-n', str(n_iter), 
                        '-nw', 
                        file_in, 
                        file_out], check = True)
    except subprocess.CalledProcessError:
        sys.exit("error: surface smoothing failed!")
        