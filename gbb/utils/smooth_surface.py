# -*- coding: utf-8 -*-

# python standard library inputs
import os
import sys
import subprocess

# local inputs
from gbb.io.get_filename import get_filename


def smooth_surface(file_in, file_out, n_iter):
    """ Smooth surface
    
    This function smoothes a surface mesh using freesurfer.

    Parameters
    ----------
    file_in : str
        Filename of input surface.
    file_out : str
        Filename of output surface.
    n_iter : int
        Number of smoothing iterations.

    Returns
    -------
    None.

    Notes
    -------
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
        