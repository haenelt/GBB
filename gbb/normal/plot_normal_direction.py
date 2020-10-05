# -*- coding: utf-8 -*-

# external inputs
import numpy as np
import nibabel as nb
from nibabel.freesurfer.io import read_geometry

# local inputs
from gbb.normal import get_normal
    

def plot_normal_direction(input_surf, axis=2):
    """
    This function plots the direction from the white surface towards the pial 
    surface based on the white surface vertex normals along one axis.
    Inputs:
        *input_surf (str): filename of source mesh (white surface).
        *axis (int): axis for distance calculation in ras space (0,1,2).
        
    created by Daniel Haenelt
    Date created: 13-12-2019         
    Last modified: 05-10-2020
    """

    # fixed parameter
    line_threshold = 0.05 # if direction is along one axis, omit line if length is below threshold

    # load geometry
    vtx_source, fac_source = read_geometry(input_surf)

    # get surface normals per vertex
    norm = get_normal(vtx_source, fac_source)
    
    # get distance along one axis
    r_dist = norm[:,axis].copy()
    
    # get directions
    r_dist[r_dist > line_threshold] = 1
    r_dist[r_dist < -line_threshold] = -1
    r_dist[np.abs(r_dist) != 1] = 0

    # write output
    header = nb.freesurfer.mghformat.MGHHeader()
    output = nb.freesurfer.mghformat.MGHImage(r_dist, np.eye(4), header)
    nb.save(output, input_surf+"_plot_normal_dir"+str(axis)+".mgh")
    