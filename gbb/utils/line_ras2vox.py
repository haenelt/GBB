# -*- coding: utf-8 -*-

# python standard library inputs
import sys

# external inputs
import numpy as np
from nibabel.affines import apply_affine
    
    
def line_ras2vox(line_dir, ras2vox_tkr):
    """
    This function transforms the line direction from ras coordinates into its 
    orientation in voxel space.
    Inputs:
        *line_dir (int): line direction in ras space (0,1,2).
        *ras2vox_tkr (arr): transformation matrix from ras to voxel space.
    Outputs:
        *line_vox (int): line direction in voxel space.
        
    created by Daniel Haenelt
    Date created: 28-12-2019
    Last modified: 05-10-2020
    """

    # get line direction in ras space
    if line_dir == 0:
        pt_ras = [1,0,0]
    elif line_dir == 1:
        pt_ras = [0,1,0]
    elif line_dir == 2:
        pt_ras = [0,0,1]
    else:
        sys.exit("error: invalid axis direction in gradient calculation!")

    # get unit vector in voxel space
    pt0_vox = apply_affine(ras2vox_tkr, [0,0,0])
    pt_vox = apply_affine(ras2vox_tkr, pt_ras)

    pt_vox = np.abs( pt_vox - pt0_vox )
    pt_vox = pt_vox / np.max(pt_vox)
    pt_vox = np.round(pt_vox).astype(int)

    # updated line direction in voxel space
    line_vox = np.where(pt_vox == 1)[0][0]

    return line_vox
