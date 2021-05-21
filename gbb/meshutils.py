# -*- coding: utf-8 -*-
"""Generation of test objects (sphere, cube) and download of test data."""

# python standard library inputs
import os

# external inputs
import numpy as np
import pyvista as pv
import nibabel as nb
import gdown
from numpy.fft import fftn, ifftn, fftshift
from numpy.random import normal
from nibabel.freesurfer.io import write_geometry

# local inputs
from config import MATRIX_SIZE, VOXEL_RES
from gbb.normal import get_normal
from interpolation import linear_interpolation3d

__all__ = ['Mesh']


class Mesh:
    pass

    # read surf or filename
    # get normal
    # plot_normal
    # plot_normal_direction
    # adjm
    # remove vertex
    # smooth surface (laplacian)
    # remesh