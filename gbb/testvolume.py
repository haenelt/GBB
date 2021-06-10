# -*- coding: utf-8 -*-
"""Generation of test volumes (sphere, cube)."""

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb
from numpy.fft import fftn, ifftn, fftshift
from numpy.random import normal

# local inputs
from .config import MATRIX_SIZE, VOXEL_RES

__all__ = ['SphereVolume', 'CubeVolume']


class SphereVolume:
    """
    Brain with spherical shape.

    3D array of spherical test object containing three compartment (WM: 1, GM:
    2, CSF: 3). The object is placed in the center of the array. White noise and
    a lowpass filter can be added to test various noise levels.

    Parameters
    ----------
    scale : float
        Radius of WM segment in px.
    cortical_thickness : float
        Thickness of GM segment in px.

    Attributes
    ----------
    arr_brain : (N,N,N) np.ndarray
        Array of test object.

    """

    def __init__(self, scale, cortical_thickness):

        self.scale = scale
        self.cortical_thickness = cortical_thickness

        d = np.arange(-MATRIX_SIZE/2, MATRIX_SIZE/2, dtype=float)
        x, y, z = np.meshgrid(d, d, d, indexing='ij')
        arr_dist = np.sqrt(x ** 2 + y ** 2 + z ** 2)

        self.arr_brain = np.zeros_like(arr_dist) + 3  # CSF
        self.arr_brain[arr_dist < scale + cortical_thickness] = 2  # GM
        self.arr_brain[arr_dist < scale] = 1  # WM

    def add_noise(self, sigma=0.5, arr_noise=None):
        """Add white noise to array.

        Parameters
        ----------
        sigma : float, optional
            Standard deviation of noise in px.
        arr_noise : (N,N,N), optional
            Use an already existing white noise array instead of computing a new
            one.

        Raises
        ------
            ValueError
                If `arr_noise` is used as argument and has not the right shape.

        """

        correct_shape = np.shape(arr_noise) == (
            MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE)
        if arr_noise and correct_shape:
            pass
        elif arr_noise and not correct_shape:
            raise ValueError("Passed noise array is not of the right shape.")
        else:
            arr_noise = normal(0, sigma, MATRIX_SIZE ** 3)
            arr_noise = np.reshape(arr_noise,
                                   (MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE))

        self.arr_brain += arr_noise

    def smooth(self, sigma=0.2):
        """Apply lowpass filter to data.

        Parameters
        ----------
        sigma : float, optional
            Standard deviation of noise in cycles/px.

        """

        # define k-space
        k_max = MATRIX_SIZE / (2 * MATRIX_SIZE * VOXEL_RES)
        k = np.linspace(-k_max, k_max, MATRIX_SIZE)
        kx_grid, ky_grid, kz_grid = np.meshgrid(k, k, k, indexing='ij')
        arr_kr = np.sqrt(kx_grid ** 2 + ky_grid ** 2 + kz_grid ** 2)

        # gaussian filter
        beta = 1
        arr_f = beta * np.exp(-arr_kr ** 2 / (2 * sigma ** 2))
        arr_f = fftshift(arr_f)

        self.arr_brain = np.abs(ifftn(fftn(self.arr_brain) * arr_f))

    def save(self, file_out):
        """Save volume to disk.

        Parameters
        ----------
        file_out : str
            Filename of volume array in nifti format.

        """

        dir_out = os.path.dirname(file_out)
        if not os.path.exists(dir_out):
            os.makedirs(dir_out)

        output = nb.Nifti1Image(self.arr_brain, np.eye(4))
        nb.save(output, file_out)

    @property
    def scale(self):
        return self._scale

    @scale.setter
    def scale(self, s):
        self._scale = s

    @property
    def cortical_thickness(self):
        return self._cortical_thickness

    @cortical_thickness.setter
    def cortical_thickness(self, c):
        self._cortical_thickness = c


class CubeVolume(SphereVolume):
    """
    Brain with cubic shape.

    3D array of cubic test object containing three compartment (WM: 1, GM: 2,
    CSF: 3). The object is placed in the center of the array. White noise and a
    lowpass filter can be added to test various noise levels.

    Parameters
    ----------
    scale : float
        Half edge length of WM segment in px.
    cortical_thickness : float
        Thickness of GM segment in px.

    Attributes
    ----------
    arr_brain : (N,N,N) np.ndarray
        Array of test object.

    """

    def __init__(self, scale, cortical_thickness):
        super().__init__(scale, cortical_thickness)

        d = np.arange(-MATRIX_SIZE/2, MATRIX_SIZE/2, dtype=float)
        x, y, z = np.abs(np.meshgrid(d, d, d, indexing='ij'))

        wm = np.ones_like(x)
        wm[x > self.scale] = 0
        wm[y > self.scale] = 0
        wm[z > self.scale] = 0

        gm = np.ones_like(x)
        gm[x > self.scale + self.cortical_thickness] = 0
        gm[y > self.scale + self.cortical_thickness] = 0
        gm[z > self.scale + self.cortical_thickness] = 0

        self.arr_brain = np.zeros_like(x) + 3
        self.arr_brain[gm == 1] = 2
        self.arr_brain[wm == 1] = 1
