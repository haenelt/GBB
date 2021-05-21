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

__all__ = ['SphereMesh', 'CubeMesh', 'SphereVolume', 'CubeVolume',
           'download_data']


class SphereMesh:
    """
    Triangle mesh of icosahedron.

    Computation of a spherical mesh based on an icosahedron as primitive. The
    icosahedron consists of 20 equilateral triangles. The primitive can be
    further refined to a sphere by subdivision of triangles. The center point of
    the sphere is in the origin of the coordinate system. The code for
    subdivision is largely taken from [1]_.

    Parameters
    ----------
    scale : float
        Radius of sphere in px.

    Attributes
    ----------
    subdiv : int
        Number of triangle subdivisions.
    arr_f : (N,N,N) np.ndarray
        Noise array.
    middle_point_cache : dict
        Cache to prevent creation of duplicated vertices.
    vtx : (N,3) np.ndarray
        Array of vertex coordinates.
    fac : (N,3) np.ndarray
        Array of corresponding faces.

    References
    ----------
    .. [1] https://sinestesia.co/blog/tutorials/python-icospheres/

    """

    def __init__(self, scale):

        self.scale = scale
        self.subdiv = 0
        self.arr_f = None
        self.middle_point_cache = {}

        t = (1 + np.sqrt(5)) / 2  # golden ratio
        self.vtx = np.array([self._vertex(-1, t, 0),
                             self._vertex(1, t, 0),
                             self._vertex(-1, -t, 0),
                             self._vertex(1, -t, 0),

                             self._vertex(0, -1, t),
                             self._vertex(0, 1, t),
                             self._vertex(0, -1, -t),
                             self._vertex(0, 1, -t),

                             self._vertex(t, 0, -1),
                             self._vertex(t, 0, 1),
                             self._vertex(-t, 0, -1),
                             self._vertex(-t, 0, 1),
                             ])

        self.fac = np.array([
            # 5 faces around point 0
            [0, 11, 5],
            [0, 5, 1],
            [0, 1, 7],
            [0, 7, 10],
            [0, 10, 11],

            # adjacent faces
            [1, 5, 9],
            [5, 11, 4],
            [11, 10, 2],
            [10, 7, 6],
            [7, 1, 8],

            # 5 faces around 3
            [3, 9, 4],
            [3, 4, 2],
            [3, 2, 6],
            [3, 6, 8],
            [3, 8, 9],

            # adjacent faces
            [4, 9, 5],
            [2, 4, 11],
            [6, 2, 10],
            [8, 6, 7],
            [9, 8, 1],
        ])

    def _vertex(self, x, y, z):
        """Normalize vertex coordinates and scale.

        Parameters
        ----------
        x : float
            x-coordinate.
        y : float
            y-coordinate.
        z : float
            z-coordinate.

        Returns
        -------
        (3,) list
            Scaled coordinates.

        """

        length = np.sqrt(x ** 2 + y ** 2 + z ** 2)

        return [(x * self.scale) / length for x in (x, y, z)]

    def _middle_point(self, ind1, ind2):
        """Find a middle point between two vertices and project to unit sphere.

        Parameters
        ----------
        ind1 : int
            Index of vertex 1.
        ind2 : int
            Index of vertex 2.

        Returns
        -------
        index : int
            Index of created middle point.

        """

        # We check if we have already cut this edge first to avoid duplicated
        # vertices
        smaller_index = min(ind1, ind2)
        greater_index = max(ind1, ind2)

        key = '{0}-{1}'.format(smaller_index, greater_index)

        if key in self.middle_point_cache:
            return self.middle_point_cache[key]

        # If it's not in cache, then we can cut it
        vert_1 = self.vtx[ind1, :]
        vert_2 = self.vtx[ind2, :]

        middle = np.mean([vert_1, vert_2], axis=0)
        self.vtx = np.vstack((self.vtx, self._vertex(*middle)))

        index = len(self.vtx) - 1
        self.middle_point_cache[key] = index

        return index

    def subdivide(self, subdiv):
        """Subdivide icosahedron by splitting each triangle into 4 smaller
        triangles.

        Parameters
        ----------
        subdiv : int
            Number of subdivision iterations.

        """

        self.subdiv += subdiv
        self.middle_point_cache = {}

        for i in range(self.subdiv):
            faces_subdiv = []

            for tri in self.fac:
                v1 = self._middle_point(tri[0], tri[1])
                v2 = self._middle_point(tri[1], tri[2])
                v3 = self._middle_point(tri[2], tri[0])

                faces_subdiv.append([tri[0], v1, v3])
                faces_subdiv.append([tri[1], v2, v1])
                faces_subdiv.append([tri[2], v3, v2])
                faces_subdiv.append([v1, v2, v3])

            self.fac = np.array(faces_subdiv)

    def add_noise(self, amplitude=10, rho=0.001, sigma=0.05, negative=True,
                  arr_noise=None):
        """Add noise to vertex coordinates.

        A deformation field is created by applying a gaussian bandpass filter
        to white noise. The mesh is then shifted in normal direction according
        to the deformation field.

        Parameters
        ----------
        amplitude : float, optional
            Amplitude of deformation field.
        rho : float, optional
            Spatial center frequency of gaussian bandpass filter.
        sigma : float, optional
            Standard deviation in cycles/px of gaussian bandpass filter.
        negative : bool, optional
            Only allow deformations in outward direction.
        arr_noise : (N,N,N) np.ndarray, optional
            Use an already existing white noise array instead of computing a new
            one.

        Raises
        ------
            ValueError
                If `arr_noise` is used as argument and has not the right shape.

        """

        # define k-space
        k_max = MATRIX_SIZE / (2 * MATRIX_SIZE * VOXEL_RES)
        k = np.linspace(-k_max, k_max, MATRIX_SIZE)
        kx_grid, ky_grid, kz_grid = np.meshgrid(k, k, k, indexing='ij')
        arr_kr = np.sqrt(kx_grid ** 2 + ky_grid ** 2 + kz_grid ** 2)

        # gaussian filter
        gaussian = np.exp(-(arr_kr - rho) ** 2 / (2 * sigma ** 2))
        gaussian = fftshift(gaussian)

        # gaussian noise
        correct_shape = np.shape(arr_noise) == (
            MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE)
        if arr_noise and correct_shape:
            pass
        elif arr_noise and not correct_shape:
            raise ValueError("Passed noise array is not of the right shape.")
        else:
            arr_noise = normal(0, 1, MATRIX_SIZE ** 3)
            arr_noise = np.reshape(arr_noise,
                                   (MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE))

        self.arr_f = np.real(ifftn(fftn(arr_noise) * gaussian))
        self.arr_f /= np.max(self.arr_f)

        if not negative:
            self.arr_f[self.arr_f < 0] = 0

        # sample shifts from deformation map
        n = get_normal(self.vtx, self.fac)
        vtx_shift = linear_interpolation3d(self.vtx[:, 0],
                                           self.vtx[:, 1],
                                           self.vtx[:, 2],
                                           self.arr_f)

        self.vtx = self.vtx + amplitude * n * vtx_shift[:, np.newaxis]

    def save_img(self, file_out, rotation=(0, 0, 0)):
        """Save image of mesh to disk.

        Parameters
        ----------
        file_out : str
            Filename of saved image file.
        rotation : (3,) tuple, optional
            Apply rotation in (x, y, z) directions to the mesh before saving the
            image.
        """

        dir_out = os.path.dirname(file_out)
        if not os.path.exists(dir_out):
            os.makedirs(dir_out)

        fac_conv = np.zeros((len(self.fac), 4), dtype=int)
        fac_conv[:, 0] = 3
        fac_conv[:, 1:] = self.fac

        polygon = pv.PolyData()
        polygon.points = self.vtx
        polygon.faces = fac_conv

        # rotate mesh
        polygon.rotate_x(rotation[0])
        polygon.rotate_y(rotation[1])
        polygon.rotate_z(rotation[2])

        # create plot
        plotter = pv.Plotter(off_screen=True, window_size=[1024, 1024])
        plotter.add_mesh(polygon, show_edges=True, color="tan")
        plotter.show(screenshot=file_out)

    def save(self, file_out):
        """Save mesh to disk.

        Parameters
        ----------
        file_out : str
            Filename of mesh file in freesurfer format.

        """

        dir_out = os.path.dirname(file_out)
        if not os.path.exists(dir_out):
            os.makedirs(dir_out)

        write_geometry(file_out, self.vtx, self.fac)

    @property
    def scale(self):
        return self._scale

    @scale.setter
    def scale(self, s):
        self._scale = s


class CubeMesh(SphereMesh):
    """
    Triangle mesh of cube.

    Computation of a cube primitive which consists of 12 equilateral triangles.
    The primitive can be further refined by subdivision of triangles. The center
    point of the cube is in the origin of the coordinate system. The code for
    subdivision is largely taken from [1]_.

    Parameters
    ----------
    scale : float
        Half edge length of cube in px.

    Attributes
    ----------
    vtx : (N,3) np.ndarray
        Array of vertex coordinates.
    fac : (N,3) np.ndarray
        Array of corresponding faces.

    References
    ----------
    .. [1] https://sinestesia.co/blog/tutorials/python-icospheres/

    """

    def __init__(self, scale):
        super().__init__(scale)

        self.vtx = np.array([[1.0, -1.0, 1.0],
                             [1.0, -1.0, -1.0],
                             [1.0, 1.0, -1.0],
                             [1.0, 1.0, 1.0],
                             [-1.0, -1.0, 1.0],
                             [-1.0, -1.0, -1.0],
                             [-1.0, 1.0, -1.0],
                             [-1.0, 1.0, 1.0],
                             ])
        self.vtx *= self.scale

        self.fac = np.array([[4, 0, 3],
                             [4, 3, 7],
                             [0, 1, 2],
                             [0, 2, 3],
                             [1, 5, 6],
                             [1, 6, 2],
                             [5, 4, 7],
                             [5, 7, 6],
                             [7, 3, 2],
                             [7, 2, 6],
                             [0, 5, 1],
                             [0, 4, 5],
                             ])

    def _middle_point(self, ind1, ind2):
        """Find a middle point between two vertices.

        Parameters
        ----------
        ind1 : int
            Index of vertex 1.
        ind2 : int
            Index of vertex 2.

        Returns
        -------
        index : int
            Index of created middle point.

        """

        # We check if we have already cut this edge first
        # to avoid duplicated verts
        smaller_index = min(ind1, ind2)
        greater_index = max(ind1, ind2)

        key = '{0}-{1}'.format(smaller_index, greater_index)

        if key in self.middle_point_cache:
            return self.middle_point_cache[key]

        # If it's not in cache, then we can cut it
        vert_1 = self.vtx[ind1, :]
        vert_2 = self.vtx[ind2, :]

        middle = np.mean([vert_1, vert_2], axis=0)
        self.vtx = np.vstack((self.vtx, middle))

        index = len(self.vtx) - 1
        self.middle_point_cache[key] = index

        return index


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


def download_data(dir_out):
    """ Download data.

    This function downloads the example data from my personal google drive [1]_
    which enables you to directly execute the code of this repository. Data is
    taken from a human 7 T fMRI study. Further information about the dataset
    can be found in the provided documentation.

    Parameters
    ----------
    dir_out : str
        Directory in which downloaded files are stored.

    Raises
    ------
    FileExistsError
        If `dir_out` already exists.

    Returns
    -------
    dict
        Dictionary collecting the output under the following keys

        * control_points (str) : Path to control points.
        * ignore (str) : Path to ignore mask.
        * mean_epi_enhanced (str) : Path to mean epi with enhanced contrast.
        * mean_epi (str) : Path to mean epi.
        * vein (str) : Path to vein mask.
        * pial (str) : path to pial surface.
        * white (str) : path to white surface.

    References
    ----------
    .. [1] https://drive.google.com/drive/folders/1dwra3FDRIQ-W3iO81vfxvYRHKqGgqSWn

    """

    # make output folder
    if not os.path.exists(dir_out):
        os.makedirs(dir_out)

    # base url to google drive folder
    url = 'https://drive.google.com/uc?export=download&id='

    # file IDs of source files
    file_id = ['1Ljotc8rgbJ33Uial708OSCPA5jzosGiR',
               '1w-iZlmKk6133uLxpl38LoVOkJel2kpWl',
               '1x1KO_ZQsTiht0IERLck4uNDRweTFiQYy',
               '1r6yUQpxmnALp6MWfB-pCvytB11off5gN',
               '1j02nHLrLX-sO5NnU5QoNt39sytdhPprh',
               '1ZY2Ld8VXocSfbxJIxIqm9hX1rBDB_99X',
               '1us0lNqNZe-YuO8SodsLbSe5gOmUE52cU',
               '1_FR17jOaEZpCtYL36nG0TMe9SYhK52Dg',
               ]

    # file names of output files
    filename = ['control_points.dat',
                'ignore.nii',
                'mean_epi_enhanced.nii',
                'mean_epi.nii',
                'vein.nii',
                'lh.pial',
                'lh.white',
                'Readme.md',
                ]

    file_sources = [url + id_ for id_ in file_id]
    file_targets = [os.path.join(dir_out, file) for file in filename]

    for source, target in zip(file_sources, file_targets):

        # check if file already exists
        if os.path.exists(target):
            raise FileExistsError("The file " + target + " already exists!")
        else:
            gdown.download(source, target, quiet=False)

    return dict(control_points=file_targets[0],
                ignore=file_targets[1],
                mean_epi_enhanced=file_targets[2],
                mean_epi=file_targets[3],
                vein=file_targets[4],
                pial=file_targets[5],
                white=file_targets[6])
