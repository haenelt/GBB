# -*- coding: utf-8 -*-
"""Generation of test objects (sphere, cube) and download of test data."""

# python standard library inputs
import os
import itertools

# external inputs
import numpy as np
import pyvista as pv
import nibabel as nb
from scipy.sparse import csr_matrix
from numpy.fft import fftn, ifftn, fftshift
from numpy.random import normal
from nibabel.freesurfer.io import write_geometry

# local inputs
from .config import MATRIX_SIZE, VOXEL_RES
from .interpolation import linear_interpolation3d

__all__ = ['Mesh', 'SphereMesh', 'CubeMesh']

# TODO: inherit this parent class in SphereMesh and CubeMesh?!?
# TODO: http://graphics.stanford.edu/courses/cs468-12-spring/LectureSlides/13_Remeshing1.pdf
"""
Remeshing

- consists of edge collapse, edge split, edge flip, vertex shift
- we only consider vertex shifts since we want to preserve topology
- vertex shift -> local "spring" relaxation (uniform Laplacian smoothing)
- bary-center of one-ring neighbor:
    c_i = 1 / valence(v_i) \sum_{j\in N(v_i)} p_j
- keep vertex (approximately) on surface
    p_i \leftarrow c_i - n_in_i^T(c_i-p_i)
- at the moment, no edge will be preserved since brain surface is smooth
- projection derivation
    p = \frac{a\cdot b}{a\cdot a}
      = \frac{1}{a\cdot a}a(a\cdot b)
      = \frac{1}{a^Ta}a(a^Tb)
      = \frac{1}{a^Ta}(aa^T)b
"""


class Mesh:
    """
    # save_normals
    # save_normal_dir
    # save_img
    # adjm
    # vertex_normals
    # smooth surface (laplacian)
    # save
    # remove vertex
    # remesh

    """

    def __init__(self, vtx, fac):

        self.vtx = vtx
        self.fac = fac

    @property
    def adjm(self):
        """Get adjm.
    
        This function computes a sparse adjacency matrix for a triangular surface
        mesh. The matrix has the size (nvertex,nvertex). Each matrix entry with 
        value 1 stands for an edge of the surface mesh.
    
        Returns
        -------
        sparse_adjm : obj
            Sparse adjacency matrix.
        
        """

        # get number of vertices and faces
        nvtx = len(self.vtx)
        nfac = len(self.fac)

        # initialise
        row = []
        col = []

        # get rows and columns of edges
        row.extend([self.fac[i, 0] for i in range(nfac)])
        col.extend([self.fac[i, 1] for i in range(nfac)])

        row.extend([self.fac[i, 1] for i in range(nfac)])
        col.extend([self.fac[i, 2] for i in range(nfac)])

        row.extend([self.fac[i, 2] for i in range(nfac)])
        col.extend([self.fac[i, 0] for i in range(nfac)])

        # make sure that all edges are symmetric
        row.extend([self.fac[i, 1] for i in range(nfac)])
        col.extend([self.fac[i, 0] for i in range(nfac)])

        row.extend([self.fac[i, 2] for i in range(nfac)])
        col.extend([self.fac[i, 1] for i in range(nfac)])

        row.extend([self.fac[i, 0] for i in range(nfac)])
        col.extend([self.fac[i, 2] for i in range(nfac)])

        # adjacency entries get value 1
        data = np.ones(len(row), dtype=np.int8)

        # write sparse adjacency matrix
        sparse_adjm = csr_matrix((data, (row, col)), shape=(nvtx, nvtx))

        return sparse_adjm

    @property
    def vertex_normals(self):
        """ Get normal

        This function computes the surfaces normals per vertex from an input surface
        mesh. The code is taken from [1] and adapted to my own purposes.

        Returns
        -------
        norm : ndarray
            Normal per vertex.

        References
        -------
        .. [1] https://sites.google.com/site/dlampetest/python/calculating-normals-
        of-a-triangle-mesh-using-numpy

        Notes
        -------
        created by Daniel Haenelt
        Date created: 13-12-2019
        Last modified: 05-10-2020

        """

        # indexed view into the vertex array
        tris = self.vtx[self.fac]

        # calculate the normal for all triangles by taking the cross product of
        # the vectors v1-v0 and v2-v0 in each triangle
        n = np.cross(tris[::, 1] - tris[::, 0], tris[::, 2] - tris[::, 0])

        # calculate vertex-wise normals and normalize
        normal = self._f2v(len(self.vtx), self.fac, n)
        normal = self._normalize_v3(normal)

        return normal

    def remove_vertex(self, ind_keep):
        """Remove vertex

        This function removes vertices from an array of vertex points and updates
        the corresponding face array.

        Parameters
        ----------
        ind_keep : list
            Index list of vertices to keep.

        Returns
        -------
        vtx : ndarray
            Remaining vertices.
        fac : ndarray
            Remaining faces.
        ind_keep : ndarray
            Updated index list of vertices to keep after vertex cleaning.

        """

        # get indices which will be removed
        ind_tmp = np.arange(len(self.vtx))
        ind_remove = list(set(ind_tmp) - set(ind_keep))
        ind_remove = sorted(ind_remove, reverse=True)

        # get new vertices
        self.vtx = self.vtx[ind_keep, :]

        # get new faces
        fac_keep = np.zeros(len(self.fac))
        fac_keep += np.in1d(self.fac[:, 0], ind_keep)
        fac_keep += np.in1d(self.fac[:, 1], ind_keep)
        fac_keep += np.in1d(self.fac[:, 2], ind_keep)
        self.fac = self.fac[fac_keep == 3, :]

        # reindex faces
        loop_status = 0
        loop_length = len(ind_remove)
        for i in range(loop_length):

            # print status
            counter = np.floor(i / loop_length * 100)
            if counter != loop_status:
                print("sort faces: " + str(counter) + " %")
                loop_status = counter

            tmp = self.fac[self.fac >= ind_remove[i]] - 1
            self.fac[self.fac >= ind_remove[i]] = tmp

        # get indices which will be cleaned
        ind_vtx = np.arange(len(self.vtx))
        ind_fac = list(itertools.chain(*self.fac))
        ind_fac = list(set(ind_fac))
        ind_remove = list(set(ind_vtx) - set(ind_fac))
        ind_remove = sorted(ind_remove, reverse=True)

        # remove singularities (vertices without faces)
        loop_status = 0
        loop_length = len(ind_remove)
        for i in range(loop_length):

            # print status
            counter = np.floor(i / loop_length * 100)
            if counter != loop_status:
                print("clean faces: " + str(counter) + " %")
                loop_status = counter

            # remove vertex and index
            self.vtx = np.delete(self.vtx, ind_remove[i], 0)
            ind_keep = np.delete(ind_keep, ind_remove[i], 0)

            # sort faces
            tmp = self.fac[self.fac >= ind_remove[i]] - 1
            self.fac[self.fac >= ind_remove[i]] = tmp

        return self.vtx, self.fac

    def _smooth(self, niter):
        vtx_copy = self.vtx.copy()
        adjm_ind = self.adjm.indices
        adjm_ptr = self.adjm.indptr
        for i in range(niter):
            print(i)
            vtx_copy = [
                np.mean(vtx_copy[adjm_ind[adjm_ptr[j]:adjm_ptr[j + 1]], :], 0)
                for j, _ in enumerate(vtx_copy)]
            vtx_copy = np.asarray(vtx_copy)

        return vtx_copy

    def smooth(self, niter):
        return self.smooth(niter), self.fac

    def remesh(self, niter):
        for i in range(niter):
            vtx_smooth = self._smooth(1)
            vtx_diff = vtx_smooth - self.vtx
            n = self.vertex_normals

            self.vtx = [vtx_smooth[j] - np.dot(np.outer(v, v), vtx_diff[j])
                        for j, v in enumerate(n)]

        return self.vtx, self.fac

    def save_normals(self, file_out):
        """ Plot normal

        This function generates lines to visualize outward directed surface
        normals of an input surface mesh.

        Parameters
        ----------
        file_out : str
            Filename of output surface mesh.

        Returns
        -------
        None.

        Notes
        -------
        created by Daniel Haenelt
        Date created: 13-12-2019
        Last modified: 05-10-2020

        """

        # initialise faces for specific shape
        fac_new = [[0, 1, 0]]
        fac_iter = 2

        vtx_res = []
        fac_res = []
        normals = self.vertex_normals
        nvtx = len(self.vtx)

        for i in range(nvtx):
            A = list(self.vtx[i, :])
            B = list(self.vtx[i, :] - normals[i, :])
            vtx_new = [A, B]

            # update faces
            if i > 0:
                for j in range(len(fac_new)):
                    fac_new[j] = [x + fac_iter for x in fac_new[j]]

            # update resulting vertex and face list
            vtx_res.extend(vtx_new)
            fac_res.extend(fac_new)

        # vertices and faces as array
        vtx_res = np.array(vtx_res)
        fac_res = np.array(fac_res)

        # write output geometry
        write_geometry(file_out, vtx_res, fac_res)

    def save_normal_dir(self, file_out, axis=2):
        """ Plot normal direction

        This function plots the direction from the white surface towards the pial
        surface based on the white surface vertex normals along one axis.

        Parameters
        ----------
        file_out : str
            Filename of source mesh (white surface).
        axis : int, optional
            Axis for distance calculation in ras space (0,1,2). The default is 2.

        Returns
        -------
        None.

        """

        # fixed parameter
        line_threshold = 0.05  # if direction is along one axis, omit line if length is below threshold

        # get distance along one axis
        r_dist = self.vertex_normals[:, axis].copy()

        # get directions
        r_dist[r_dist > line_threshold] = 1
        r_dist[r_dist < -line_threshold] = -1
        r_dist[np.abs(r_dist) != 1] = 0

        # write output
        header = nb.freesurfer.mghformat.MGHHeader()
        output = nb.freesurfer.mghformat.MGHImage(r_dist, np.eye(4), header)
        nb.save(output, file_out)

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

    @staticmethod
    def _normalize_v3(arr):
        """Normalize a numpy array of 3 component vectors shape=(n,3) """
        lens = np.sqrt(arr[:, 0] ** 2 + arr[:, 1] ** 2 + arr[:, 2] ** 2)
        res = np.zeros_like(arr)
        res[:, 0] = arr[:, 0] / lens
        res[:, 1] = arr[:, 1] / lens
        res[:, 2] = arr[:, 2] / lens
        return res

    @staticmethod
    def _f2v(length, fac, nt):
        """get average vertex-wise normal by adding up all face-wise normals
        around vertex """
        nv = np.zeros((length, 3))
        for i in range(len(fac)):
            nv[fac[i, 0], :] += nt[i, :]
            nv[fac[i, 1], :] += nt[i, :]
            nv[fac[i, 2], :] += nt[i, :]

        return nv

    @property
    def vtx(self):
        return self._vtx

    @vtx.setter
    def vtx(self, v):
        v = np.asarray(v)
        if v.ndim != 2 or np.shape(v)[1] != 3:
            raise ValueError("Vertices have wrong shape!")

        self._vtx = v

    @property
    def fac(self):
        return self._fac

    @fac.setter
    def fac(self, f):
        f = np.asarray(f)
        if f.ndim != 2 or np.shape(f)[1] != 3:
            raise ValueError("Vertices have wrong shape!")

        if np.max(f) != len(self.vtx) - 1:
            raise ValueError("Faces do not match vertex array!")

        self._fac = f


class SphereMesh(Mesh):
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

    def __init__(self, scale=1):

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
        n = self.vertex_normals
        vtx_shift = linear_interpolation3d(self.vtx[:, 0],
                                           self.vtx[:, 1],
                                           self.vtx[:, 2],
                                           self.arr_f)

        self.vtx = self.vtx + amplitude * n * vtx_shift[:, np.newaxis]

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
