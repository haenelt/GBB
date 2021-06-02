# -*- coding: utf-8 -*-
"""Generation of test objects (sphere, cube) and download of test data."""

# python standard library inputs
import os
import itertools

# external inputs
import numpy as np
import pyvista as pv
import nibabel as nb
import gdown
from scipy.sparse import csr_matrix
from numpy.fft import fftn, ifftn, fftshift
from numpy.random import normal
from nibabel.freesurfer.io import write_geometry

# local inputs
#from .config import MATRIX_SIZE, VOXEL_RES
#from gbb.normal import get_normal
#from .neighbor import nn_2d
#from .interpolation import linear_interpolation3d

__all__ = ['Mesh']


class Mesh:
    
    def __init__(self, vtx, fac):
        
        self.vtx = vtx
        self.fac = fac



    @property
    def adjm(self):
        """Get adjm
    
        This function computes a sparse adjacency matrix for a triangular surface
        mesh. The matrix has the size (nvertex,nvertex). Each matrix entry with 
        value 1 stands for an edge of the surface mesh.    
    
        Parameters
        ----------
        vtx : ndarray
            Array of vertices.
        fac : ndarray
            Array of faces.
    
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
        row.extend( [self.fac[i,0] for i in range(nfac)] )
        col.extend( [self.fac[i,1] for i in range(nfac)] )
        
        row.extend( [self.fac[i,1] for i in range(nfac)] )
        col.extend( [self.fac[i,2] for i in range(nfac)] )
        
        row.extend( [self.fac[i,2] for i in range(nfac)] )
        col.extend( [self.fac[i,0] for i in range(nfac)] )
            
        # make sure that all edges are symmetric
        row.extend( [self.fac[i,1] for i in range(nfac)] )
        col.extend( [self.fac[i,0] for i in range(nfac)] )
        
        row.extend( [self.fac[i,2] for i in range(nfac)] )
        col.extend( [self.fac[i,1] for i in range(nfac)] )
        
        row.extend( [self.fac[i,0] for i in range(nfac)] )
        col.extend( [self.fac[i,2] for i in range(nfac)] )
        
        # adjacency entries get value 1
        data = np.ones(len(row), dtype=np.int8)
        
        # write sparse adjacency matrix
        sparse_adjm = csr_matrix((data, (row, col)), shape=(nvtx,nvtx))
    
        return sparse_adjm

    @property
    def vertex_normals(self):
        """ Get normal

        This function computes the surfaces normals per vertex from an input surface
        mesh. The code is taken from [1] and adapted to my own purposes.

        Parameters
        ----------
        vtx : ndarray
            Vertex array.
        fac : ndarray
            Face array.

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

        # initialize array with same type and shape as vtx
        norm = np.zeros(self.vtx.shape, dtype=self.vtx.dtype)

        # indexed view into the vertex array
        tris = self.vtx[self.fac]

        # calculate the normal for all triangles by taking the cross product of
        # the vectors v1-v0 and v2-v0 in each triangle
        n = np.cross(tris[::, 1] - tris[::, 0], tris[::, 2] - tris[::, 0])

        # calculate vertex-wise normals and normalize
        self.normal = self._f2v(len(self.vtx), self.fac, n)
        self.normal = self._normalize_v3(self.normal)

        return self.normal

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

    def smooth(self, niter):
        vtx_copy = self.vtx.copy()
        adjm_ind = self.adjm.indices
        adjm_ptr = self.adjm.indptr
        for i in range(niter):
            print(i)
            vtx_copy = [np.mean(vtx_copy[adjm_ind[adjm_ptr[j]:adjm_ptr[j+1]], :], 0) for j in range(len(vtx_copy))]
            vtx_copy = np.asarray(vtx_copy)

        return vtx_copy

    def save_normals(self, file_out):
        """ Plot normal

        This function generates lines to visualize outward directed surface normals
        of an input surface mesh.

        Parameters
        ----------
        vtx : ndarray
            Array of vertex points.
        fac : ndarray
            Corresponding face array.
        adjm : obj
            Adjacency matrix.
        file_out : str
            Filename of output surface mesh.
        step_size : int, optional
            Subset of vertices. The default is 100.
        shape : str, optional
            line, triangle, prism. The default is "line".

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
            A = list(self.vtx[i,:])
            B = list(self.vtx[i,:] - normals[i,:])
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
        input_surf : str
            Filename of source mesh (white surface).
        axis : int, optional
            Axis for distance calculation in ras space (0,1,2). The default is 2.

        Returns
        -------
        None.

        Notes
        -------
        created by Daniel Haenelt
        Date created: 13-12-2019
        Last modified: 05-10-2020

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
        """ Normalize a numpy array of 3 component vectors shape=(n,3) """
        lens = np.sqrt(arr[:, 0] ** 2 + arr[:, 1] ** 2 + arr[:, 2] ** 2)
        res = np.zeros_like(arr)
        res[:, 0] = arr[:, 0] / lens
        res[:, 1] = arr[:, 1] / lens
        res[:, 2] = arr[:, 2] / lens
        return res

    @staticmethod
    def _f2v(length, fac, nt):
        """ get average vertex-wise normal by adding up all face-wise normals
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


    # save_normals
    # save_normal_dir
    # save_img
    # adjm
    # vertex_normals
    # smooth surface (laplacian)
    # save
    # remove vertex
    # remesh


#from nibabel.freesurfer.io import read_geometry
#vtx, fac = read_geometry("/home/daniel/Schreibtisch/test/data/sphere")
#mesh = Mesh(vtx, fac)