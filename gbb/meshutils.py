# -*- coding: utf-8 -*-
"""Generation of test objects (sphere, cube) and download of test data."""

# python standard library inputs
import os

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
from .config import MATRIX_SIZE, VOXEL_RES
from gbb.normal import get_normal
from .neighbor import nn_2d
from .interpolation import linear_interpolation3d

__all__ = ['Mesh']


class Mesh:
    
    def __init__(self, vtx, fac):
        
        self.vtx = vtx
        self.fac = fac
        self.adjm = self._get_adjm()
    
    def _get_adjm(self):
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
    
    def smooth(self, niter):
        vtx_copy = self.vtx.copy()
        adjm_ind = self.adjm.indices
        adjm_ptr = self.adjm.indptr
        for i in range(niter):
            print(i)
            vtx_copy = [np.mean(vtx_copy[adjm_ind[adjm_ptr[j]:adjm_ptr[j+1]], :], 0) for j in range(len(vtx_copy))]
            vtx_copy = np.asarray(vtx_copy)

        return vtx_copy
    


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


    # get normal
    # plot_normal
    # plot_normal_direction
    # adjm
    # smooth surface (laplacian)
    # remove vertex
    # remesh
    # save_img
    # save
