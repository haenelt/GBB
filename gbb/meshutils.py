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
from config import MATRIX_SIZE, VOXEL_RES
from gbb.normal import get_normal
from .neighbor import nn_2d
from interpolation import linear_interpolation3d

__all__ = ['Mesh']


class Mesh:
    
    def __init__(self, vtx, fac):
        
        self.vtx = vtx
        self.fac = fac
        self.adjm = self._adjm(self)
    
    def _adjm(self):
        """ Get adjm
    
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
    
    def smooth(self):
        pass
    
    
    @property
    def vtx(self):
        return self._v
    
    @vtx.setter
    def vtx(self, v):
        v = np.asarray(v)
        if v.ndim != 2 or np.shape(v)[1] != 3:
            raise ValueError("Vertices have wrong shape!")
            
        self._vtx = v


    @property
    def fac(self):
        return self._f
    
    @fac.setter
    def fac(self, f):
        f = np.asarray(f)
        if f.ndim != 2 or np.shape(f)[1] != 3:
            raise ValueError("Vertices have wrong shape!")
        
        if np.max(f) != len(self.vtx):
            raise ValueError("Faces do not match vertex array!")
            
        self._fac = f


    # get normal
    # plot_normal
    # plot_normal_direction
    # adjm
    # remove vertex
    # smooth surface (laplacian)
    # remesh
    # save_img
    # save
