# -*- coding: utf-8 -*-

# python standard library inputs
import os
import sys

# external inputs
import numpy as np
from numpy.linalg import norm
from nibabel.freesurfer.io import read_geometry, write_geometry
from nibabel.affines import apply_affine

# local inputs
from gbb.neighbor.nn_2d import nn_2d
from gbb.utils.update_mesh import update_mesh
from gbb.utils.smooth_surface import smooth_surface


def run_devein(vtx, fac, vtx_norm, arr_vein, arr_ignore, adjm, ras2vox_tkr, 
               n_neighbor=20, line_dir=2, smooth_iter=30, max_iterations=1000):
    """
    This function finds vertex points which are located within marked veins and 
    shift these and their neighborhood until the mesh is free of trapped 
    vertices (or the maximum number of iterations is reached). Shifts are 
    computed by averaging vertex coordinates within a neighborhood around a 
    single vertex along one direction or along all axes. All shifts are applied 
    in inward direction. Optionally, the output mesh can be smoothed.    
    Inputs.
        *vtx (arr): array of vertex points.
        *fac (arr): array of corresponding faces.
        *vtx_norm (arr): array of corresponding vertex normals.
        *arr_vein (arr): vein mask.
        *arr_ignore (arr): binary mask where deveining is omitted.
        *adjm (obj): adjacecy matrix.
        *ras2vox_tkr (arr): transformation matrix from ras to voxel space.
        *n_neighbor (int): neighborhood size.
        *line_dir (int): shift direction in ras conventions (0,1,2,3).
        *smooth_iter (int): number of smoothing iterations of final mesh.
        *max_iterations (int): maximum number of deveining iterations.
    Outputs:
        *vtx (arr): shifted array of vertex points.
        *counter (int): number of deveining iterations.
    
    created by Daniel Haenelt
    Date created: 06-02-2020             
    Last modified: 08-10-2020  
    """
    
    # fix parameters
    line_threshold = 0.05 # if direction is along one axis, omit line if length is below threshold

    # check line direction
    if line_dir > 3 or line_dir < 0:
        sys.exit("error: choose a valid line direction!")

    # check vein array
    if arr_vein is None:
        sys.exit("error: no vein mask found for deveining!")

    # load arrays
    arr_vein = np.round(arr_vein).astype(int)
    
    # get image dimensions
    xdim = np.shape(arr_vein)[0]
    ydim = np.shape(arr_vein)[1]
    zdim = np.shape(arr_vein)[2]
    
    if arr_ignore is not None:
        arr_ignore = np.round(arr_ignore).astype(int)
        arr_vein[arr_ignore == 1] = 0

    # centroid
    vtx_c = np.mean(vtx, axis=0)    

    # get nearest voxel coordinates
    vtx_vox = apply_affine(ras2vox_tkr, vtx)
    vtx_vox = np.round(vtx_vox).astype(int)
    
    # mask outlier points (ignored in deveining)
    vtx_mask = np.ones(len(vtx_vox))
    
    vtx_mask[vtx_vox[:,0] > xdim - 1] = 0
    vtx_mask[vtx_vox[:,1] > ydim - 1] = 0
    vtx_mask[vtx_vox[:,2] > zdim - 1] = 0
    
    vtx_mask[vtx_vox[:,0] < 0] = 0
    vtx_mask[vtx_vox[:,1] < 0] = 0
    vtx_mask[vtx_vox[:,2] < 0] = 0

    # get vertices trapped in veins
    vein_mask = arr_vein[vtx_vox[:,0], vtx_vox[:,1], vtx_vox[:,2]] * vtx_mask
    n_veins = len(vein_mask[vein_mask == 1])

    print("start mesh initialization (deveining)")

    # loop through vertices
    counter = 0
    while n_veins > 0 and counter < max_iterations:
    
        # print current status
        print("i: "+str(counter)+", # of trapped vertices: "+str(len(vein_mask[vein_mask == 1])))
    
        # select random vein vertex
        vein_ind = np.where(vein_mask == 1)[0]
        curr_ind = vein_ind[np.random.randint(len(vein_ind))]
        nn_ind = nn_2d(curr_ind, adjm, n_neighbor)
    
        # get shift as weighted average
        if line_dir == 3:
            vtx_shift = vtx[nn_ind,:].copy()
            vtx_shift = vtx[curr_ind,:] - vtx_shift
        else:
            vtx_shift = np.zeros((len(nn_ind),3))
            vtx_shift[:,line_dir] = vtx[nn_ind,line_dir]
            vtx_shift[:,line_dir] = vtx[curr_ind,line_dir] - vtx_shift[:,line_dir]
        
        vtx_shift = np.mean(vtx_shift, axis=0)
        
        # check inward shift by comparing distance to centroid
        vtx_dist_noshift = norm(vtx[curr_ind,:] - vtx_c)
        vtx_dist_shift = norm(vtx[curr_ind,:] - vtx_shift - vtx_c)
        
        # do only inward shifts        
        if vtx_dist_shift - vtx_dist_noshift > 0:
            vtx_shift = -1 * vtx_shift
        
        # update mesh if valid inward shift
        if line_dir != 3 and np.abs(vtx_norm[curr_ind,line_dir]) < line_threshold:
            vtx_temp_vox = apply_affine(ras2vox_tkr, vtx[curr_ind])
            vtx_temp_vox = np.round(vtx_temp_vox).astype(int)
            arr_vein[vtx_temp_vox[0],vtx_temp_vox[1],vtx_temp_vox[2]] = 0
        else:
            vtx = update_mesh(vtx, vtx_shift, curr_ind, nn_ind, 1)
            vtx_vox = apply_affine(ras2vox_tkr, vtx)
            vtx_vox = np.round(vtx_vox).astype(int)    
    
        # get all vertices within vein
        vein_mask = np.zeros(len(vtx))
        vein_mask = arr_vein[vtx_vox[:,0], vtx_vox[:,1], vtx_vox[:,2]] * vtx_mask
        n_veins = len(vein_mask[vein_mask == 1])
        
        counter += 1
    
    if counter < max_iterations:
        print("deveining converged!")
    
    # smooth output
    if smooth_iter:
        tmp = np.random.randint(0, 10, 5)
        tmp_string = ''.join(str(i) for i in tmp)
        surf_temp = os.path.join(os.getcwd(),"surf_"+tmp_string)
        write_geometry(surf_temp, vtx, fac)
        smooth_surface(surf_temp, surf_temp, smooth_iter)
        vtx, _ = read_geometry(surf_temp)
        os.remove(surf_temp)
        
    return vtx, counter
