# -*- coding: utf-8 -*-
import os
import sys
import nibabel as nb
from nibabel.freesurfer.io import read_geometry
from gbb.io.read_anchor import read_anchor
from gbb.normal.get_normal import get_normal
from gbb.utils.get_adjm import get_adjm
from gbb.utils.get_gradient import get_gradient
from gbb.utils.vox2ras import vox2ras


def load_data(io_file, reg_params):
    """
    This function loads volume arrays, transformation matrices between voxel and 
    ras space, surface data and point sets.
    Inputs:
        *io_file (dict): input and output parameters.
        *reg_params (dict): registration parameters.
    Outputs:
        *volume (dict): loaded volume arrays.
        *T (dict): loaded transformation matrices.
        *surf (dict): loaded surface vertices and faces.
        *point (dict): loaded point set.
        *basename (dict): basenames.
        *path_temp (str): path of temporary folder.
        
    created by Daniel Haenelt
    Date created: 12-05-2020
    Last modified: 03-10-2020
    """
    
    # initialize dictionaries
    surf = dict()
    basename = dict()
    volume = dict()
    T = dict()
    point = dict()
    
    # make output folder
    if not os.path.exists(io_file["o_output"]):
        os.makedirs(io_file["o_output"])
    
    # load white surface
    if io_file["i_white"]:
        surf["vtx_white"], surf["fac_white"] = read_geometry(io_file["i_white"])
        surf["vtx_white_archive"] = surf["vtx_white"].copy()
        
        # get normals
        surf["n_white"] = get_normal(surf["vtx_white"], surf["fac_white"])
        
        # adjacency matrix
        T["adjm"] = get_adjm(io_file["i_white"])
        
        # basename
        basename["white"] = os.path.basename(io_file["i_white"])
    else:
        sys.exit("White surface is not found!")
    
    # load pial surface
    if io_file["i_pial"]:
        surf["vtx_pial"], surf["fac_pial"] = read_geometry(io_file["i_pial"])
        
        # basename
        basename["pial"] = os.path.basename(io_file["i_pial"])
    else:
        print("No pial surface is loaded")
        surf["vtx_pial"] = None
        surf["fac_pial"] = None
        basename["pial"] = None
    
    # load reference volume
    if io_file["i_ref"]:
        volume["instance"] = nb.load(io_file["i_ref"])
        volume["ref"] = volume["instance"].get_fdata()
    
        # get transformation
        T["vox2ras"], T["ras2vox"] = vox2ras(io_file["i_ref"])
    
        # get gradient
        volume["gradient"] = get_gradient(io_file["i_ref"], 
                                          T["ras2vox"], 
                                          reg_params["line_dir"], 
                                          reg_params["gradient_sigma"], 
                                          reg_params["gradient_kernel"], 
                                          reg_params["gradient_write"], 
                                          io_file["o_output"])
    else:
        sys.exit("Reference volume is not found!")
    
    # load vein mask
    if io_file["i_vein"]:
        volume["vein"] = nb.load(io_file["i_vein"]).get_fdata()
    else:
        print("No vein mask is loaded")
        volume["vein"] = None
    
    # load ignore volume
    if io_file["i_ignore"]:
        volume["ignore"] = nb.load(io_file["i_ignore"]).get_fdata()
    else:
        print("No ignore mask is loaded")
        volume["ignore"] = None
    
    # load anchor points
    if io_file["i_anchor"]:
        point["anchor"] = read_anchor(io_file["i_anchor"])
    else:
        print("No anchor points are loaded")
        point["anchor"] = None
    
    return volume, T, surf, point, basename
