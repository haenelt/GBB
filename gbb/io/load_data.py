# -*- coding: utf-8 -*-

# python standard library inputs
import os
import sys

# external inputs
import nibabel as nb
from nibabel.freesurfer.io import read_geometry

# local inputs
from gbb.io.read_anchor import read_anchor
from gbb.normal.get_normal import get_normal
from gbb.utils.get_adjm import get_adjm
from gbb.utils.get_gradient import get_gradient
from gbb.utils.vox2ras import vox2ras


def load_data(io_file, gbb_params):
    """
    This function loads volume arrays, transformation matrices between voxel and 
    ras space, surface data and point sets.
    Inputs:
        *io_file (dict): input and output parameters.
        *gbb_params (dict): registration parameters.
    Outputs:
        *volume (dict): loaded volume arrays.
        *T (dict): loaded transformation matrices.
        *surf (dict): loaded surface vertices and faces.
        *point (dict): loaded point set.
        *basename (dict): basenames.
        
    created by Daniel Haenelt
    Date created: 12-05-2020
    Last modified: 11-10-2020
    """
    
    # initialize dictionaries
    surf = dict()
    basename = dict()
    volume = dict()
    T = dict()
    point = dict()
    
    # make output folder    
    try:
        if not os.path.exists(io_file["o_output"]):
            os.makedirs(io_file["o_output"])
    except TypeError:
        sys.exit("error: output directory not defined!")
    
    # load white surface
    if io_file["i_white"]:
        try:
            surf["vtx_white"], surf["fac_white"] = read_geometry(io_file["i_white"])
        except FileNotFoundError:
            sys.exit("error: white surface not found!")

        surf["vtx_white_archive"] = surf["vtx_white"].copy()
        
        # get normals
        surf["n_white"] = get_normal(surf["vtx_white"], surf["fac_white"])
        
        # adjacency matrix
        T["adjm"] = get_adjm(surf["vtx_white"], surf["fac_white"])
        
        # basename
        basename["white"] = os.path.basename(io_file["i_white"])
    else:
        sys.exit("error: white surface not defined!")

    # load reference volume
    if io_file["i_ref"]:
        try:
            volume["instance"] = nb.load(io_file["i_ref"])
        except FileNotFoundError:
            sys.exit("error: reference volume not found!")

        volume["ref"] = volume["instance"].get_fdata()
    
        # get transformation
        T["vox2ras"], T["ras2vox"] = vox2ras(io_file["i_ref"])
    
        # get gradient
        volume["gradient"] = get_gradient(io_file["i_ref"], 
                                          T["ras2vox"], 
                                          gbb_params["line_dir"], 
                                          gbb_params["gradient_sigma"], 
                                          gbb_params["gradient_kernel"], 
                                          gbb_params["gradient_write"], 
                                          io_file["o_output"],
                                          basename["white"]+"_gradient")
    else:
        sys.exit("error: reference volume not defined!")
    
    # load pial surface
    if io_file["i_pial"]:
        try:
            surf["vtx_pial"], surf["fac_pial"] = read_geometry(io_file["i_pial"])
        except FileNotFoundError:
            sys.exit("error: pial surface not found!")
        
        # basename
        basename["pial"] = os.path.basename(io_file["i_pial"])
    else:
        print("no pial surface loaded")
        surf["vtx_pial"] = None
        surf["fac_pial"] = None
        basename["pial"] = None
        
    # load vein mask
    if io_file["i_vein"]:
        try:
            volume["vein"] = nb.load(io_file["i_vein"]).get_fdata()
        except FileNotFoundError:
            sys.exit("error: vein mask not found!")
    else:
        print("no vein mask loaded")
        volume["vein"] = None
    
    # load ignore volume
    if io_file["i_ignore"]:
        try:
            volume["ignore"] = nb.load(io_file["i_ignore"]).get_fdata()
        except FileNotFoundError:
            sys.exit("error: ignore mask not found!")
    else:
        print("no ignore mask loaded")
        volume["ignore"] = None
    
    # load anchor points
    if io_file["i_anchor"]:
        try:
            point["anchor"] = read_anchor(io_file["i_anchor"])
        except FileNotFoundError:
            sys.exit("error: anchor points not found!")
    else:
        print("no anchor points loaded")
        point["anchor"] = None
    
    return volume, T, surf, point, basename
