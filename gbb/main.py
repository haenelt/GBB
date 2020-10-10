# -*- coding: utf-8 -*-

# python standard library inputs
import os
import sys
import importlib.util

# local inputs
from gbb.io.get_filename import get_filename
from gbb.config import config
from gbb.io.load_data import load_data
from gbb.io.write_json import write_json
from gbb.process.run_devein import run_devein
from gbb.process.run_anchor import run_anchor
from gbb.process.run_gbb import run_gbb
from gbb.utils.deformation_field import deformation_field
from gbb.utils.apply_deformation import apply_deformation


def main(file_white, file_ref, path_output, file_pial=None, file_vein=None, 
         file_ignore=None, file_anchor=None, file_config=None):
    """
    This function executes the gradient-based boundary (GBB) surface refinement.
    Inputs:
        *file_white (str): filename of gm/wm boundary surface mesh.
        *file_ref (str): filename of reference volume.
        *path_output (str): path where output is written.
        *file_pial (str): filename of gm/csf boundary surface mesh (optional).
        *file_vein (str): filename of vein mask (optional).
        *file_ignore (str): filename of ignore mask (optional).
        *file_anchor (str): filename of anchor points (optional).
        *file_config (str): filename of a custom configuration file (optional).
    
    created by Daniel Haenelt
    Date created: 03-10-2020
    Last modified: 10-10-2020
    """

    # input and output parameters
    io_file = dict()
    io_file["i_white"] = file_white
    io_file["i_pial"] = file_pial
    io_file["i_ref"] = file_ref
    io_file["i_vein"] = file_vein
    io_file["i_ignore"] = file_ignore
    io_file["i_anchor"] = file_anchor
    io_file["o_output"] = path_output

    # initialize dictionary for performance parameters
    res_params = dict()
    res_params["devein"] = dict()
    res_params["anchor"] = dict()
    res_params["gbb"] = dict()

    # load configurations
    io_params = config.io_params
    devein_params = config.devein_params
    anchor_params = config.anchor_params
    gbb_params = config.gbb_params
    bbr_params = config.bbr_params

    # load custom configurations
    if file_config:
        path_config, name_config, _ = get_filename(file_config)
        spec = importlib.util.spec_from_file_location(name_config, file_config)
        config_custom = importlib.util.module_from_spec(spec)
        sys.modules[name_config] = config_custom
        spec.loader.exec_module(config_custom)
        
        io_params = config_custom.io_params
        devein_params = config_custom.devein_params
        anchor_params = config_custom.anchor_params
        gbb_params = config_custom.gbb_params
        bbr_params = config_custom.bbr_params

    # make output folder    
    try:
        if not os.path.exists(io_file["o_output"]):
            os.makedirs(io_file["o_output"])
    except TypeError:
        sys.exit("error: output directory not defined!")
    
    # load input
    volume, T, surf, point, basename = load_data(io_file, gbb_params)
    
    # run deveining
    if devein_params["run"]:        
        surf["vtx_white"], res_params["devein"]["niter"] = run_devein(surf["vtx_white"], 
                                                                      surf["fac_white"], 
                                                                      surf["n_white"],
                                                                      volume["vein"], 
                                                                      volume["ignore"],  
                                                                      T["adjm"], 
                                                                      T["ras2vox"],
                                                                      devein_params["n_neighbor"], 
                                                                      3, # along all directions
                                                                      devein_params["n_smooth"], 
                                                                      devein_params["max_iter"])   
    
    # run anchoring
    ind_control = []
    if anchor_params["run"]:
        surf["vtx_white"], ind_control = run_anchor(surf["vtx_white"], 
                                                    surf["fac_white"], 
                                                    T["adjm"], 
                                                    point["anchor"],
                                                    anchor_params["n_neighbor"], 
                                                    anchor_params["n_smooth"])
        
        res_params["anchor"]["niter"] = len(ind_control)
    
    # run gbb
    if gbb_params["run"]:
        
        # do not consider control points in gbb
        if gbb_params["overwrite_control"]:
            ind_control = []
        
        surf["vtx_white"], res_params["gbb"] = run_gbb(surf["vtx_white"],
                                                       surf["fac_white"], 
                                                       surf["n_white"],
                                                       ind_control, 
                                                       volume["ref"], 
                                                       volume["gradient"], 
                                                       volume["vein"], 
                                                       volume["ignore"], 
                                                       bbr_params["t2s"], 
                                                       T["vox2ras"], 
                                                       T["ras2vox"], 
                                                       gbb_params["line_dir"], 
                                                       gbb_params["line_length"], 
                                                       gbb_params["r_size"], 
                                                       gbb_params["l_rate"], 
                                                       gbb_params["max_iter"], 
                                                       gbb_params["cost_threshold"], 
                                                       gbb_params["cost_step"], 
                                                       gbb_params["cost_sample"], 
                                                       io_file["o_output"], 
                                                       gbb_params["show_cost"], 
                                                       gbb_params["show_slope"], 
                                                       gbb_params["intermediate_write"],
                                                       bbr_params["Q0"],
                                                       bbr_params["M"],
                                                       bbr_params["h"],
                                                       bbr_params["s"])
    
    # write deformation field
    print("write deformation field")
    arr_deformation = deformation_field(surf["vtx_white_archive"],
                                        surf["vtx_white"],
                                        io_file["i_ref"],
                                        T["ras2vox"],
                                        io_params["o_sigma"],
                                        io_file["o_output"],
                                        basename["white"]+"_deformation",
                                        True)
          
    # apply deformation to white and pial surface
    print("apply deformation")
    _, _ = apply_deformation(surf["vtx_white_archive"],
                             surf["fac_white"],
                             arr_deformation,
                             T["vox2ras"],
                             T["ras2vox"],
                             io_file["o_output"],
                             basename["white"]+"_refined", 
                             True)
    
    if surf["vtx_pial"] is not None:
        _, _ = apply_deformation(surf["vtx_pial"], 
                                 surf["fac_pial"], 
                                 arr_deformation,
                                 T["vox2ras"],
                                 T["ras2vox"],
                                 io_file["o_output"], 
                                 basename["pial"]+"_refined", 
                                 True)
    
    # write json   
    if file_config:
        io_params["config_source"] = config_custom.__file__
    else:
        io_params["config_source"] = config.__file__
   
    write_json(io_file,
               io_params,
               devein_params,
               anchor_params,
               gbb_params,
               bbr_params,
               res_params,
               io_file["o_output"],
               basename["white"])
