# -*- coding: utf-8 -*-

# python standard library inputs
import os
import sys
import json
import datetime


def write_json(io_file, io_params, devein_params, anchor_params, gbb_params, 
               bbr_params, res_params, path_output, prefix):
    """ Write JSON

    This function writes a json file to summarize input and output parameters.     

    Parameters
    ----------
    io_file : dict
        Input and output filenames.
    io_params : dict
        Input and output parameters.
    devein_params : dict
        Deveining parameters.
    anchor_params : dict
        Anchoring parameters.
    gbb_params : dict
        Registration parameters.
    bbr_params : dict
        Cost function parameters.
    res_params : dict
        Performance parameters.
    path_output : str
        Path where output is saved.
    prefix : str
        Prefix of output file.

    Returns
    -------
    None.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 03-10-2020
    Last modified: 10-10-2020

    """
    
    # make output folder    
    try:
        if not os.path.exists(path_output):
            os.makedirs(path_output)
    except TypeError:
        sys.exit("error: output directory not defined!")

    # get date string for moved files
    date = datetime.datetime.now().strftime("%Y%m%d%H%M")

    # filename
    file_out = os.path.join(path_output, prefix+"_info_"+date+".json")

    # sum all parameters in one dictionary
    dump_params = dict()
    dump_params["io_file"] = io_file
    dump_params["io_params"] = io_params
    dump_params["devein_params"] = devein_params
    dump_params["anchor_params"] = anchor_params
    dump_params["gbb_params"] = gbb_params
    dump_params["bbr_params"] = bbr_params
    dump_params["res_params"] = res_params
    
    # json
    with open(file_out, 'w', encoding='utf-8') as f:
        json.dump(dump_params, f, ensure_ascii=False, indent=4)
