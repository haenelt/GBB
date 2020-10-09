# -*- coding: utf-8 -*-

# python standard library inputs
import os
import sys
import json
import datetime


def write_json(io_file, io_params, devein_params, anchor_params, reg_params, 
               bbr_params, iter_params, gbb_params, path_output, prefix):
    """
    This function writes a json file to summarize input and output parameters. 
    Inputs.
        *io_file (dict): input and output filenames.
        *io_params (dict): input and output parameters.
        *devein_params (dict): deveining parameters.
        *anchor_params (dict): anchoring parameters.
        *reg_params (dict): registration parameters.
        *bbr_params (dist): cost function parameters.
        *iter_params (dict): deveining and anchoring performance parameters.
        *gbb_params (dict): gbb performance parameters.
        *path_output (str): path where output is saved.
        *name_output (str): prefix of output file.

    created by Daniel Haenelt
    Date created: 03-10-2020
    Last modified: 08-10-2020
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

    # json
    with open(file_out, 'w', encoding='utf-8') as f:
        json.dump(io_file, f, ensure_ascii=False, indent=4)
        json.dump(io_params, f, ensure_ascii=False, indent=4)
        json.dump(devein_params, f, ensure_ascii=False, indent=4)
        json.dump(anchor_params, f, ensure_ascii=False, indent=4)
        json.dump(reg_params, f, ensure_ascii=False, indent=4)
        json.dump(iter_params, f, ensure_ascii=False, indent=4)
        json.dump(bbr_params, f, ensure_ascii=False, indent=4)
        json.dump(gbb_params, f, ensure_ascii=False, indent=4)
