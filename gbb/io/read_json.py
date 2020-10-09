# -*- coding: utf-8 -*-

# python standard library inputs
import json


def read_json(file_in):
    """
    This function reads a json file and returns the content in a dictionary.
    Inputs.
        *file_in (str): filename of json file.
    Outputs:
        *data (dict): content of json file. 

    created by Daniel Haenelt
    Date created: 09-10-2020
    Last modified: 09-10-2020
    """
    
    with open(file_in) as json_file: 
        data = json.load(json_file)
        
    return data
