# -*- coding: utf-8 -*-

# python standard library inputs
import json


def read_json(file_in):
    """ Read JSON

    This function reads a json file and returns the content in a dictionary.    

    Parameters
    ----------
    file_in : str
        Filename of json file.

    Returns
    -------
    data : dict
        Content of json file. 

    Notes
    -------
    created by Daniel Haenelt
    Date created: 09-10-2020
    Last modified: 09-10-2020

    """
    
    with open(file_in) as json_file: 
        data = json.load(json_file)
        
    return data
