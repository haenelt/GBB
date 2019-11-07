"""
White2Pial lines

The purpose of the following script is to generate lines between corresponding vertices at the white
and pial surface to visualize the shift between matched vertices caused by realigning surfaces 
independently. You can either construct prisms, triangles or lines.

created by Daniel Haenelt
Date created: 07-11-2019
Last modified: 07-11-2019
"""
import os
import numpy as np
from nibabel.freesurfer.io import read_geometry, write_geometry

input_white = "/home/daniel/Schreibtisch/p4/lh.white"
input_pial = "/home/daniel/Schreibtisch/p4/lh.pial"
path_output = "/home/daniel/Schreibtisch/"
hemi = "lh"
basename_output = "bla"
shape = "line" # line, triangle or prism
step_size = 100

""" do not edit below """

# read geometry
vtx_white, fac_white, header_white = read_geometry(input_white, read_metadata=True)
vtx_pial, fac_pial = read_geometry(input_pial)

# array containing a list of considered vertices
t = np.arange(0,len(vtx_white),step_size)

# initialise faces for specific shape
if shape is "prism":
    fac_new = [[0, 1, 2],
               [3, 4, 5],
               [0, 1, 4],
               [0, 3, 4],
               [1, 2, 5],
               [1, 4, 5],
               [0, 2, 5],
               [0, 3, 5]]
    fac_iter = 6
elif shape is "triangle":
    fac_new = [[0,1,2]]
    fac_iter = 3
elif shape is "line":
    fac_new = [[0,1,0]]
    fac_iter = 2

vtx_res = []
fac_res = []
for i in range(len(t)):
    
    # get triangle from nearest neighbour point of a given vertex
    neighbour, _ = np.where(fac_white == t[i])
    neighbour = fac_white[neighbour,:]
    neighbour = np.concatenate(neighbour)
    neighbour = neighbour[neighbour != t[i]]
    neighbour = neighbour[:2]
    
    # get all vertex points for specific shape
    if shape is "prism":
        A = list(vtx_white[t[i]])
        B = list(vtx_white[neighbour[0]])
        C = list(vtx_white[neighbour[1]])
        D = list(vtx_pial[t[i]])
        E = list(vtx_pial[neighbour[0]])
        F = list(vtx_pial[neighbour[1]])
        vtx_new = [A,B,C,D,E,F]
    elif shape is "triangle":
        A = list(vtx_white[t[i]])
        B = list(vtx_white[neighbour[0]])
        C = list(vtx_pial[t[i]])
        vtx_new = [A,B,C]
    elif shape is "line":
        A = list(vtx_white[t[i]])
        B = list(vtx_pial[t[i]])
        vtx_new = [A,B]

    # update faces
    if i > 0:
        for j in range(len(fac_new)):
            fac_new[j] = [x+fac_iter for x in fac_new[j]]
    
    # update resulting vertex and face list
    vtx_res.extend(vtx_new)
    fac_res.extend(fac_new)

# vertices and faces as array
vtx_res = np.array(vtx_res)
fac_res = np.array(fac_res)

# write output geometry
write_geometry(os.path.join(path_output,hemi+"."+basename_output+".line"),
               vtx_res, fac_res, volume_info=header_white)
