# -*- coding: utf-8 -*-
import argparse
import gbb
from gbb.run_gbb import run_gbb


"""
GBB

Run the GBB package using the command line via python gbb <args>.

created by Daniel Haenelt
Date created: 03-10-2020
Last modified: 03-10-2020
"""

# parse arguments from command line
parser_description = "Gradient-based boundary (GBB) surface refinement."
parser = argparse.ArgumentParser(description=parser_description)
parser.add_argument('-w','--white', type=str, help='input white surface')
parser.add_argument('-r','--ref', type=str, help='input reference volume')
parser.add_argument('-o','--output', type=str, help='output path')
parser.add_argument('-p','--pial', type=str, help='input pial surface', default=None)
parser.add_argument('-v','--vein', type=str, help='input vein mask', default=None)
parser.add_argument('-i','--ignore', type=str, help='input mask', default=None)
parser.add_argument('-a','--anchor', type=str, help='input anchor points', default=None)
args = parser.parse_args()

# check if mandatory arguments are defined
if not args.white:
    parser.error("-w is not defined!")

if not args.ref:
    parser.error("-r is not defined!")

if not args.output:
    parser.error("-o is not defined!")

# run main module
print("-----------------------------------------------------------------------")
print("GBB "+"(v"+str(gbb.__version__)+")")
print("author: "+str(gbb.__author__))
print("-----------------------------------------------------------------------")

run_gbb(args.white, 
        args.ref, 
        args.output, 
        args.pial, 
        args.vein, 
        args.ignore, 
        args.anchor)
