# -*- coding: utf-8 -*-

# python standard library inputs
import sys
import subprocess
import argparse

# local inputs
import gbb
from gbb.main import main


"""
GBB

Run the GBB package using the command line via python -m gbb <args>.

created by Daniel Haenelt
Date created: 03-10-2020
Last modified: 05-10-2020
"""

# description
parser_description = "Gradient-based boundary (GBB) surface refinement. A " \
                     "GM/WM boundary surface (-w), a reference and an output " \
                     "directory (-o) are mandatory arguments. Optionally, a " \
                     "configuration file (-c) can be defined. Depending on " \
                     "the parameters in the configuration file, a vein mask " \
                     "(-v) and/or anchoring points (-a) are mandatory inputs " \
                     "as well. Additionally, a mask (-i) can be used to " \
                     "exlude regions from processing. A second surface (-p) " \
                     "can be used to apply the same deformation onto another " \
                     "surface mesh."

# parse arguments from command line
parser = argparse.ArgumentParser(description=parser_description)
parser.add_argument('-p','--pial', type=str, help='input pial surface', default=None)
parser.add_argument('-i','--ignore', type=str, help='input ignore mask', default=None)
parser.add_argument('-c','--config', type=str, help='configuration file', default=None)
parser.add_argument('-a','--anchor', type=str, help='input anchor points', default=None)
parser.add_argument('-v','--vein', type=str, help='input vein mask', default=None)
requiredNames = parser.add_argument_group('mandatory arguments')
requiredNames.add_argument('-w', '--white', help='input white surface')
requiredNames.add_argument('-r','--ref', type=str, help='input reference volume')
requiredNames.add_argument('-o','--output', type=str, help='output path')
args = parser.parse_args()

# check freesurfer installation
try:
    freesurfer_version = subprocess.check_output(['recon-all', '--{}'.format("version")]).decode()
    freesurfer_version = freesurfer_version.rstrip()
except FileNotFoundError:
    sys.exit("No freesurfer installation found!")

# run main module
print("-----------------------------------------------------------------------")
print("GBB "+"(v"+str(gbb.__version__)+")")
print("author: "+str(gbb.__author__))
print("use: "+str(freesurfer_version))
print("-----------------------------------------------------------------------")

main(args.white, 
     args.ref, 
     args.output, 
     args.pial, 
     args.vein, 
     args.ignore, 
     args.anchor,
     args.config)
