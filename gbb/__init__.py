# -*- coding: utf-8 -*-

# python standard library inputs
import os

# local inputs
from .config import *
from .interpolation import *
from .io import *
from .neighbor import *
from .normal import *
from .plot import *
from .process import *
from .utils import *


"""
Python package for gradient-based boundary (GBB) surface refinement.

created by Daniel Haenelt
Date created: 05-10-2020
Last modified: 05-10-2020
"""

path_root = os.path.dirname(os.path.dirname(__file__))
with open(os.path.join(path_root, 'VERSION')) as version_file:
    version = version_file.read().strip()

__author__ = "Daniel Haenelt"
__license__ = "GPL v3"
__version__ = version
__status__ = "Development"
