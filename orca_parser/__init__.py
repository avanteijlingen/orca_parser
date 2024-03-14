# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 17:42:17 2022

@author: avtei
"""
import numpy as np
import os, glob, sys, pandas, ase



# import submodules
from .version import __version__
from orca_parser.ORCAParse import *
from orca_parser.HessianTools import *
from orca_parser.parse_engrad import *
from orca_parser.ams_parse import *

# 1 Bohr = 0.52917724900001 Angstrom

a="""

                                           _ __                                          
 ██████╗ ██████╗  ██████╗ █████╗         ██████╗  █████╗ ██████╗ ███████╗███████╗██████╗ 
██╔═══██╗██╔══██╗██╔════╝██╔══██╗        ██╔══██╗██╔══██╗██╔══██╗██╔════╝██╔════╝██╔══██╗
██║   ██║██████╔╝██║     ███████║        ██████╔╝███████║██████╔╝███████╗█████╗  ██████╔╝
██║   ██║██╔══██╗██║     ██╔══██║        ██╔═══╝ ██╔══██║██╔══██╗╚════██║██╔══╝  ██╔══██╗
╚██████╔╝██║  ██║╚██████╗██║  ██║███████╗██║     ██║  ██║██║  ██║███████║███████╗██║  ██║
 ╚═════╝ ╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝╚══════╝╚═╝     ╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝╚══════╝╚═╝  ╚═╝
                                                                                         
"""