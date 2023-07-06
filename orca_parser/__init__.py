# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 17:42:17 2022

@author: avtei
"""
import numpy as np
import os, glob, sys, pandas, ase


# import submodules
from orca_parser.ORCAParse import *
from orca_parser.HessianTools import *

# 1 Bohr = 0.52917724900001 Angstrom




        
        
if __name__ == "__main__":
    Hess = HessianTools("Test-cases/Coo/Coo.hess")
    print(Hess.normalmodes)

        
