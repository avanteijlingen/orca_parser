# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 20:22:03 2024

@author: Alex
"""
import numpy as np
from orca_parser.ORCAParse import *

def parse_engrad(engrad):
    engrad = tricky_readin(engrad)
    
    natoms = engrad.split("# The current total energy in Eh")[0]
    natoms = natoms.strip()
    natoms = int([x for x in natoms.split("\n") if x[0] != "#"][-1].strip())
    #print(natoms)
    
    Eh = engrad.split("# The current total energy in Eh")[1].split("# The current gradient in Eh/bohr")[0]
    Eh = Eh.strip()
    Eh = float([x for x in Eh.split("\n") if x[0] != "#"][-1].strip())
    #print(Eh)
    
    Forces = engrad.split("# The current gradient in Eh/bohr")[1].split("# The atomic numbers and current coordinates in Bohr")[0]
    Forces = Forces.strip()
    Forces = np.array([float(x) for x in Forces.split("\n") if x[0] != "#"])
    Forces = Forces.reshape(-1, 3)
    
    coords_ans = engrad.split("# The atomic numbers and current coordinates in Bohr")[1].strip()
    coords_ans = "\n".join([x for x in coords_ans.split("\n") if x[0] != "#"]).strip()
    coords = []
    atomic_numbers = []
    for line in coords_ans.split("\n"):
        line = line.strip().split()
        if len(line) != 4:
            continue
        atomic_numbers.append(int(line[0]))
        coords.append([float(x) for x in line[1:]])
    atomic_numbers = np.array(atomic_numbers)
    coords = np.array(coords)
    #print(atomic_numbers)
    #print(coords)
    return natoms, Eh, atomic_numbers, coords, Forces