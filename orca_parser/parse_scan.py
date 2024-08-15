#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 16:55:00 2024

@author: bwb16179
"""


from orca_parser.ORCAParse import *
import numpy as np

class parse_scan(ORCAParse):
    def __init__(self, file):
        self.op = ORCAParse(file)
        self.op.ValidateOutput()
        
        self.scan_atoms = []
        
    def parse_scan_coords(self):
        content = self.op.raw
        self.num_steps = int(content.split('%geom Scan', 1)[1].splitlines()[1].split()[-1])
        
        frames =  content.split("*** FINAL ENERGY EVALUATION AT THE STATIONARY POINT ***")[1:]
        
        # assert self.num_steps == frames, "The number of scan steps should equal the number of optimized coordinates. This job has not finished correctly"
        
        for i,frame in enumerate(frames):
            positions = []
            
            part = frame.split("CARTESIAN COORDINATES (ANGSTROEM)")[1:]
            part = frame.split("CARTESIAN COORDINATES (A.U.)")[0]
            for line in part.split("\n"):
                line = line.split()
                if len(line) != 4:
                    continue
                positions.append([float(x) for x in line[1:]])
                if i == 0:
                    self.scan_atoms.append(line[0])
            if i == 0:
                self.scan_coords = np.array(positions).reshape(1, -1, 3)
            else:
                self.scan_coords = np.vstack((self.scan_coords, np.array(positions).reshape(1, -1, 3)))
    