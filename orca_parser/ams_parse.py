# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 17:11:35 2024

@author: Alex
"""
from . import ORCAParse
import numpy as np

#### Module for reading Amsterdam Modeling Suite outputs
class ams_parse(ORCAParse):
    def ValidateOutput(self):
        self.valid = False
        if "NORMAL TERMINATION" in self.raw:
            if "Energy (hartree)" in self.raw:
                self.valid = True
    def parse_energies(self):
        self.energies = np.ndarray((0,), np.float64)
        for part in self.raw.split("Energy (hartree)   ")[1:]:
            part = part.split("\n")[0].strip()
            part = float(part)
            self.energies = np.hstack((self.energies, [part]))
            
    def parse_coords(self):
        self.coords = []
        self.atoms = []
        self.masses = []
        frames = self.raw.split("Geometry")[1:]
        for i,frame in enumerate(frames):
            positions = []
            frame = frame.split("\n\n")[0]
            for line in frame.split("\n"):
                line = line.split()
                if len(line) != 5:
                    continue
                positions.append([float(x) for x in line[2:]])
                if i == 0:
                    self.atoms.append(line[1])
            if i == 0:
                self.coords = np.array(positions).reshape(1, -1, 3)
            else:
                self.coords = np.vstack((self.coords, np.array(positions).reshape(1,-1,3)))
        for x in self.atoms:
            self.masses.append(self.Masses[x])
        
    def __init__(self, fname, verbose = False):
        super().__init__(fname, verbose = False)
        self.ValidateOutput()