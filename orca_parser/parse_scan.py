#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 16:55:00 2024

@author: bwb16179
"""


from . import ORCAParse
import numpy as np

class parse_scan(ORCAParse):
    def __init__(self, file):
        self.op = ORCAParse(file)
        self.op.ValidateOutput()
        
    def parse_scan_coords(self):
        """
        Function to parse a surface scan output from ORCA.

        Returns
        -------
        self.scan_atoms: list, the atomic symbols from atoms in the system
        self.scan_coords: numpy.array, scanned coords from the output of every relaxed optimization calculation in the file. shape (Num_Frames,Num_atoms,3)

        """
        
        self.scan_atoms = []
        self.num_steps = int(self.op.raw.split('%geom Scan', 1)[1].splitlines()[1].split()[-1])
        
        frames =  self.op.raw.split("*** FINAL ENERGY EVALUATION AT THE STATIONARY POINT ***")[1:]
                
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
                
    def parse_scan_energies(self):
        """
        Function to return the 'Final Single Point Energy' from each sep of a relaxed surface scan in ORCA.
        
        units: Hartree

        Returns
        -------
        self.scan_energies: list, list of 'Actual' energies from final energy evaluations at stationary points

        """
        self.scan_energies = []
        energies_list = self.op.raw.split("*** OPTIMIZATION RUN DONE ***")[:-1]
        for energy in energies_list:
            E = energy.split("\n")[-4]
            self.scan_energies.append(float(E.split()[-1]))
                       