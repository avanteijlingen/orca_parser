#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 15:20:02 2024

@author: bwb16179
"""

import numpy as np
import sys
from . import ORCAParse
import numpy as np

class GaussianParse(ORCAParse):
    def ValidateOutput(self):
        self.valid = False
        if "Normal termination of Gaussian" in self.raw:
            if "SCF Done:" in self.raw:
                self.valid = True
                
    def __init__(self, filepath, verbose=False):
        super().__init__(filepath, verbose = verbose)
        self.filepath = filepath

        self.lines = self._read_lines()
        if self.validate_output():
            if self.verbose:
                print("File terminated Normally")
        else:
            if self.verbose:
                print("The Gaussian output file did not terminate normally or is not valid.")

    def _read_lines(self):
        """Read the file and store lines for further processing."""
        with open(self.filepath, 'r') as file:
            lines = file.readlines()
        return lines
    
    def parse_atoms(self):
        atms, line_skip = False, 0
        for line in self.lines:
            if 'Charge =  ' in line:
                atms = True
                continue
            if " Redundant internal coordinates found in file." in line:
                atms = False
                continue
            if atms:
                parts = line.strip().split()
                if len(parts) == 0:
                    atms = False
                    continue
                atm = parts[0].split("(PDBName")[0]
                self.atoms.append(atm)
        for atom in self.atoms:
            self.masses.append(self.Masses[atom])
            

    def validate_output(self):
        """Validate the Gaussian output file."""
        return any('Normal termination of Gaussian' in line for line in self.lines)

    def parse_energies(self):
        """Parse the SCF energies from the output file."""
        self.energies = np.array([])
        for line in self.lines:
            if 'SCF Done:' in line:
                energy = float(line.split()[4])
                self.energies = np.append(self.energies, energy)
                continue

    def parse_coords(self):
        """Parse the coordinates from each step of the optimization."""
        reading_coordinates, skip_lines = False, 0
        self.coords = np.ndarray((0,3))
        
        coord_blocks = self.raw.split('Standard orientation:')[1:]
        for block in coord_blocks:
            natoms = 0
            block = block.split("---------------------------------------------------------------------")[2]
            for line in block.split("\n"):
                line = line.split()
                if len(line) != 6:
                    continue
                self.coords = np.vstack((  self.coords, [float(x) for x in line[3:]] ))                      
                natoms += 1
        self.coords = self.coords.reshape((-1, natoms, 3))
        # Sometimes the last one is doubled up in the file
        if self.coords.shape[0] > 1:
            if (self.coords[-1] == self.coords[-2]).all():
                self.coords = np.delete(self.coords, [-1], axis=0)
        