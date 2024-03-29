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
        step_coords = []
        
        for line in self.lines:
            if skip_lines > 0:
                skip_lines -= 1
                continue
            
            if 'Standard orientation:' in line:
                if step_coords:
                    if type(self.coords) is list:
                        self.coords = np.array(step_coords).reshape(1, -1, 3)
                    else:
                        #print(self.coords.shape, np.array(step_coords).reshape(1, -1, 3).shape)
                        self.coords = np.vstack((  self.coords, np.array(step_coords).reshape(1, -1, 3) ))                       
                    step_coords = []
                skip_lines = 4
                reading_coordinates = True
                continue
    
            if reading_coordinates:
                if "---------------------------------------------------------------------" in line:
                    reading_coordinates = False
                    continue
    
                parts = line.strip().split()
                if len(parts) >= 6: 

                    step_coords.append(list(map(float, parts[3:6])))