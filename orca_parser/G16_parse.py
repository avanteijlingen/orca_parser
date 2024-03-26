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
    def __init__(self, filepath, verbose):
        super().__init__(filepath, verbose = verbose)
        self.filepath = filepath
        self.atoms = []  # Atom symbols, populated once
        self.masses = []  # Atomic masses, populated once
        self.energies = np.array([])  # Energies for each step
        self.lines = self._read_lines()
        if self.validate_output():
            self.parse_energies()
            self.parse_coords()
            self.parse_atoms()
        else:
            print("The Gaussian output file did not terminate normally or is not valid.")

    def _read_lines(self):
        """Read the file and store lines for further processing."""
        with open(self.filepath, 'r') as file:
            lines = file.readlines()
        return lines
    
    def parse_atoms(self):
        parse_atoms, line_skip = False, 0
        for line in self.lines:
            if 'Charge =  ' in line:
                parse_atoms = True
                continue
            if " Redundant internal coordinates found in file.  (old form)." in line:
                parse_atoms = False
                continue
            if parse_atoms:
                parts = line.strip().split()
                if len(parts) == 0:
                    parse_atoms = False
                    continue
                print(parts)
                self.atoms.append(parts[0])
        for atom in self.atoms:
            self.masses.append(self.Masses[atom])
            

    def validate_output(self):
        """Validate the Gaussian output file."""
        return any('Normal termination of Gaussian' in line for line in self.lines)

    def parse_energies(self):
        """Parse the SCF energies from the output file."""
        for line in self.lines:
            if 'SCF Done:' in line:
                energy = float(line.split()[4])
                self.energies = np.append(self.energies, energy)

    def parse_coords(self):
        """Parse the coordinates from each step of the optimization."""
        reading_coordinates, skip_lines = False, 0
        step_coords = []  # Temporary list to hold coordinates for the current step
        self.coordinates = []  # Final list of steps, where each step is a list of coordinates
        
        for line in self.lines:
            if skip_lines > 0:
                skip_lines -= 1
                continue
            
            if 'Input orientation:' in line:
                if step_coords:  # Check if there are coordinates from a previous step
                    self.coordinates.append(np.array(step_coords).reshape(1, -1, 3))  # Add them to the coordinates list
                    step_coords = []  # Reset for the next step
                skip_lines = 4  # Skip the next 4 lines (headers)
                reading_coordinates = True
                continue
    
            if reading_coordinates:
                if "---------------------------------------------------------------------" in line:
                    reading_coordinates = False  # Stop reading coordinates
                    continue  # Skip to the next line in the file
    
                parts = line.strip().split()
                if len(parts) >= 6:  # Valid coordinate line
                    # Extract and store coordinates
                    step_coords.append(list(map(float, parts[3:6])))
        


# Example usage:
parser = GaussianParse('../Test_Gaussian_Output/but-2-ene.log')
# print("Energies:", parser.energies)
# print("Coordinates:", parser.coords)  # A list of numpy arrays, each representing a step
# print("Atoms:", parser.atoms)
# print("Masses:", parser.masses)

