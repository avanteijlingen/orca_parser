#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 15:20:02 2024

@author: bwb16179
"""

import numpy as np
import sys, ase
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
                
    def parse_input(self):
        self.Z = int(self.raw.split("Charge = ")[1].split("Multiplicity")[0])
        self.Multiplicity = int(self.raw.split("Multiplicity = ")[1].split("\n")[0].strip())
        if "Optimization completed" in self.raw :
            Job = "OPT"
        elif "IRC-IRC-IRC" in self.raw:
            Job = "IRC"
        else:
            Job = "SP"
        
        if "Solvent              :" in self.raw:
            Solv = self.raw.split("Solvent              :")[1].strip().split()[0].replace(",", "")
        else:
            Solv = "Gas"
            
        inp_dict = {"Job": Job,
                    "Multiplicity": self.Multiplicity,
                    "Charge": self.Z,
                    "software": "gaussian",
                    "Freq": "Harmonic frequencies" in self.raw,
                    "Functional": self.raw.split("Ex+Corr=")[1].split(" ExCW")[0],
                    "version": self.raw.split("*********************\n Gaussian")[1].split(":")[0].split()[0],
                    "BasisSet": self.raw.split("Standard basis:")[1].split("\n")[0].strip().split()[0],
                    "Solvation": Solv,
                    "def2J": None,
                    "RIJCOSX": None,
                    }
        
        return inp_dict

    def _read_lines(self):
        """Read the file and store lines for further processing."""
        with open(self.filepath, 'r') as file:
            lines = file.readlines()
        return lines
    

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
                
    def parse_forces(self):
        """Parse the SCF energies from the output file."""
        self.parse_coords()
        self.forces = np.ndarray((0, self.coords.shape[1], 3))
        reading_forces = False
        for line in self.lines:
            if not reading_forces:
                if "Forces (Hartrees/Bohr)" in line:
                    reading_forces = True
                    F = []
                    continue
                    
            if reading_forces:
                if "Cartesian Forces:" in line:
                    reading_forces = False
                    self.forces = np.vstack((self.forces, np.array(F).reshape(1,-1,3)))
                    continue
                elif "--------------" in line:
                    continue
                elif "Number" in line:
                    continue
                elif len(line.split()) == 5:
                    line = line.split()
                    x,y,z = line[2:]
                    xyz = [float(x), float(y), float(z)]
                    F.append(xyz)
                
    
    def parse_free_energy(self):
        self.AllGibbs = {}
        self.entropies = {}
        self.enthalpies = {}
        if "Sum of electronic and thermal Free Energies=" not in self.raw:
            return None
        
        for i,block in enumerate(self.coord_blocks):
            if "Sum of electronic and thermal Free Energies=" in block:
                G = float(block.split("Sum of electronic and thermal Free Energies=")[1].split("\n")[0])
                self.AllGibbs[i-1] = G

        self.Gibbs = list(self.AllGibbs.values())[-1]
        
        
    def parse_coords(self):
        """Parse the coordinates from each step of the optimization."""
        reading_coordinates, skip_lines = False, 0
        self.coords = np.ndarray((0,3))
        
        self.split_str = "Standard orientation:" if "Standard orientation:" in self.raw else "Input orientation:"
        
        protons = {v: k for k, v in ase.data.atomic_numbers.items()} 
        
        self.coord_blocks = self.raw.split(self.split_str)[1:]
        for block in self.coord_blocks:
            self.atoms = []
            self.atomic_numbers = []
            block = block.split("---------------------------------------------------------------------")[2]
            for line in block.split("\n"):
                line = line.split()
                if len(line) != 6:
                    continue
                self.coords = np.vstack((  self.coords, [float(x) for x in line[3:]] ))           
                self.atomic_numbers.append(int(line[1]))
                self.atoms.append(protons[int(line[1])])
                
        
        natoms = len(self.atoms)
        self.coords = self.coords.reshape((-1, natoms, 3))
        # Sometimes the last one is doubled up in the file
        if self.coords.shape[0] > 1:
            if (self.coords[-1] == self.coords[-2]).all():
                self.coords = np.delete(self.coords, [-1], axis=0)
        