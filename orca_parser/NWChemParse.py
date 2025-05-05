#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on [Current Date]

@author: [Your Name]
"""

import numpy as np
import ase
from . import ORCAParse
from .ORCAParse import tricky_readin


class NWChemParse(ORCAParse):
    def ValidateOutput(self):
        """Validate the NWChem output file."""
        self.valid = False
        if "Total times" in self.raw:  # NWChem typically ends with timing information
            if "Task times" in self.raw:
                self.valid = True

    def __init__(self, fname, verbose=False):
        super().__init__(fname, verbose=verbose)
        self.filepath = fname
        self.raw = tricky_readin(fname)
        if self.validate_output():
            self.valid = True
            if self.verbose:
                print("File terminated Normally")
        else:
            self.valid = False
            if self.verbose:
                print(
                    "The NWChem output file did not terminate normally or is not valid."
                )

    def validate_output(self):
        """Validate the NWChem output file."""
        return "Total times" in self.raw

    def parse_coords(self):
        self.parse()

    def parse(self):
        """
        Parse the coordinates from each step of the optimization.
        In an optimization the first 2 "frames" are the same
        and the last 2 "frames" are the same
        """
        protons = {v: k for k, v in ase.data.atomic_numbers.items()}
        self.coords = []
        self.atoms = []
        self.energies = []
        self.atomic_numbers = []

        # Get energy
        for line in self.raw.split("\n"):
            if len(line) == 0:
                continue

            if line[0] == "@":
                if line.split()[2] == "Energy" or line.split()[2] == "----------------":
                    continue
                E = float(line.split()[2])
                self.energies.append(E)

        blocks = self.raw.split("Output coordinates in angstroms")
        for iblock, block in enumerate(blocks[1:]):
            positions = np.ndarray((0, 3))
            block = (
                block.split("Atomic Mass")[0].split("--------------")[-1].split("\n")
            )
            for line in block:
                line = line.split()
                if len(line) < 4:
                    continue
                index, element, mass, x, y, z = line
                if iblock == 0:
                    self.atoms.append(element)
                    self.atomic_numbers.append(mass)
                xyz = np.array([float(x), float(y), float(z)]).reshape(1, 3)
                positions = np.vstack((positions, xyz))
            self.coords.append(positions)
        self.coords = np.array(self.coords)
        if (self.coords[0] == self.coords[1]).all():
            self.coords = self.coords[1:]
        if (self.coords[-2] == self.coords[-1]).all():
            self.coords = self.coords[:-1]
            self.energies = self.energies[:-1]
