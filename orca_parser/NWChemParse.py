#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import ase
from . import ORCAParse
from .ORCAParse import tricky_readin


class NWChemParse(ORCAParse):
    def ValidateOutput(self):
        """Validate the NWChem output file."""
        self.valid = True
        self.clean_stop = False
        if "Total times" in self.raw:  # NWChem typically ends with timing information
            if "Task times" in self.raw:
                self.clean_stop = True

    def __init__(self, fname, verbose=False):
        super().__init__(fname, verbose=verbose)
        self.filepath = fname
        self.raw = tricky_readin(fname)
        self.net_charge = None
        self.ValidateOutput()

    def parse_coords(self):
        self.parse()

    def parse_COSMO(self):
        self.cosmo_energies = []
        for line in self.raw.split("\n"):
            if "COSMO energy =" in line:
                self.cosmo_energies.append(
                    float(line.split("COSMO energy =")[1].split("\n")[0])
                )

    def parse_net_charge(self):
        return int(self.raw.split("Charge           :")[1].split("\n")[0].strip())

    def parse_multiplicity(self):
        return int(self.raw.split("Spin multiplicity:")[1].split("\n")[0].strip())

    def isSMD(self):
        return "COSMO-SMD solvation results" in self.raw

    def get_functional(self):
        if "MN15 Method XC Functional" in self.raw:
            return "MN15"
        elif "PBE0 Method XC Functional" in self.raw:
            return "PBE0"
        else:
            return None

    def get_basis_sets(self):
        self.RaBasis = None
        self.orgBasis = None
        if "Ra                        crenbl" in self.raw:
            self.RaBasis = "crenbl_ecp"
        if "   def2-SVP  " in self.raw:
            self.orgBasis = "def2-SVP"
        if "   def2-TZVPP  " in self.raw:
            self.orgBasis = "def2-TZVPP"

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
        self.parse_COSMO()
        self.net_charge = self.parse_net_charge()
        self.multiplicity = self.parse_multiplicity()
        self.SMD = self.isSMD()
        self.functional = self.get_functional()
        self.get_basis_sets()

        if "@" in self.raw:
            # This is an optimization

            # Get energy
            for line in self.raw.split("\n"):
                if len(line) == 0:
                    continue

                if line[0] == "@":
                    if (
                        line.split()[2] == "Energy"
                        or line.split()[2] == "----------------"
                    ):
                        continue
                    E = float(line.split()[2])
                    self.energies.append(E)

            blocks = self.raw.split("Output coordinates in angstroms")
            for iblock, block in enumerate(blocks[1:]):
                positions = np.ndarray((0, 3))
                block = (
                    block.split("Atomic Mass")[0]
                    .split("--------------")[-1]
                    .split("\n")
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
            if len(self.coords) > 1:
                if (self.coords[0] == self.coords[1]).all():
                    self.coords = self.coords[1:]
                if (self.coords[-2] == self.coords[-1]).all():
                    self.coords = self.coords[:-1]
                    self.energies = self.energies[:-1]
        else:
            # This might be a single point energy
            self.energies.append(
                float(self.raw.split("Total DFT energy =")[1].split("\n")[0])
            )
            positions = np.ndarray((0, 3))
            for line in (
                self.raw.split("XYZ format geometry")[1].split("Basis")[0].split("\n")
            ):
                line = line.strip().split()
                if len(line) != 4:
                    continue
                element, x, y, z = line
                self.atoms.append(element)
                xyz = np.array([float(x), float(y), float(z)]).reshape(1, 3)
                positions = np.vstack((positions, xyz))
            self.coords.append(positions)

    def get_data(self):
        if self.net_charge is None:
            self.parse()
        return {
            "charge": self.net_charge,
            "multiplicity": self.multiplicity,
            "SMD": self.SMD,
            "functional": self.functional,
            "orgBasis": self.orgBasis,  # This should be the original basis set used for the calculation
            "RaBasis": self.RaBasis,  # This should be the auxiliary basis set used for the calculation
        }
