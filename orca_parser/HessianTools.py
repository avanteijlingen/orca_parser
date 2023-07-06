# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 10:34:45 2023

@author: Alex
"""
import pandas
from orca_parser.ORCAParse import *


class HessianTools(ORCAParse):
    def parseHess(self):
        for section in self.raw.split("$")[1:]:
            section = section.split("\n")
            self.content[section[0]] = "\n".join([x for x in section[1:] if len(x) != 0])
            
        self.getAtoms()
        self.getNormalModes()
        self.getSpectra()
    #def getASEmol(self):
        #pass
    def getAtoms(self):
        self.Natoms = int(self.content["atoms"].split("\n")[0][0])
        self.atoms = []
        self.positions = np.ndarray((0, 3))
        for line in self.content["atoms"].split("\n")[1:]:
            line = line.split()
            self.atoms.append(line[0])
            self.positions = np.vstack((self.positions, [float(x) for x in line[2:]]))
        self.positions *= 0.52917724900001 # Bohr -> Angstrom
        self.mol = ase.Atoms(self.atoms, self.positions)
        
    def getNormalModes(self):
        basecol = 0
        self.Nnormalmodes = self.Natoms = int(self.content["normal_modes"].split("\n")[0].split()[0])
        #self.normalmodes = pandas.DataFrame()
        self.normalmodes = np.ones((self.Nnormalmodes, self.Nnormalmodes))
        for line in self.content["normal_modes"].split("\n"):
            if line[0] == "#":
                continue
            line = line.split()
            try:
                [int(x) for x in line]
                continue
            except:
                pass
            
            row = int(line[0])
            #print(row, basecol, line)
            for col in range(0, len(line[1:])):
                self.normalmodes[row, col+basecol] = float(line[col+1])*0.52917724900001 # Bohr -> Angstrom
            
            if row == self.Nnormalmodes-1:
                basecol += 5
        self.normalmodes = pandas.DataFrame(self.normalmodes)
    def getSpectra(self):
        self.IR = pandas.DataFrame(columns="wavenumber eps Int TX TY TZ".split())
        i = 0
        for line in self.content["ir_spectrum"].split("\n"):
            line=line.split()
            if len(line) < len(self.IR.columns):
                continue
            line = [float(x) for x in line]
            self.IR.loc[i] = line
            i+=1
    def WriteMode(self, fname, mode, steps=10, amplitude=1):
        self.mol.write(fname, append=False)
        disp = self.mol.copy()
        #yield steps < 1
        if steps == 1:
            steps = [1]
        else:
            steps = np.linspace(0, 1, steps)
        for direction in [1, -1]:
            for step in steps:
                disp.positions = self.mol.positions + (self.normalmodes[mode].values.reshape(-1,3) * step * direction * amplitude)
                disp.write(fname, append=True)
            # Reverse direction to make it smoothly back and forth
            for step in np.flip(steps):
                disp.positions = self.mol.positions + (self.normalmodes[mode].values.reshape(-1,3) * step * direction * amplitude)
                disp.write(fname, append=True)
        
    def __init__(self, fname):
        self.raw = readin(fname)
        self.content = {}
        self.parseHess()