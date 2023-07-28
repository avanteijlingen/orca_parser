# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 10:37:13 2023

@author: Alex
"""
import numpy as np
import ase, pandas
from ase.io import read

def readin(fname):
    # ORCA sometimes makes output that is hard to parse, we will try to read it the quick way first
    try:
        f = open(fname)
        content = f.read()
        f.close()
        return content
    except:
        content = ""
        with open(fname, 'rb') as f:
            for line in f:
                content = content + line.decode("utf-8", errors='ignore') + "\n"
        return content


def fit_rms(ref_c,c):
    # move geometric center to the origin
    ref_trans = np.average(ref_c, axis=0)
    ref_c = ref_c - ref_trans
    c_trans = np.average(c, axis=0)
    c = c - c_trans

    # covariance matrix
    C = np.dot(c.T, ref_c)

    # Singular Value Decomposition
    (r1, s, r2) = np.linalg.svd(C)

    # compute sign (remove mirroring)
    if np.linalg.det(C) < 0:
        r2[2,:] *= -1.0
    U = np.dot(r1, r2)
    return (c_trans, U, ref_trans)

def calc_rmsd(c1, c2):
    rmsd = 0.0
    c_trans, U, ref_trans = fit_rms(c1, c2)
    new_c2 = np.dot(c2 - c_trans, U) + ref_trans
    rmsd = np.sqrt( np.average( np.sum( ( c1 - new_c2 )**2, axis=1 ) ) )
    return rmsd


class ORCAParse:
    def ValidateOutput(self):
        if "TOTAL RUN TIME:" in self.raw:
            self.time = self.raw.split("TOTAL RUN TIME:")[1]
            T = self.time.strip().split()
            self.time  = int(T[0]) * 24 *60 * 60
            self.time += int(T[2]) * 60 * 60
            self.time += int(T[4]) * 60
            self.time += int(T[6]) 
            self.time += int(T[8]) / 1000
            
        if not "***ORCA TERMINATED NORMALLY***" in self.raw:
            if self.verbose:
                print(self.fname, "orca did not terminate normally!")
            self.valid = False
        elif "The optimization did not converge but reached the maximum number of" in self.raw:
            if self.verbose:
                print(self.fname, "hit geom MaxIter!")
            self.valid = False
        else:
            self.valid = True

        if "THE OPTIMIZATION HAS CONVERGED" in self.raw:  
            self.CONVERGED = True
        else:
            self.CONVERGED = False
            
    def thermodynamics(self):
        #Search for the INNER ENERGY section
        pass
    
    def parse_energies(self):
        self.energies = np.ndarray((0,), np.float64)
        self.energy_warnings = np.ndarray((0,), np.bool_)
        for part in self.raw.split("FINAL SINGLE POINT ENERGY")[1:]:
            part = part.split("\n")[0].strip()
            if "(Wavefunction not fully converged!)" in part:
                part = part.split()[0]
                self.energy_warnings = np.hstack((self.energy_warnings, [True]))
            else:
                self.energy_warnings = np.hstack((self.energy_warnings, [False]))
            part = float(part)
            self.energies = np.hstack((self.energies, [part]))
        #self.r_energies = self.energies - self.energies.min()
    def parse_dispersion(self):
        splits = self.raw.split("\nDispersion correction")[1:]
        self.dispersions = np.ndarray((len(splits),))
        for i in range(len(splits)):
            E_disp = splits[i].split("\n")[0].strip()
            if "Starting D4" in E_disp: #Not a results line
                continue
            self.dispersions[i] = float(E_disp)
        
    def parse_coords(self):
        self.coords = []
        self.atoms = []
        frames = self.raw.split("CARTESIAN COORDINATES (ANGSTROEM)")[1:]
        for i,frame in enumerate(frames):
            positions = []
            frame = frame.split("CARTESIAN COORDINATES (A.U.)")[0]
            for line in frame.split("\n"):
                line = line.split()
                if len(line) != 4:
                    continue
                positions.append([float(x) for x in line[1:]])
                if i == 0:
                    self.atoms.append(line[0])
            if i == 0:
                self.coords = np.array(positions).reshape(1, -1, 3)
            else:
                self.coords = np.vstack((self.coords, np.array(positions).reshape(1,-1,3)))
    
    def scan_bond(self, a0, a1):
        distances = np.ndarray((self.coords.shape[0]))
        for i in range(self.coords.shape[0]):
            positions = self.coords[i]
            d = np.linalg.norm(positions[a0] - positions[a1])
            distances[i] = d
        return distances
            
    def parse_freqs(self):
        self.frequencies = []
        frames = self.raw.split("VIBRATIONAL FREQUENCIES")[1:]
        for i,frame in enumerate(frames):
            frequencies = []
            frame = frame.split("NORMAL MODES")[0]
            for line in frame.split("\n"):
                if "cm" not in line:
                    continue
                line = line.split()
                freq = float(line[1])
                frequencies.append(freq)
            if i == 0:
                self.frequencies = np.array(frequencies).reshape(1, -1)
            else:
                self.frequencies = np.vstack((self.frequencies, np.array(frequencies).reshape(1,-1)))
    
    def parse_IR(self):
        self.IR = pandas.DataFrame(columns=['freq', 'eps', 'Int', 'T**2', 'TX', 'TY', 'TZ'])
        ir = self.raw.split("IR SPECTRUM")[-1].split("--------------------------")[2]
        for line in ir.split("\n"):
            line = line.replace(":", "").replace("(", "").replace(")", "")
            line=line.split()
            if len(line) != 8:
                continue
            #print(len(line), line)
            try:
                self.IR.loc[line[0]] = [float(x) for x in line[1:]]
            except ValueError:
                pass
        
    def parse_free_energy(self):
        if "out the resulting rotational entropy values for" in self.raw:
            symmetric_number = self.raw.split("Final Gibbs free energy         ...")[-2].split("out the resulting rotational entropy values for")[1]
            symmetric_number = symmetric_number.split("--------------------------------------------------------")[1]
            self.symmetry_number = symmetric_number.count("sn=")
        else:
            self.symmetry_number = 2
        
        self.AllGibbs = [float(x.split("\n")[0].strip().split()[0]) for x in self.raw.split("Final Gibbs free energy         ...")[1:]]
        self.Gibbs = float(self.raw.split("Final Gibbs free energy         ...")[-1].split("\n")[0].strip().split()[0])
        self.entropies = [float(x.split("\n")[0].strip().split()[0]) for x in self.raw.split("Total enthalpy                    ...")[1:]]
        self.enthalpies = [float(x.split("\n")[0].strip().split()[0]) for x in self.raw.split("Total entropy correction          ...")[1:]]
        
        
        
    def seconds(self):
        time_str = self.raw.split("TOTAL RUN TIME:")[1].strip().split()
        days = int(time_str[0])
        hours = int(time_str[2])
        minutes = int(time_str[4])
        seconds = int(time_str[6])
        miliseconds = int(time_str[8])
        
        hours = hours + (days*24)
        minutes = minutes + (hours*60)
        seconds = seconds + (minutes*60)
        seconds = seconds + (miliseconds/1000)
        return seconds
    
    def parse_input(self):
        self.Z = int(round(float(self.raw.split("Sum of atomic charges:")[1].split("\n")[0])))
        self.Multiplicity = int(round(float(self.raw.split("* xyz")[1].split("\n")[0].split()[1])))
        inp = self.raw.split("|  1>")[1].split("\n")[0]
        inp = inp.upper()
        print(inp)
        inp = inp.replace("!", "")
        inp = inp.split()
        inp_dict = {"Job": None, 
                    "Freq": False,
                    "Solvation": "Gas",
                    "Dispersion": None,
                    "Charge": self.Z,
                    "defgrid": None,
                    "def2/J": None,
                    "Multiplicity": self.Multiplicity}
        i = 0
        while i < len(inp):
            if inp[i] in ["SP", "OPT"]:
                inp_dict["Job"] = inp[i]
                del inp[i]
                continue
            elif "FREQ" in inp[i]:
                inp_dict["Freq"] = True
                del inp[i]
                continue
            elif inp[i] in ["B3LYP", "PBE"] or "WB9" in inp[i]:
                inp_dict["Functional"] = inp[i]
                del inp[i]
                continue
            elif "CPCM" in inp[i] or "SMD" in inp[i]:
                inp_dict["Solvation"] = inp[i]
                del inp[i]
                continue
            elif "DEF2-" in inp[i] or "SVP" in inp[i] or "VP" in inp[i]:
                inp_dict["BasisSet"] = inp[i]
                del inp[i]
                continue
            elif "D3" in inp[i] or "D3BJ" in inp[i] or "D4" in inp[i]:
                inp_dict["Dispersion"] = inp[i]
                del inp[i]
                continue
            elif "DEFGRID" in inp[i]:
                inp_dict["defgrid"] = inp[i]
                del inp[i]
                continue
            elif "DEF2/J" in inp[i]:
                inp_dict["def2/J"] = True
                del inp[i]
                continue
            
            
            i+= 1
        return inp_dict
    
    def convergence(self):
        self.conv = dict()
        self.tol = dict()
        if "Geometry convergence" not in self.raw:
            return 0
        for segment in self.raw.split("Geometry convergence")[1:]:
            segment = segment.split("---------------------------------------------------------------------")[1]
            for line in segment.split("\n"):
                if ".................." in line:
                    continue
                line = line.strip().split()
                if len(line) == 0:
                    continue
                
                if "(" in line[0]:
                    if line[0] not in list(self.conv.keys()):
                        self.conv[line[0]] = []
                    if line[2] not in list(self.conv.keys()):
                        self.conv[line[2]] = []
                    
                    self.conv[line[0]].append(float(line[1]))
                    self.conv[line[2]].append(float(line[3]))
                else:
                    key = " ".join(line[:2])
                    if key not in list(self.conv.keys()):
                        self.conv[key] = []
                        self.tol[key] = []
                    self.conv[key].append(float(line[2]))
                    self.tol[key].append(float(line[3]))
    
                #print(line)
            
    def parse_absorption(self):
        txt = self.raw.split("ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS")[1]
        txt = txt.split("-----------------------------------------------------------------------------")[2]
        self.wavelengths = []
        self.fosc = []
        for line in txt.split("\n"):
            line=line.split()
            if len(line) > 4:
                self.wavelengths.append(float(line[2]))
                self.fosc.append(float(line[3]))
    
    def parse_CD(self):
        txt = self.raw.split("CD SPECTRUM")[1]
        txt = txt.split("-------------------------------------------------------------------")[2]
        self.CD = []
        self.R = []
        for line in txt.split("\n"):
            line=line.split()
            if len(line) > 4:
                self.CD.append(float(line[2]))
                self.R.append(float(line[3]))
    
    
    def __init__(self, fname, verbose = False):
        self.fname = fname
        self.verbose = verbose
        self.raw = readin(fname)
        self.ValidateOutput()
        self.convergence()
        self.TDDFT = False