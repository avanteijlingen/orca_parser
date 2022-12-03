# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 17:42:17 2022

@author: avtei
"""
import numpy as np
import os, glob, sys, pandas
import matplotlib.pyplot as plt
from ase.io import read


# 1 Bohr = 0.52917724900001 Angstrom

def readin(fname):
    content = ""
    print(fname)
    with open(fname, 'rb') as f:
        for line in f:
            content = content + line.decode(errors='ignore') + "\n"
    
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
        for part in self.raw.split("FINAL SINGLE POINT ENERGY")[1:]:
            part = float(part.split("\n")[0].strip())
            self.energies = np.hstack((self.energies, [part]))
        #self.r_energies = self.energies - self.energies.min()
    def parse_dispersion(self):
        splits = self.raw.split("Dispersion correction")[1:]
        self.dispersions = np.ndarray((len(splits),))
        for i in range(len(splits)):
            E_disp = splits[i].split("\n")[0].strip()
            self.dispersions[i] = float(E_disp)
        
    def get_coords(self):
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
            
    def get_freqs(self):
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
    def parse_free_energy(self):
        if "out the resulting rotational entropy values for" in self.raw:
            symmetric_number = self.raw.split("Final Gibbs free energy         ...")[-2].split("out the resulting rotational entropy values for")[1]
            symmetric_number = symmetric_number.split("--------------------------------------------------------")[1]
            self.symmetry_number = symmetric_number.count("sn=")
        else:
            self.symmetry_number = 2
        
        self.AllGibbs = [float(x.split("\n")[0].strip().split()[0]) for x in self.raw.split("Final Gibbs free energy         ...")[1:]]
        self.Gibbs = float(self.raw.split("Final Gibbs free energy         ...")[-1].split("\n")[0].strip().split()[0])
        
    def time(self):
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
    
    def get_input(self):
        return self.raw.split("|  1>")[1].split("\n")[0]
    
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
            
    
    def __init__(self, fname, verbose = False):
        self.fname = fname
        self.verbose = verbose
        self.raw = readin(fname)
        self.ValidateOutput()
        self.convergence()


    
if __name__ == "__main__":
    #fname = "K-DHP/TS/Insertion/ORCA_SCAN/Insertion_TS.out"
    #fname = "K-DHP/TS/Insertion/ORCA_TS_fromDNN/DNN_Berny.out"
    #fname = "K-DHP/TS/Deprotonation/ORCA_SCAN/Deprotonation_ScanTS.out"
    
    #files = glob.glob("K-DHP_impTHF/TS/Trapping_tBu.ENBA/TZVPP_D4/IRC_forward/IRC_op*.out")
    #files = glob.glob("K-DHP_impTHF/TS/Hydride*/TZVPP_D4/*.out")
    #files = glob.glob("K-DHP/TS/*/OtherTS/*.out")
    #files = glob.glob("K-DHP/TS/*/ORCA_SCAN/*.out")
    #files = glob.glob("K-DHP/Opt/*/*.out")
    #files = glob.glob("K-DHP_impTHF/Opt/TZVPP_D4/KH/*.out")
    #files = ["/users/rkb19187/Desktop/Roitberg/ScanTS_OptTS/K-DHP_impTHF/TS/Deprotonation_tBu/TZVPP_D4/IRC_forward/IRC_opt.out"]
    #files = ["/users/rkb19187/Desktop/Roitberg/ScanTS_OptTS/Cs-DHP_impTHF/TS/Trapping_tBu.ENBA/TZVPP_D4/IRC_forward/IRC_opt.out"]
    files = ["/users/rkb19187/Desktop/Roitberg/ScanTS_OptTS/K-DHP_impTHF/TS/SurrogateHydrideTransfer_tBu/TZVPP_D4/SurrogateHydrideTransfer_tBu.out"]
    
    print(files)
    for TS_out in files:
        print(TS_out)

        op = ORCAParse(TS_out)
        op.get_coords()
        op.parse_energies()
        print("op.valid:", op.valid)
        #print("op.valid:", op.energies)
        #plt.plot(op.r_energies, label=os.path.basename(TS_out))
        #plt.show()
        
        
        op.get_freqs()
        freqs = op.frequencies[-1]
        print("Negative Freqs:", np.where(freqs < 0)[0].shape[0])
        
        op.parse_free_energy()
        print(op.Gibbs)
        
        print("CONVERGED:", op.CONVERGED)
        
        print(op.conv)
        print(op.tol)
        
# =============================================================================
#         fig, axs = plt.subplots(int(np.ceil(len(list(op.conv.keys())) / 2)),2)
#         fig.suptitle('Vertically stacked subplots')
#         for i,key in enumerate(list(op.conv.keys())):
#             index = int(np.floor(i / 2))
#             indey = i%2
#             axs[index][indey].set_title(key)
#             axs[index][indey].plot(op.conv[key])
# =============================================================================

    for key in list(op.conv.keys()):
        if key in list(op.tol.keys()):
            print(key)
            
            X = np.arange(len(op.conv[key]))
            plt.plot(X, op.conv[key])
            plt.plot(X, op.tol[key])
            plt.title(key)
            plt.show()
        
        

        
