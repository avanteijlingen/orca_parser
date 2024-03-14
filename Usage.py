# -*- coding: utf-8 -*-
"""
Created on Sat Dec  3 13:59:58 2022

@author: avtei
"""
import matplotlib.pyplot as plt
from orca_parser import *
import ase, sys


ams = ams_parse("Test-cases/AMS/133885_xyz.out")
ams.parse_energies()
ams.parse_coords()
print(ams.valid)
print(ams.energies)
print(ams.coords)
ams.asemol.write("Test-cases/AMS/133885_xyz.xyz")

sys.exit()

print("Test-cases/K_partial_crown.out")
K_partial_crown = ORCAParse("Test-cases/K_partial_crown.out", verbose=True)
K_partial_crown.parse_charges()
print(K_partial_crown.charges)

assert K_partial_crown.charges["Mulliken"][3] == -0.741417
assert K_partial_crown.charges["Mayer"][14] == 0.1226

print("Test-cases/Min_140.out")
Min_140 = ORCAParse("Test-cases/Min_140.out", verbose=True)
assert Min_140.valid, "This file is invalid"
print(Min_140.parse_input())

print("H2O")
H2O = ORCAParse("Test-cases/H2O.out")
assert H2O.valid, "This file is invalid"
print(H2O.parse_input())


print("Phenol")
Phenol = ORCAParse("Test-cases/Phenol/Opt.out")
print(Phenol.parse_input())

print("carbene CPCM_OPT")
CPCM_OPT = ORCAParse("Test-cases/cpcm_opt.out")
assert CPCM_OPT.valid, "This file is invalid"
print(CPCM_OPT.parse_input())

print("carbene gas-SP")
gas_SP = ORCAParse("Test-cases/gasSP-3.out")
assert gas_SP.valid, "This file is invalid"
print(gas_SP.parse_input())

print("Cl2")
Cl2 = ORCAParse("Test-cases/Cl2.out")
assert Cl2.valid, "This file is invalid"
print(Cl2.parse_input())

print("Coordination_0")
Coordination_0 = ORCAParse("Test-cases/Coordination_0.out")
assert Coordination_0.valid, "This file is invalid"
print(Coordination_0.parse_input())


H2O.parse_energies()
H2O.parse_free_energy()

print(H2O.energies)
print(H2O.AllGibbs)
print(H2O.entropies)
print(H2O.enthalpies)



#Hess = HessianTools("Test-cases/COO/COO.hess")
Hess = HessianTools("Test-cases/Coordination_0.hess")
#print(Hess.normalmodes)
#Hess.normalmodes.to_csv("Test.csv")
print(Hess.IR.iloc[20])
Hess.WriteMode("NormalMode_20.xyz", 20, steps=50)



print("Phenol")
Optimization = ORCAParse("Test-cases/Phenol/Opt.out")

print("ORCA exited normally:", Optimization.valid)
print("Job took:", Optimization.seconds(), "seconds")
print("Job input line:", Optimization.parse_input())



Optimization.parse_coords()
print("Atoms:", Optimization.atoms)
print("Final coordinates:")
print(Optimization.coords[-1])


Optimization.parse_energies()
print("Energy at each step:", Optimization.energies)
print("Energy at convergence:", Optimization.energies[-1], "Ha")

TDDFT = ORCAParse("Test-cases/Phenol/TDDFT.out")

TDDFT.parse_absorption()
print("Wavelengths:", TDDFT.wavelengths)
TDDFT.parse_CD()
print("CD Wavelengths:", TDDFT.CD)
X, Y = [],[]
for nm,r in zip(TDDFT.CD, TDDFT.R):
    X.append(nm)
    Y.append(0)
    X.append(nm)
    Y.append(r)
    X.append(nm)
    Y.append(0)
plt.plot(X, Y)
plt.xlabel("Wavelength (nm)")

print("\n")
print("Meisenheimer Complex")
Optimization = ORCAParse("Test-cases/MeisenheimerComplex.out")
Optimization.parse_dispersion()
print("Meisenheimer Complex D4 energy:", Optimization.dispersions[-1])
Optimization.parse_free_energy()
print("Meisenheimer Complex Gibbs free energy:", Optimization.Gibbs)
Optimization.parse_freqs()
print("Meisenheimer Complex Frequencies:", Optimization.frequencies)
print("Meisenheimer Complex entropies:", Optimization.entropies)
print("Meisenheimer Complex enthalpies:", Optimization.enthalpies)



### We can use energy warnings to check for things like a wavefunction not being fully converged:
notconv = ORCAParse("Test-cases/cpcm_opt.out")
notconv.parse_energies()
print(notconv.energies.shape[0])
if any(notconv.energy_warnings):
    print("Unconverged energies found, removing them")
print(notconv.energies[~notconv.energy_warnings].shape[0])



print("Parse Dipole Moment:")
Dipole = ORCAParse("Test-cases/Dipole.out")
Dipole.parse_dipole()



print("ASE mol:")

print(Dipole.asemol)