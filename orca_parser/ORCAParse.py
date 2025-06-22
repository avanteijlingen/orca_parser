# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 10:37:13 2023

@author: Alex
"""
import numpy as np
import ase, pandas, os
from ase.io import read
from ase import Atoms


def tricky_readin(fname):
    # First try with the repr function
    assert os.path.exists(fname), f"{fname} not found"

    # =============================================================================
    #     with open(fname, "r", encoding="utf8", errors='ignore') as f:
    #         content = f.read()
    #     print(content[-100:])
    # =============================================================================
    # first try with defaults
    try:
        with open(fname, "r") as f:
            content = f.read()
        assert "Frank Neese" in content
    except:
        # next try a linux / archie output
        try:
            with open(fname, "rb") as f:
                content = f.read().decode("utf-8")
            # print("Encoded as UTF-8")
        except:
            # next try a winwdows 10/11 encoding
            try:
                with open(fname, "rb") as f:
                    content = f.read().decode("utf-16")
                # print("Encoded as UTF-16")
                # Need to change \r\n to \n
                content = content.replace("\r\n", "\n")
                content = content.replace("\r", "")
            except:
                content = f"{fname}: Couldnt decode output file as utf-8 or utf-16"
                print(content)
    return content.replace("\r\n", "\n")


def fit_rms(ref_c, c):
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
        r2[2, :] *= -1.0
    U = np.dot(r1, r2)
    return (c_trans, U, ref_trans)


def calc_rmsd(c1, c2):
    rmsd = 0.0
    c_trans, U, ref_trans = fit_rms(c1, c2)
    new_c2 = np.dot(c2 - c_trans, U) + ref_trans
    rmsd = np.sqrt(np.average(np.sum((c1 - new_c2) ** 2, axis=1)))
    return rmsd


class ORCAParse:
    @property
    def asemol(self):
        """
        Returns
        -------
        ase.Atoms
            An ase molecule with the coordinates of the last conformer in the orca output.

        """
        self.parse_coords()
        # print(len(self.atoms), len(self.coords[-1]))
        return self.makeASE(-1)

    def makeASE(self, frame: int):
        return Atoms(self.atoms, self.coords[frame])

    def ValidateOutput(self):
        if "TOTAL RUN TIME:" in self.raw:
            self.time = self.raw.split("TOTAL RUN TIME:")[1]
            T = self.time.strip().split()
            self.time = int(T[0]) * 24 * 60 * 60
            self.time += int(T[2]) * 60 * 60
            self.time += int(T[4]) * 60
            self.time += int(T[6])
            self.time += int(T[8]) / 1000

        if not "***ORCA TERMINATED NORMALLY***" in self.raw:
            if self.verbose:
                print(self.fname, "orca did not terminate normally!")
            self.valid = False
        elif (
            "The optimization did not converge but reached the maximum number of"
            in self.raw
        ):
            if self.verbose:
                print(self.fname, "hit geom MaxIter!")
            self.valid = False
        else:
            self.valid = True

        if "ORCA ab initio Molecular Dynamics Module" in self.raw:
            if self.verbose:
                print("Orca-parser does not do ab-initio jobs")
            self.valid = False

        if "this file is used internally by ORCA" in self.raw:
            self.valid = False
        elif (
            "this file is use internally by ORCA" in self.raw
        ):  # There is a spell mistake we need to account for in older versions of ORCA
            self.valid = False

        if "THE OPTIMIZATION HAS CONVERGED" in self.raw:
            self.CONVERGED = True
        else:
            self.CONVERGED = False

    def thermodynamics(self):
        # Search for the INNER ENERGY section
        pass

    def parse(self):
        if "Global Optimization Algorithm" in self.raw:
            self.GOAT = True
        else:
            self.GOAT = False
        if "xtb is free software" in self.raw:
            self.XTB = True
        else:
            self.XTB = False

        if self.XTB:
            self.Z = float(self.raw.split(":: total charge")[1].split("e")[0])
            self.orca_version = self.raw.split("Program Version ")[1].split()[0]
        else:
            if "END OF INPUT" in self.raw:
                self.input_dict = self.parse_input()

        if "DIPOLE MOMENT" in self.raw:
            self.parse_dipole()
        if "FINAL SINGLE POINT ENERGY" in self.raw:
            self.parse_energies()
        if "Dispersion correction" in self.raw:
            self.parse_dispersion()
        if "CARTESIAN COORDINATES" in self.raw:
            self.parse_coords()
        if "VIBRATIONAL FREQUENCIES" in self.raw:
            self.parse_freqs()
        if "IR SPECTRUM" in self.raw:
            self.parse_IR()
        if "Gibbs free energy" in self.raw:
            self.parse_free_energy()
        if "ORBITAL ENERGIES" in self.raw:
            self.parse_HOMO_LUMO()

        if "convergence" in self.raw:
            self.convergence()
        if "ABSORPTION SPECTRUM VIA TRANSITION" in self.raw:
            self.parse_absorption()
        if "CD SPECTRUM" in self.raw:
            self.parse_CD()
        if "MULLIKEN ATOMIC CHARGES" in self.raw:
            self.parse_charges()

    def parse_dipole(self):
        assert "DIPOLE MOMENT" in self.raw
        dipoles = []
        for section in self.raw.split("DIPOLE MOMENT")[1:]:
            section = section.split("Rotational spectrum ")[0]
            dipole = {}

            for part in [
                "Electronic contribution",
                "Nuclear contribution",
                "Total Dipole Moment",
                "Magnitude (a.u.)",
                "Magnitude (Debye)",
            ]:
                try:
                    dipole[part] = (
                        section.split(part)[1].split(":")[1].split("\n")[0].split()
                    )
                    dipole[part] = [float(x) for x in dipole[part]]
                except:
                    dipole[part] = None

            dipoles.append(dipole)

        return dipoles

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
        # self.r_energies = self.energies - self.energies.min()

    def parse_dispersion(self):
        splits = self.raw.split("\nDispersion correction")[1:]
        self.dispersions = np.ndarray((len(splits),))
        for i in range(len(splits)):
            E_disp = splits[i].split("\n")[0].strip()
            if "Starting D4" in E_disp:  # Not a results line
                continue
            elif "... done" in E_disp:  # Not a results line
                continue
            self.dispersions[i] = float(E_disp)

    def parse_coords(self):
        self.coords = []
        self.atoms = []
        frames = self.raw.split("CARTESIAN COORDINATES (ANGSTROEM)")[1:]
        for i, frame in enumerate(frames):
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
                self.coords = np.vstack(
                    (self.coords, np.array(positions).reshape(1, -1, 3))
                )
        for x in self.atoms:
            self.masses.append(self.Masses[x])

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
        for i, frame in enumerate(frames):
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
                self.frequencies = np.vstack(
                    (self.frequencies, np.array(frequencies).reshape(1, -1))
                )

    def parse_IR(self):
        self.IR = pandas.DataFrame(
            columns=["freq", "eps", "Int", "T**2", "TX", "TY", "TZ"]
        )
        ir = self.raw.split("IR SPECTRUM")[-1].split("--------------------------")[2]
        for line in ir.split("\n"):
            line = line.replace(":", "").replace("(", "").replace(")", "")
            line = line.split()
            if len(line) != 8:
                continue
            # print(len(line), line)
            try:
                self.IR.loc[line[0]] = [float(x) for x in line[1:]]
            except ValueError:
                pass

    def parse_free_energy(self):
        if "out the resulting rotational entropy values for" in self.raw:
            symmetric_number = self.raw.split("Final Gibbs free energy         ...")[
                -2
            ].split("out the resulting rotational entropy values for")[1]
            symmetric_number = symmetric_number.split(
                "--------------------------------------------------------"
            )[1]
            self.symmetry_number = symmetric_number.count("sn=")
        else:
            self.symmetry_number = 2

        self.AllGibbs = {}
        self.entropies = {}
        self.enthalpies = {}

        if "Gibbs free energy" not in self.raw:
            return None

        frames = self.raw.split("GEOMETRY OPTIMIZATION CYCLE")

        for i, frame in enumerate(frames):
            if "Final Gibbs free energy" in frame:
                G = frame.split("Final Gibbs free energy         ...")[1].split("Eh")[0]
                H = frame.split("Total enthalpy                    ...")[1].split("Eh")[
                    0
                ]
                S = frame.split("Total entropy correction          ...")[1].split("Eh")[
                    0
                ]
                self.AllGibbs[i] = float(G)
                self.enthalpies[i] = float(H)
                self.entropies[i] = float(S)
        self.Gibbs = list(self.AllGibbs.values())[-1]

    def parse_HOMO_LUMO(self):
        frames = self.raw.split("ORBITAL ENERGIES")[1:]

        self.all_HOMO = []
        self.all_LUMO = []

        for i, frame in enumerate(frames):
            orbitals_occ = []
            orbitals_unocc = []
            frame = frame.split("                    * MULLIKEN POPULATION ANALYSIS *")[
                0
            ]
            for line in frame.split("\n"):
                line = line.split()
                if len(line) != 4:
                    continue
                if "NO" in line:
                    continue
                if "0.0000" in line:
                    orbitals_unocc.append(
                        [int(line[0]), float(line[1]), float(line[2]), float(line[3])]
                    )
                    continue
                else:
                    orbitals_occ.append(
                        [int(line[0]), float(line[1]), float(line[2]), float(line[3])]
                    )
                    continue
            self.all_HOMO.append(max(orbitals_occ, key=lambda sublist: sublist[0]))
            self.all_LUMO.append(min(orbitals_unocc, key=lambda sublist: sublist[0]))

        self.HOMO_LUMO_gap = self.all_HOMO[-1][2] - self.all_LUMO[-1][2]

    def seconds(self):
        time_str = self.raw.split("TOTAL RUN TIME:")[1].strip().split()
        days = int(time_str[0])
        hours = int(time_str[2])
        minutes = int(time_str[4])
        seconds = int(time_str[6])
        miliseconds = int(time_str[8])

        hours = hours + (days * 24)
        minutes = minutes + (hours * 60)
        seconds = seconds + (minutes * 60)
        seconds = seconds + (miliseconds / 1000)
        return seconds

    def parse_input(self):
        self.Z = int(
            self.raw.split("Total Charge           Charge          ....")[1].split(
                "\n"
            )[0]
        )
        self.Multiplicity = int(
            self.raw.split("Multiplicity           Mult            ....")[1].split(
                "\n"
            )[0]
        )
        self.orca_version = self.raw.split("Program Version ")[1].split()[0]
        inp = self.raw.split("INPUT FILE")[1].split("****END OF INPUT****")[0]
        inp = inp.split("> !")[1].split("\n")[0]
        inp = inp.upper()
        inp = inp.replace("!", "")
        inp_dict = {
            "Job": "OPT" if "Geometry Optimization Run" in self.raw else "SP",
            "BasisSet": self.raw.split("Your calculation utilizes the basis:")[1]
            .split("\n")[0]
            .replace(",", "")
            .rstrip()
            .replace(" ", ""),
            "Freq": "VIBRATIONAL FREQUENCIES" in self.raw,
            "RIJCOSX": "RIJ-COSX (HFX calculated with COS-X)).... on" in self.raw,
            "def2J": "Your calculation utilizes the auxiliary basis: def2/J"
            in self.raw,
            "Solvation": "Gas",
            "Charge": self.Z,
            "Multiplicity": self.Multiplicity,
            "version": self.orca_version,
            "software": "ORCA",
            "UKS": "UKS" in self.raw.split("|  1>")[1].split("\n")[0],
            "ECP": "ECP gradient" in self.raw,
        }

        # Finicky functionals
        if "HF-3C" in inp:
            inp_dict["Functional"] = "HF"
        elif (
            "Exchange Functional    Exchange        .... B88" in self.raw
            and "Correlation Functional Correlation     .... LYP" in self.raw
        ):
            inp_dict["Functional"] = "B3LYP"
        else:
            inp_dict["Functional"] = (
                self.raw.split("Exchange Functional    Exchange        ....")[1]
                .split("\n")[0]
                .rstrip()
                .replace(" ", "")
            )

        # Dispersion
        if "DFT DISPERSION CORRECTION" in self.raw:
            if "USING Becke-Johnson damping" in self.raw:
                inp_dict["Dispersion"] = "D3BJ"
            elif "Active option DFTDOPT                   ...         2" in self.raw:
                inp_dict["Dispersion"] = "D2"
            else:
                dispersion = self.raw.split("DFT DISPERSION CORRECTION")[1].split("\n")[
                    2
                ]
                dispersion = dispersion.split("DFT")[1].split(" ")[0]
                inp_dict["Dispersion"] = dispersion

        # Defgrid
        if "DEFGRID1" in inp:
            inp_dict["defgrid"] = "DEFGRID1"
        elif "DEFGRID3" in inp:
            inp_dict["defgrid"] = "DEFGRID3"
        else:
            inp_dict["defgrid"] = "DEFGRID2"

        # slowconv
        if "SLOWCONV" in inp:
            if "VERYSLOWCONV" in inp:
                inp_dict["SlowConv"] = "VERYSLOWCONV"
            else:
                inp_dict["SlowConv"] = "SLOWCONV"

        # SCF convergence
        if "SCF" in inp:
            scf = inp.split()
            scf = [x for x in scf if "SCF" in x][0]
            inp_dict["SCF_conv_tol"] = scf

        # Solvation
        if "CPCM SOLVATION MODEL" in self.raw:
            inp_dict["Solvation"] = "CPCM"
            eps = float(
                self.raw.split("Epsilon                                         ...")[1]
                .split("\n")[0]
                .strip()
            )
            refrac = float(
                self.raw.split("Refrac                                          ...")[1]
                .split("\n")[0]
                .strip()
            )
            if eps == 80.4 and refrac == 1.33:
                inp_dict["Solvation"] = "CPCM(WATER)"
            elif eps == 7.25 and refrac == 1.407:
                inp_dict["Solvation"] = "CPCM(THF)"
            elif eps == 47.2 and refrac == 1.479:
                inp_dict["Solvation"] = "CPCM(DMSO)"
            elif eps == 36.6 and refrac == 1.344:
                inp_dict["Solvation"] = "CPCM(Acetonitrile)"
            else:
                print(
                    f"Couldnt assign CPCM solvent from eps = {eps} and refrac = {refrac}"
                )

        return inp_dict

    def convergence(self):
        self.conv = dict()
        self.tol = dict()
        if "Geometry convergence" not in self.raw:
            return 0
        for segment in self.raw.split("Geometry convergence")[1:]:
            segment = segment.split(
                "---------------------------------------------------------------------"
            )[1]
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

                # print(line)

    def parse_absorption(self):
        txt = self.raw.split(
            "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"
        )[1]
        txt = txt.split(
            "-----------------------------------------------------------------------------"
        )[2]
        self.wavelengths = []
        self.fosc = []
        for line in txt.split("\n"):
            line = line.split()
            if len(line) > 4:
                self.wavelengths.append(float(line[2]))
                self.fosc.append(float(line[3]))

    def parse_CD(self):
        txt = self.raw.split("CD SPECTRUM")[1]
        txt = txt.split(
            "-------------------------------------------------------------------"
        )[2]
        self.CD = []
        self.R = []
        for line in txt.split("\n"):
            line = line.split()
            if len(line) > 4:
                self.CD.append(float(line[2]))
                self.R.append(float(line[3]))

    def parse_charges(self):
        self.charges = {"Mulliken": [], "Loewdin": [], "Mayer": []}
        txt = self.raw.split("MULLIKEN ATOMIC CHARGES")[-1]
        txt = txt.split("Sum of atomic charges")[0]
        for line in txt.split("\n"):
            if ":" in line:
                self.charges["Mulliken"].append(float(line.split()[-1]))
        txt = self.raw.split("LOEWDIN POPULATION ANALYSIS")[-1]
        txt = txt.split("LOEWDIN REDUCED ORBITAL CHARGES")[0]
        for line in txt.split("\n"):
            if ":" in line:
                self.charges["Loewdin"].append(float(line.split()[-1]))

        txt = self.raw.split("MAYER POPULATION ANALYSIS")[-1]
        txt = txt.split("Mayer bond orders larger than")[0]
        for line in txt.split("\n"):
            line = line.split()
            if len(line) == 8:
                self.charges["Mayer"].append(float(line[4]))

        self.charges["Mulliken"] = np.array(self.charges["Mulliken"])
        self.charges["Loewdin"] = np.array(self.charges["Loewdin"])
        self.charges["Mayer"] = np.array(self.charges["Mayer"])

    def __init__(self, fname, verbose=False):
        self.fname = fname
        self.verbose = verbose
        self.raw = tricky_readin(fname)
        if self.raw == "Couldnt decode output file as utf-8 or utf-16":
            print(self.raw)
            self.valid = False
            return None
        self.ValidateOutput()
        self.convergence()
        self.TDDFT = False
        self.coords = []
        self.atoms = []
        self.masses = []
        self.Masses = {
            "H": 1.008,
            "He": 4.003,
            "Li": 6.941,
            "Be": 9.012,
            "B": 10.811,
            "C": 12.011,
            "N": 14.007,
            "O": 15.999,
            "F": 18.998,
            "Ne": 20.180,
            "Na": 22.990,
            "Mg": 24.305,
            "Al": 26.982,
            "Si": 28.086,
            "P": 30.974,
            "S": 32.066,
            "Cl": 35.453,
            "Ar": 39.948,
            "K": 39.098,
            "Ca": 40.078,
            "Sc": 44.956,
            "Ti": 47.867,
            "V": 50.942,
            "Cr": 51.996,
            "Mn": 54.938,
            "Fe": 55.845,
            "Co": 58.933,
            "Ni": 58.693,
            "Cu": 63.546,
            "Zn": 65.38,
            "Ga": 69.723,
            "Ge": 72.631,
            "As": 74.922,
            "Se": 78.971,
            "Br": 79.904,
            "Kr": 84.798,
            "Rb": 84.468,
            "Sr": 87.62,
            "Y": 88.906,
            "Zr": 91.224,
            "Nb": 92.906,
            "Mo": 95.95,
            "Tc": 98.907,
            "Ru": 101.07,
            "Rh": 102.906,
            "Pd": 106.42,
            "Ag": 107.868,
            "Cd": 112.414,
            "In": 114.818,
            "Sn": 118.711,
            "Sb": 121.760,
            "Te": 126.7,
            "I": 126.904,
            "Xe": 131.294,
            "Cc": 132.905,
            "Ba": 137.328,
            "La": 138.905,
            "Ce": 140.116,
            "Pr": 140.908,
            "Nd": 144.243,
            "Pm": 144.913,
            "Sm": 150.36,
            "Eu": 151.964,
            "Gd": 157.25,
            "Tb": 158.925,
            "Dy": 162.500,
            "Ho": 164.930,
            "Er": 167.259,
            "Tm": 168.934,
            "Yb": 173.055,
            "Lu": 174.967,
            "Hf": 178.49,
            "Ta": 180.948,
            "W": 183.84,
            "Re": 186.207,
            "Os": 190.23,
            "Ir": 192.217,
            "Pt": 195.085,
            "Au": 196.967,
            "Hg": 200.592,
            "Tl": 204.383,
            "Pb": 207.2,
            "Bi": 208.980,
            "Po": 208.982,
            "At": 209.987,
            "Rn": 222.081,
            "Fr": 223.020,
            "Ra": 226.025,
            "Ac": 227.028,
            "Th": 232.038,
            "Pa": 231.036,
            "U": 238.029,
            "Np": 237,
            "Pu": 244,
            "Am": 243,
            "Cm": 247,
            "Bk": 247,
            "Ct": 251,
            "Es": 252,
            "Fm": 257,
            "Md": 258,
            "No": 259,
            "Lr": 262,
            "Rf": 261,
            "Db": 262,
            "Sg": 266,
            "Bh": 264,
            "Hs": 269,
            "Mt": 268,
            "Ds": 271,
            "Rg": 272,
            "Cn": 285,
            "Nh": 284,
            "Fl": 289,
            "Mc": 288,
            "Lv": 292,
            "Ts": 294,
            "Og": 294,
        }
