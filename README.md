# orca_parser
A python module to parse data out of ORCA output files

This module is availible at: https://pypi.org/project/orca-parser/
and can be installed *via* 

```bash
pip install orca-parser
```

The module use as:
```python
import orca_parser
Optimization = orca_parser.ORCAParse("Test-cases/Phenol/Opt.out")

print("ORCA exited normally:", Optimization.valid)
print("Job took:", Optimization.seconds(), "seconds")
print("Job input line:", Optimization.parse_input())

Optimization.parse_coords()
print("Atoms:", Optimization.atoms)
print("Final coordinates:")
print(Optimization.coords[-1])
```


### Requirements:
[ase](https://gitlab.com/ase/ase), numpy, pandas

# Orca parser contains different classes for different types of QM output:

### ORCAParse (Orca 5+)
Reads orca output streams (.out/.log/.hess)
Can read frequencies, atoms, coordinates, IR spectra, free energy (broken down into its components as well)
Can also tell you how the job finished, if it converged etc

Reads:
 - Energies & Free Energies
 - Forces
 - Coordinates
 - Spectra
 - Dipoles
 - etc

#### HessianTools
Reads orca Hessian (.hess) outputs
Can parse the atoms and coordinates, normal modes, IR spectra.
From this .xyz trajectories of normal modes can be written

### GaussianParse (Gaussian 16) 
Reads Gaussian output streams (.out, .log etc)

Reads:
 - Energies & Free Energies
 - Forces
 - Coordinates

### ams_parse (Amsterdam Modeling Suite)

Reads:
 - Energies 
 - Coordinates

### NWChem (NWChem 7+)

Introduced the .parse() function for this one so you dont have to do .parse_coords() .parse_energies() etc

Reads:
 - Energies 
 - Coordinates


## Usage:
Example usage can be found in Usage.py, along with example ORCA output files in the "Test-cases" folder.

```java
                                                                               
 ██████╗ ██████╗  ██████╗ █████╗         ██████╗  █████╗ ██████╗ ███████╗███████╗██████╗ 
██╔═══██╗██╔══██╗██╔════╝██╔══██╗        ██╔══██╗██╔══██╗██╔══██╗██╔════╝██╔════╝██╔══██╗
██║   ██║██████╔╝██║     ███████║        ██████╔╝███████║██████╔╝███████╗█████╗  ██████╔╝
██║   ██║██╔══██╗██║     ██╔══██║        ██╔═══╝ ██╔══██║██╔══██╗╚════██║██╔══╝  ██╔══██╗
╚██████╔╝██║  ██║╚██████╗██║  ██║███████╗██║     ██║  ██║██║  ██║███████║███████╗██║  ██║
 ╚═════╝ ╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝╚══════╝╚═╝     ╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝╚══════╝╚═╝  ╚═╝
                                                                                         
```