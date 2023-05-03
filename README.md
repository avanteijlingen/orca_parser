# ORCA-Parser
A python module to parse data out of ORCA output files

###Requirements:
[ase](https://gitlab.com/ase/ase), numpy, matplotlib, pandas

#Orca parser contains two classes:

####ORCAParse
Reads orca output streams (.out/.log/etc)
Can read frequencies, atoms, coordinates, IR spectra, free energy (broken down into its components as well)
Can also tell you how the job finished, if it converged etc

####HessianTools
Reads orca Hessian (.hess) outputs
Can parse the atoms and coordinates, normal modes, IR spectra.
From this .xyz trajectories of normal modes can be written

