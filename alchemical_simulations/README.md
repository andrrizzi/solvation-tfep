# Alchemical simulations
In this folder, you will find a Python script to automatically perform 
an absolute hydration free energy (AHFE) and a relative hydration free 
energy (RHFE) analysis from methane to iodomethane with two different 
software packages (`biosimspace` and `openfe`).

The first script, `bss_simulation.py`, uses the `biosimspace` simulation 
package from OpenBioSim, and performs a default RHFE simulation for 1
ns of equilibration and 4 ns of production. You may change the molecules 
to perform the analysis on by editing the SMILES in input, but avoid 
extremely flexible molecules as the mapping function is set to RMSD.

The second script, `openfe_simulation.py`, uses the `openfe` simulation
package from OpenFreeEnergy, and can perform either a AHFE or a RHFE
simulation. The length of the simulations, left as default, is left 
explicit in the main function so that it may be changed for testing.
`openfe` uses `sdf` format files as input.

Please download the `biosimspace` and `openfe` packages.

When you launch the script, directories will be created to setup the 
simulation, then at the end the results will be printed to screen.
You will be warned if the simulation was too short and the results
are not reliable.
If you decide to rerun a simulation, please backup or remove your
previous results. Be aware that the simulation may take a while to run.
