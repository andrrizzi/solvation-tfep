# Alchemical simulations
In this folder, you will find a Python script to automatically perform 
a relative free energy of solvation analysis from methane to iodomethane.

The script uses the `biosimspace` simulation package from OpenBioSim,
and performs a simulation for 40 ps of equilibration and 100 ps of
production. You may change these parameters (equilibration and production 
duration) to your liking by editing the main function in the script.
You may also change the molecules to perform the analysis on, but
avoid extremely flexible molecules as the mapping function is set to 
RMSD.

Please download the `biosimspace` package from the OpenBioSim suite.

When you launch the script, directories will be created to setup the 
simulation, then at the end the results will be printed to screen.
You will be warned if the simulation was too short and the results
are not reliable.
If you decide to rerun a simulation, please backup or remove your
previous results. Be aware that the simulation may take a while to run.
