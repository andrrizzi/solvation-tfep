# Alchemical simulations
In this folder, you will find a Python script to automatically perform 
absolute hydration free energy (AHFE) simulations for a given library
of compounds using OpenFE.

The script (`openfe_transform.py`) can print an help string to explain
usage. A library file is necessary, and a force field and a partial charge
method can be specified on input (they default to `openff-2.2.1` and 
`nagl`, but other options are available and can be printed with the 
help of the script).

A directory will be created containing the JSON format transformations 
that can be run using the `openfe` package. Refer to `openfe quickrun` 
documentation for further information.

Results of the simulations for almost (~500/650) all molecules in the 
library are available for `openff-nagl` and `gaff-am1bcc` force field/
partial charge method combinations, and are divided as follows: 
- `simulation_results.json` contains absolute hydration free energy estimates 
and uncertainties obtained from the simulations;
- `diagnostics.json` contains diagnostic information regarding the $\Delta U$
distribution between each force field in each direction (for now, 
`gaff-openff` and `openff-gaff`) at $\lambda=1.0$, such as mean, median, 
standard deviation and effective sample size ratio. These statistics are 
obtained with the functions in `diagnostic_utils.py`.
