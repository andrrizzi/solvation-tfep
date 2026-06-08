# Alchemical simulations
In this folder, you will find a Python script to automatically perform 
absolute hydration free energy (AHFE) simulations for a given library
of compounds using OpenFE.

The script (`openfe_transform.py`) can print an help string to explain
usage. A library file is necessary, and a force field and a partial charge
method can be specified (they default to `openff-2.2.1` and `nagl`, but
other options are available).

A directory will be created containing the JSON format transformations 
that can be run using the `openfe` package. Refer to `openfe` 
documentation for further information.

Results of the simulations for almost (500/650) all molecules in the 
library are available for openff-nagl and gaff-am1bcc force field/partial
charge method combinations. The results, indexed by mobley ID, consist of:
- a results section, with the absolute hydration free energy estimate and 
uncertainty for each force field and partial charge method combination 
(`openff-nagl`, `gaff-am1bcc`);
- an analysis section, with values describing the distributions of potentials 
of the simulations obtained with FEP (see andrrizzi/tfep) (`mean`, `median`,
`standard_deviation`, `confidence_interval`, `effective_sample_size`) for 
each force field and partial charge method in both directions (`openff-gaff`,
`gaff-openff`).

File snippet below:
```
{
    ...
    "mobley_1723043": {
        "results": {
            "gaff-am1bcc": {
                "estimate": {
                    "magnitude": 1.4730056952143862,
                    "unit": "kilocalorie_per_mole"
                },
                "uncertainty": {
                    "magnitude": 0.06999298380487426,
                    "unit": "kilocalorie_per_mole"
                }
            },
            "openff-nagl": {
                "estimate": {
                    "magnitude": 2.456256739808371,
                    "unit": "kilocalorie_per_mole"
                },
                "uncertainty": {
                    "magnitude": 0.012369145070170616,
                    "unit": "kilocalorie_per_mole"
                }
            }
        },
        "analysis": {
            "fep": {
                "gaff-openff": {
                    "mean": -275.28871637805076,
                    "median": -275.73644854215354,
                    "standard_deviation": 2.943703851337015,
                    "confidence_interval": {
                        "low": -279.9596987866364,
                        "high": -268.2420609650112
                    },
                    "effective_sample_size": 1.0000023594773202
                },
                "openff-gaff": {
                    "mean": -354.2557583160784,
                    "median": -354.69171148266446,
                    "standard_deviation": 3.079063592947018,
                    "confidence_interval": {
                        "low": -358.8039043236034,
                        "high": -346.92646779010386
                    },
                    "effective_sample_size": NaN
                }
            }
        }
    },
    ...
}
```

