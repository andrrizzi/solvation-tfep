# Dataset evaluation
Molecules were sampled from the FreeSolv dataset and binned by 
flexibility. The CalcNumRotatableBonds function of the RDKit package 
in the rdMolDescriptors module was employed in strict mode to count 
all rotatable bonds in the molecules. Furthermore, molecules 
containing aliphatic rings were identified using the 
CalcNumAliphaticRings function of the same module, and these entities
were then visually inspected to count "rotatable" bonds in these 
rings: this includes all bonds in aliphatic rings which do not 
involve carbon atoms in sp^2^ hybridisation state (meaning carbon 
atoms that participate in double bonds), because these atoms are 
constrained on the same plane and do not undergo meaningful rotation.

In order to store this information regarding molecule flexibility, a 
JSON file was generated with four keys:
- `molecule_ID` (string), which indicates the molecule ID attributed 
in the FreeSolv dataset;
- `only_aromatic_ring` (bool), which indicates whether the molecule 
contains aliphatic rings;
- `num_rotatable_bonds` (int), the amount of classically defined 
rotatable bonds in the molecule as obtained with the RDKit package;
- `num_aliphatic_bonds` (int), the amount of bonds in aliphatic rings 
counted by visual inspection, as defined above.

The directory also includes a manually annotated JSON with the 
numbers of cyclic aliphatic bonds, a library in SDF format and a 
script in Python to use the two former files to count rotatable bonds
in the entire dataset.
