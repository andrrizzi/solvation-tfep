import json
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

torsionJSON = []
library = [mol for mol in Chem.SDMolSupplier("library.sdf")]

with open("aliphatic.json", "r") as alinput :
    aliJSON = json.load(alinput)

for mol in library :
    molDict = {}

    molDict["molecule_ID"] = mol.GetProp('_Name')
    molDict["only_aromatic_rings"] = True if rdMolDescriptors.CalcNumAliphaticRings(mol) == 0 else False
    molDict["num_rotatable_bonds"] = rdMolDescriptors.CalcNumRotatableBonds(mol) 
    molDict["num_aliphatic_bonds"] = 0
            
    for item in aliJSON :
        if item["molecule_ID"] == molDict["molecule_ID"] :
            molDict["num_aliphatic_bonds"] = int(item["num_torsions"])
            break 

    torsionJSON.append(molDict)

with open("output.json", "w") as output :
    output.write(json.dumps(torsionJSON, indent = 2))
