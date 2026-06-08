import pathlib
import argparse

from rdkit import Chem
from openfe import (
    SmallMoleculeComponent,
    SolventComponent,
    ChemicalSystem,
    Transformation,
)
from openfe.protocols.openmm_afe import AbsoluteSolvationProtocol

def transform(input_file: str, ff: str, pcm: str) -> None:
    """ 
    Read input file, create the system and print JSON format transformations 
    to run absolute hydration free energy (AHFE) simulations for a library of compounds.
    ---
    parameters
    input_file: SDF format file containing the library of compounds.
    ff: force field for the AHFE simulations
    pcm: partial charge method for the AHFE simulations
    ---
    """
    # read sdf files in for ligands
    ligands_sdf = Chem.SDMolSupplier(input_file, removeHs=False)
    ligand_mols = [SmallMoleculeComponent(sdf) for sdf in ligands_sdf]

    # solvation
    solvent = SolventComponent(positive_ion="K", negative_ion="Cl", neutralize=True)

    # protocol - set partial charge model and small molecule forcefield
    ahfe_settings = AbsoluteSolvationProtocol.default_settings()

    ahfe_settings.partial_charge_settings.partial_charge_method = pcm
    ahfe_settings.vacuum_forcefield_settings.small_molecule_forcefield = ff
    ahfe_settings.solvent_forcefield_settings.small_molecule_forcefield = ff

    # protocol - gpu device index for parallel computing
    ahfe_settings.solvent_engine_settings.gpu_device_index = [0]
    ahfe_settings.vacuum_engine_settings.gpu_device_index = [0]

    ahfe_protocol = AbsoluteSolvationProtocol(settings=ahfe_settings)

    transformation_dir = pathlib.Path("transformations")
    transformation_dir.mkdir(exist_ok=True)

    for ligand in ligand_mols:
        # chemical systems
        system_a = ChemicalSystem(
            {"ligand": ligand, "solvent": solvent}, name=ligand.name
        )

        system_b = ChemicalSystem({"solvent": solvent})

        transformation = Transformation(
            stateA=system_a,
            stateB=system_b,
            mapping=None,
            protocol=ahfe_protocol,
            name=f"{system_a.name}",
        )

        # use openfe CLI to run simulations and gather results
        transformation.dump(transformation_dir / f"{transformation.name}.json")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate transformation JSON files to run absolute hydration"
                    + "free energy simulations for a library of compounds."
    )
    parser.add_argument("compoundLibrary", help="SDF library of compounds.", type=str)
    parser.add_argument(
        "-ff",
        "--forceField",
        type=str,
        choices=["openff-2.2.1", "gaff-2.2.20", "espaloma-0.4.0"],
        default="openff-2.2.1",
        help="Force field for the simulations.",
    )
    parser.add_argument(
        "-pcm",
        "--partialChargeMethod",
        type=str,
        choices=["nagl", "am1bcc", "espaloma"],
        default="nagl",
        help="Partial charge method for the simulations.",
    )
    args = parser.parse_args()

    transform(args.compoundLibrary, args.forceField, args.partialChargeMethod)
