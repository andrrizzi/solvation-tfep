def rhfe_simulation(input_file: str, durations) -> None:
    import os, pathlib
    from rdkit import Chem
    from openfe import SmallMoleculeComponent, SolventComponent, lomap_scorers
    from openfe.setup import LomapAtomMapper
    from openfe.setup.ligand_network_planning import generate_minimal_spanning_network
    from openfe.setup.alchemical_network_planner import RHFEAlchemicalNetworkPlanner
    from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol

    # read sdf files in for ligands
    ligands_sdf = Chem.SDMolSupplier(input_file, removeHs=False)
    ligand_mols = [SmallMoleculeComponent(sdf) for sdf in ligands_sdf]

    # solvation
    solvent = SolventComponent(positive_ion="K", negative_ion="Cl", neutralize=True)

    # protocol - set simulation times
    rhfe_settings = RelativeHybridTopologyProtocol.default_settings()
    (
        rhfe_settings.simulation_settings.equilibration_length,
        rhfe_settings.simulation_settings.production_length,
    ) = durations
    rhfe_protocol = RelativeHybridTopologyProtocol(settings=rhfe_settings)

    # create alchemical network of transformations based on mst
    alchem_planner = RHFEAlchemicalNetworkPlanner(
        name="project",
        mappers=[LomapAtomMapper()],
        mapping_scorer=lomap_scorers.default_lomap_score,
        ligand_network_planner=generate_minimal_spanning_network,
        protocol=rhfe_protocol,
    )

    alchemical_network = alchem_planner(
        ligands=ligand_mols,
        solvent=solvent,
    )

    # use openfe CLI to run simulations and gather results
    trans = list(alchemical_network.edges)

    for t in trans:
        transformation_dir = pathlib.Path("rhfe_json")
        transformation_dir.mkdir(exist_ok=True)

        t.dump(transformation_dir / f"{t.name}.json")

    # os.system("cd rhfe_json")
    # simulations = [s for s in os.listdir(".") if s.startswith("project")]
    # for filename in simulations:
    #     os.system(f"openfe quickrun {filename}")
    #
    # os.system("openfe gather . --report ddg -o results.tsv")


def ahfe_simulation(input_file: str, durations) -> None:
    """
    Generate an AHFE transformation json for each structures
    in the SDF input file.
    """
    import pathlib

    from rdkit import Chem
    from openfe import (
        SmallMoleculeComponent,
        SolventComponent,
        ChemicalSystem,
        Transformation,
    )
    from openfe.protocols.openmm_afe import AbsoluteSolvationProtocol

    # read sdf files in for ligands
    ligand_sdfs = Chem.SDMolSupplier(input_file, removeHs=False)
    ligand_mols = [SmallMoleculeComponent(sdf) for sdf in ligand_sdfs]

    # solvation
    solvent = SolventComponent(positive_ion="K", negative_ion="Cl", neutralize=True)

    # protocol - set simulations times
    ahfe_settings = AbsoluteSolvationProtocol.default_settings()
    (
        ahfe_settings.solvent_simulation_settings.equilibration_length,
        ahfe_settings.solvent_simulation_settings.production_length,
        ahfe_settings.vacuum_simulation_settings.equilibration_length,
        ahfe_settings.vacuum_simulation_settings.production_length,
    ) = durations

    ahfe_settings.protocol_repeats = 1
    ahfe_protocol = AbsoluteSolvationProtocol(settings=ahfe_settings)

    transformation_dir = pathlib.Path("ahfe_json")
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
    from openff.units import unit

    # rhfe simulation with default settings
    # rhfe_durations = { "equilibration": 1, "production": 5 }
    # rhfe_durations = (
        # rhfe_durations["equilibration"] * unit.nanosecond,
        # rhfe_durations["production"] * unit.nanosecond
    # )
    # rhfe_simulation("input_molecules.sdf", rhfe_durations)

    # ahfe simulation with default settings
    ahfe_durations = {
        "solvent": {"equilibration": 1, "production": 10},
        "vacuum": {"equilibration": 0.5, "production": 2},
    }
    # ahfe_durations = {
    #     "solvent": {"equilibration": 0.01, "production": 10},
    #     "vacuum": {"equilibration": 0.01, "production": 0.01},
    # }

    ahfe_durations = (
        ahfe_durations["solvent"]["equilibration"] * unit.nanosecond,
        ahfe_durations["solvent"]["production"] * unit.nanosecond,
        ahfe_durations["vacuum"]["equilibration"] * unit.nanosecond,
        ahfe_durations["vacuum"]["production"] * unit.nanosecond,
    )
    ahfe_simulation("input_molecules.sdf", ahfe_durations)
