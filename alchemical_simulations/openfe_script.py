def rhfe_simulation(input_file: str, durations) -> None:
    import os, pathlib
    from rdkit import Chem
    from openfe import SmallMoleculeComponent, SolventComponent, lomap_scorers
    from openfe.setup import LomapAtomMapper
    from openfe.setup.ligand_network_planning import generate_minimal_spanning_network
    from openfe.setup.alchemical_network_planner import RHFEAlchemicalNetworkPlanner
    from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol

    # read sdf files in for ligands
    ligands_sdf = Chem.SDMolSupplier(input_file, removeHs = False)
    ligand_mols = [ SmallMoleculeComponent(sdf) for sdf in ligands_sdf ]

    # solvation
    solvent = SolventComponent(positive_ion = 'K', negative_ion = 'Cl', 
                               neutralize = True)

    # protocol - set shorter simulation times
    rhfe_settings = RelativeHybridTopologyProtocol.default_settings()
    rhfe_settings.simulation_settings.equilibration_length, rhfe_settings.simulation_settings.production_length = durations
    rhfe_protocol = RelativeHybridTopologyProtocol(settings = rhfe_settings)

    # create alchemical network of transformations based on mst
    alchem_planner = RHFEAlchemicalNetworkPlanner(
        name = "project",
        mappers = [ LomapAtomMapper() ],
        mapping_scorer = lomap_scorers.default_lomap_score,
        ligand_network_planner = generate_minimal_spanning_network,
        protocol = rhfe_protocol,
    )

    alchemical_network = alchem_planner(
        ligands = ligand_mols,
        solvent = solvent,
    )

    # use openfe CLI to run simulations and gather results
    trans = [t for t in alchemical_network.edges]

    for t in trans :    
        transformation_dir = pathlib.Path("rhfe_json")
        transformation_dir.mkdir(exist_ok=True)

        t.dump(transformation_dir / f"{t.name}.json")

    os.system("cd rhfe_json")
    simulations = [s for s in os.listdir(".") if s.startswith("project")]
    for filename in simulations:
        os.system(f"openfe quickrun {filename}")

    os.system("openfe gather . --report ddg -o results.tsv")

def ahfe_simulation(input_file: str, durations) -> None:
    from rdkit import Chem 
    import os, pathlib
    from openfe import SmallMoleculeComponent, SolventComponent, ChemicalSystem, Transformation
    from openfe.protocols.openmm_afe import AbsoluteSolvationProtocol
    
    # read sdf files in for ligands
    ligands_sdf = Chem.SDMolSupplier(input_file, removeHs = False)
    ligand_mols = [ SmallMoleculeComponent(sdf) for sdf in ligands_sdf ]
    
    # solvation
    solvent = SolventComponent(positive_ion = 'K', negative_ion = 'Cl', 
                               neutralize = True)

    # chemical systems
    systemA = ChemicalSystem({ 'ligand': ligand_mols[0],
                                      'solvent': solvent}, name = ligand_mols[0].name)
    systemB = ChemicalSystem({'solvent': solvent})

    # protocol - set shorter simulations times
    ahfe_settings = AbsoluteSolvationProtocol.default_settings()
    ahfe_settings.solvent_simulation_settings.equilibration_length, ahfe_settings.solvent_simulation_settings.production_length = durations
    ahfe_settings.vacuum_simulation_settings.equilibration_length, ahfe_settings.vacuum_simulation_settings.production_length = durations
    ahfe_protocol = AbsoluteSolvationProtocol(settings = ahfe_settings)

    transformation = Transformation(
        stateA = systemA,
        stateB = systemB,
        mapping = None,
        protocol = ahfe_protocol,
        name = f"{systemA.name}"
    )

    # use openfe CLI to run simulations and gather results
    transformation_dir = pathlib.Path("ahfe_json")
    transformation_dir.mkdir(exist_ok=True)

    transformation.dump(transformation_dir / f"{transformation.name}.json")
    
    os.system("cd ahfe_json")
    simulations = [s for s in os.listdir(".") if s.startswith("project")]
    for filename in simulations:
        os.system(f"openfe quickrun {filename}")

    os.system("openfe gather . --report dg -o results.tsv")

if __name__ == '__main__':
    from openff.units import unit

    durations = { "equilibration": 10, "production": 10 }
    durations = (durations["equilibration"] * unit.picosecond, durations["production"] * unit.picosecond)

    # comment/uncomment based on what simulation you want to run
    ahfe_simulation("input_molecules.sdf", durations)
    # rhfe_simulation("input_molecules.sdf", durations)
