import math
import json
import openmm
import pint
import openmm.unit as unit
import numpy as np
import pandas as pd
from pathlib import Path
from itertools import permutations
from openmm.app import ForceField, Modeller, PDBFile
from netCDF4 import Dataset

UNITS = pint.UnitRegistry()
TEMPERATURE = 298.15 * UNITS.kelvin
RT = (TEMPERATURE * UNITS.molar_gas_constant).to(UNITS.kJ / UNITS.mol).magnitude


def get_sim_info(input_dir: Path) -> tuple:
    """
    Extract info from netCDF4 concerning potentials, positions and box vectors.
    ---
    parameters
    input_dir (str): directory of the molecule to analyse
    ---
    returns
    lambda1_potentials (list): potentials for simulation at lambda state 1.0
    lambda1_positions (list): positions for simulation at lambda state 1.0
    lambda1_box_vectors(list): box vectors for simulation at lambda state 1.0
    """
    solvent_simulation_files = sorted(
        [
            z
            for z in input_dir.iterdir()
            if z.name.startswith("shared_AHFESolventSimUnit")
        ]
    )

    lambda1_potentials, lambda1_positions, lambda1_box_vectors = [], [], []
    for filename in solvent_simulation_files:
        simulation_data = Dataset(filename / "solvent.nc")
        position_data = Dataset(filename / "solvent_checkpoint.nc")

        # energies[iteration][replica][state], (4001, 14, 14)
        # states[iteration][replica], (4001, 14)
        potentials, states = (
            simulation_data.variables["energies"][:],
            simulation_data.variables["states"][:],
        )

        # positions[iteration][replica][atom][spatial], (4001, 14, 16, 3)
        # box_vectors[iteration][replica][i][j], (4001, 14, 3, 3)
        positions, box_vectors = (
            position_data.variables["positions"][:],
            position_data.variables["box_vectors"][:],
        )

        state_index = [int(np.where(x == 0)[0][0]) for x in states]

        for ts in range(0, len(potentials), 400):
            lambda1_potentials.append(potentials[ts][state_index[ts]][0])

        for ts in range(len(positions)):
            lambda1_positions.append(positions[ts][state_index[ts * 400]])
            lambda1_box_vectors.append(
                box_vectors[ts][state_index[ts * 400]].data.tolist()
            )

    return lambda1_potentials, lambda1_positions, lambda1_box_vectors


def create_target_system(
    simulation_path: Path, target_forcefield: str, number_of_positions: int
) -> openmm.System:
    """
    Generate a system in the target force-field.
    ---
    parameters
    simulation_path (str): directory containing simulation data
    forcefield (str): target force field
    ---
    returns
    target_system (openmm.System): OpenMM system in target force field
    """

    solvent_setup_directories = sorted(
        [
            x
            for x in simulation_path.iterdir()
            if x.name.startswith("shared_AHFESolventSetupUnit")
        ]
    )

    vacuum_setup_directories = sorted(
        [
            y
            for y in simulation_path.iterdir()
            if y.name.startswith("shared_AHFEVacuumSetupUnit")
        ]
    )

    forcefield_input_path = solvent_setup_directories[0] / "db.json"
    forcefield_output_path = simulation_path / f"{target_forcefield}.xml"
    pdb_input_path = vacuum_setup_directories[0] / "hybrid_system.pdb"

    with forcefield_input_path.open(mode="r", encoding="utf-8") as db:
        ff = json.load(db)

    forcefield_output_path.write_text(ff[target_forcefield]["1"]["ffxml"], encoding="utf-8")

    forcefield = ForceField(
        forcefield_output_path,
        "amber/tip3p_standard.xml",
        "amber/tip3p_HFE_multivalent.xml",
    )

    pdb = PDBFile(str(pdb_input_path))
    num_mol_atoms = len(pdb.getPositions())
    solvent_atoms = number_of_positions - num_mol_atoms - 4
    num_solvent_molecules = int((solvent_atoms / 3) + 4)

    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addSolvent(
        forcefield,
        model="tip3p",
        boxShape="cube",
        numAdded=num_solvent_molecules,
        positiveIon="K+",
        negativeIon="Cl-",
        ionicStrength=0.15 * unit.molar,
        neutralize=True,
    )

    target_system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=openmm.app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=openmm.app.HBonds,
        rigidWater=True,
        removeCMMotion=False,
        hydrogenMass=3.0 * unit.amu,
    )

    return target_system


def reevaluate_in_target_system(
    positions: list, box_vectors: list, target_system: openmm.System
) -> list:
    """
    Re-evaluate potential energies of simulation in target force field.
    ---
    parameters
    positions (list): list of small molecule atom positions (x,y,z)
    target_system (openmm.System): OpenMM system in target force field
    ---
    returns
    target_energies (list): potentials evaluated in the target force field.
    """

    integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
    platform = openmm.Platform.getPlatformByName("CPU")
    context = openmm.Context(target_system, integrator, platform)

    target_energies = []
    for ts, pos in enumerate(positions):
        context.setPeriodicBoxVectors(
            box_vectors[ts][0] * unit.nanometer,
            box_vectors[ts][1] * unit.nanometer,
            box_vectors[ts][2] * unit.nanometer,
        )
        context.setPositions(pos)
        state = context.getState(getPositions=True, getEnergy=True)
        potential_in_kt = (
            state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole) / RT
        )
        target_energies.append(potential_in_kt)

    return target_energies


def get_statistics(work: np.array) -> dict:
    """
    Retrieve diagnostics for the simulation results.
    ---
    parameters
    ff1: lambda=1.0 potentials for force field 1
    ff2: lambda=1.0 potentials for force field 2
    ---
    returns
    dict containing relevant diagnostic statistics such as mean,
    median, standard deviation and effective sample size ratio
    of work
    """

    def essr(delta_u: np.array) -> float:
        """
        Effective sample size ratio.
        ---
        parameters
        delta_u: np.array containing the difference between potentials
        of force field 1 and force field 2.
        ---
        returns
        float of effective sample size ratio
        """

        weights = [math.e ** (-u) for u in delta_u if not np.isnan(u)]
        ess = (sum(weights)) ** 2 / sum(w**2 for w in weights)

        return ess / len(weights)

    statistics = {
        "mean": float(np.nanmean(work)),
        "median": float(np.nanmedian(work)),
        "standard_deviation": float(np.nanstd(work)),
        "effective_sample_size_ratio": float(essr(work)),
        "unit": "kT",
    }

    return statistics


def get_work(
    mobley_id: str, reference_forcefield: str, target_forcefield: str
) -> np.array:
    """
    Calculate work between potentials in reference force field and
    potentials re-evaluated in target forcefield.
    ---
    parameters
    mobley_id (str): mobley_id of molecule to analyse
    reference_forcefield (str): original force field of simulation
    target_forcefield (str): target force field to re-evaluate in
    ---
    returns
    work (torch.Tensor): tensor containing the difference in potentials
                         between the force fields
    """
    def remove_outliers(energies: list) -> np.array:
        """
        Remove outlier energies that mess up the statistics using 
        interquartile distance from the median.
        ---
        parameters
        energies (np.array): np.array that contains the energy values
        ---
        returns
        np.array with outliers replaced with NaNs
        """
        df = pd.DataFrame({'energies': energies})
        df_sub = df.loc[:, 'energies']

        iqr = df_sub.quantile(0.75) - df_sub.quantile(0.25)
        lim = np.abs((df_sub - df_sub.median()) / iqr) < 2.22

        df.loc[:, 'energies'] = df_sub.where(lim, np.nan)

        return np.array(df.loc[:, 'energies'])

    mol_ref_run_path = Path(f"{reference_forcefield}/runs/{mobley_id}_run/")
    mol_tar_run_path = Path(f"{target_forcefield}/runs/{mobley_id}_run/")
    reference_energies, positions, box_vectors = get_sim_info(mol_ref_run_path)

    system = create_target_system(
        mol_tar_run_path, target_forcefield, len(positions[0])
    )

    target_energies = reevaluate_in_target_system(positions, box_vectors, system)

    target_energies = remove_outliers(target_energies)
    reference_energies = remove_outliers(reference_energies)

    return target_energies - reference_energies


if __name__ == "__main__":
    with open("simulation_results.json", "r", encoding="utf-8") as s:
        results_dict = json.load(s)

    mobley_ids = []
    for key in results_dict.keys():
        try:
            x = results_dict[key]["ahfe_results"]["openff-nagl"]
            y = results_dict[key]["ahfe_results"]["gaff-am1bcc"]
            mobley_ids.append(key)
        except:
            continue

    diagnostics_dict = {}
    forcefields = ["openff-2.2.1", "gaff-2.2.20"]
    for ref, tar in list(permutations(forcefields, 2)):
        for mid in mobley_ids:
            work = get_work(mid, ref, tar)
            diagnostics_dict.setdefault(mid, {})
            diagnostics_dict[mid][f"{ref}-{tar}"] = (
                get_statistics(work)
            )

    print(diagnostics_dict)

    with open("diagnostics_v2.json", "w", encoding="utf-8") as o:
        json.dump(diagnostics_dict, o, indent=4)
