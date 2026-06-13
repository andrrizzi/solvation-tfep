import torch
import math
import json
import os
import openmm
import openmm.unit as unit
import numpy as np
from openmm.app import PDBFile, Modeller, ForceField
from netCDF4 import Dataset
import pint

UNITS = pint.UnitRegistry()
TEMPERATURE = 298.15 * UNITS.kelvin  # Temperature of the system.
RT = (TEMPERATURE * UNITS.molar_gas_constant).to(UNITS.kJ / UNITS.mol).magnitude


def get_reference_potentials(input_file: str) -> list:
    """
    Extract reference potential energies from netCDF4 for lambda state 1.0.
    ---
    parameters
    input_file: name of input file
    ---
    returns
    list of potentials for lambda state 1.0
    """
    simulation_data = Dataset(input_file)

    potentials, states = (
        simulation_data.variables["energies"][:],
        simulation_data.variables["states"][:],
    )

    state_index = [int(np.where(x == 13)[0][0]) for x in states]

    lambda1_potentials = []
    for timestep in range(0, len(potentials), 40):
        lambda1_potentials.append(potentials[timestep][state_index[timestep]][13])

    return lambda1_potentials


def get_positions(input_file: str) -> list:
    """
    Extract small molecules positions from netCDF4 for lambda state 1.0.
    ---
    parameters
    input_file: name of input file
    ---
    returns
    list of positions (x,y,z) for lambda state 1.0
    """
    simulation_data = Dataset(input_file)

    positions, states = (
        simulation_data.variables["positions"][:],
        simulation_data.variables["states"][:],
    )

    state_index = [int(np.where(x == 13)[0][0]) for x in states]

    lambda1_positions = []
    for ts in range(0, len(positions), 40):
        lambda1_positions.append(positions[ts][state_index[ts]][:-4])

    return lambda1_positions


def get_target_system(simulation_path: str, forcefield: str) -> openmm.System:
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
    solvent_setup_directories = [
        x
        for x in os.listdir(simulation_path)
        if x.startswith("shared_AHFESolventSetupUnit")
    ]
    vacuum_setup_directories = [
        y
        for y in os.listdir(simulation_path)
        if y.startswith("shared_AHFEVacuumSetupUnit")
    ]

    pdb = PDBFile(f"{simulation_path}/{vacuum_setup_directories[0]}/hybrid_system.pdb")

    with open(
        f"{simulation_path}/{solvent_setup_directories[0]}/db.json",
        "r",
        encoding="utf-8",
    ) as c:
        ff_cache = json.load(c)

    with open(f"{simulation_path}/{forcefield}.xml", "w", encoding="utf-8") as ffxml:
        ffxml.write(ff_cache[forcefield]["1"]["ffxml"])

    ff = ForceField(f"{simulation_path}/{forcefield}.xml")

    modeller = Modeller(pdb.topology, pdb.positions)

    target_system = ff.createSystem(
        modeller.topology,
    )

    return target_system


def reevaluate_in_target_system(positions: list, target_system: openmm.System) -> list:
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
    for _, pos in enumerate(positions):
        context.setPositions(pos)
        state = context.getState(getEnergy=True)
        potential_in_kt = (
            state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole) / RT
        )
        target_energies.append(potential_in_kt)

    return target_energies


def get_statistics(work: torch.Tensor) -> dict:
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

        weights = [math.e ** (-u) for u in delta_u]
        ess = (sum(weights)) ** 2 / sum(w**2 for w in weights)

        return ess / len(weights)

    statistics = {
        "mean": float(torch.mean(work)),
        "median": float(torch.median(work)),
        "standard_deviation": float(torch.std(work)),
        "effective_sample_size_ratio": float(essr(work)),
        "unit": "kT",
    }

    return statistics


def get_work(
    mobley_id: str, reference_forcefield: str, target_forcefield: str
) -> torch.Tensor:
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
    solvent_sim_unit = sorted(
        [
            z
            for z in os.listdir(f"{reference_forcefield}/runs/{mobley_id}_run/")
            if z.startswith("shared_AHFESolventSimUnit")
        ]
    )
    positions = (
        get_positions(
            f"{reference_forcefield}/runs/{mobley_id}_run/{solvent_sim_unit[0]}/solvent.nc"
        )
        + get_positions(
            f"{reference_forcefield}/runs/{mobley_id}_run/{solvent_sim_unit[1]}/solvent.nc"
        )
        + get_positions(
            f"{reference_forcefield}/runs/{mobley_id}_run/{solvent_sim_unit[2]}/solvent.nc"
        )
    )

    system = get_target_system(
        f"{target_forcefield}/runs/{mobley_id}_run",
        f"{target_forcefield}",
    )

    target_energies = reevaluate_in_target_system(positions, system)
    reference_energies = (
        get_reference_potentials(
            f"{reference_forcefield}/runs/{mobley_id}_run/{solvent_sim_unit[0]}/solvent.nc"
        )
        + get_reference_potentials(
            f"{reference_forcefield}/runs/{mobley_id}_run/{solvent_sim_unit[1]}/solvent.nc"
        )
        + get_reference_potentials(
            f"{reference_forcefield}/runs/{mobley_id}_run/{solvent_sim_unit[2]}/solvent.nc"
        )
    )

    return torch.tensor(np.array(target_energies) - np.array(reference_energies))


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

    reference_forcefield, target_forcefield = "openff-2.2.1", "gaff-2.2.20"
    for mid in mobley_ids:
        work = get_work(mid, reference_forcefield, target_forcefield)
        diagnostics_dict[mid][f"{reference_forcefield}-{target_forcefield}"] = (
            get_statistics(work)
        )

    reference_forcefield, target_forcefield = "gaff-2.2.20", "openff-2.2.1"
    for mid in mobley_ids:
        work = get_work(mid, reference_forcefield, target_forcefield)
        diagnostics_dict[mid][f"{reference_forcefield}-{target_forcefield}"] = (
            get_statistics(work)
        )

    with open("diagnostics.json", "w", encoding="utf-8") as o:
        json.dump(diagnostics_dict, o)
