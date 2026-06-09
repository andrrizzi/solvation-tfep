import torch
import numpy as np
from netCDF4 import Dataset
import math


def get_potentials(input_file: str) -> list:
    """
    Extract potential energies from netCDF4 for lambda state 1.0.
    ---
    parameters
    input_file: name of input file
    ---
    returns
    list of potentials for lambda state 1.0
    """
    data = Dataset(input_file)

    potentials, states = (
        data.variables["energies"][:],
        data.variables["states"][:],
    )

    state_index = [int(np.where(x == 13)[0][0]) for x in states]

    lambda1_potentials = [
        float(z[state_index[y]][13]) for y, z in enumerate(potentials)
    ]

    return lambda1_potentials


def get_statistics(ff1: list, ff2: list) -> dict:
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

    ff1, ff2 = np.array(ff1), np.array(ff2)
    work = torch.tensor(ff2 - ff1)
    statistics = {
        "mean": float(torch.mean(work)),
        "median": float(torch.median(work)),
        "standard_deviation": float(torch.std(work)),
        "effective_sample_size_ratio": float(essr(work)),
        "unit": "kT",
    }

    return statistics
