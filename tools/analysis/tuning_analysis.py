"""
A collection of functions for analyzing tuning properties of the cells.
"""

import numpy as np

def spatial_summation(sims, rc=[0.5, 0.5]):
    """
    calculates the average response of a cell
    in response to a patch stimulus.

    Parameters
    ----------
    sim : Simulation
        Simulation object
    rc: array_like, optional
        normalized receptive field center position.

    Returns
    -------
    dict
        mean response of all cells in one simulation and
        patch stimulus diameter.

    """
    from collections import defaultdict
    responses = defaultdict(list)
    for sim in sims:
        responses["d"].append(getattr(getattr(sim, "stimulus"), "mask_size"))
        for cell_type in sim.cell_types:
            resp = sim.single_cell_temporal_response(cell_type, rc)
            mean_resp = np.mean(resp)
            responses[str(cell_type)].append(mean_resp)

    return responses



def average_response_vs_weight(sims, key, rc=[0.5, 0.5]):
    """
    calculates the average response of a cell
    in response to a patch stimulus.

    Parameters
    ----------
    sim : Simulation
        Simulation object
    rc: array_like, optional
        normalized receptive field center position.

    Returns
    -------
    dict
        mean response of all cells in one simulation and
        patch stimulus diameter.

    """
    from collections import defaultdict
    responses = defaultdict(list)
    for sim in sims:
        responses[key].append(sim.get_attribute(key))
        for cell_type in sim.cell_types:
            resp = sim.single_cell_temporal_response(cell_type, rc)
            mean_resp = np.mean(resp)
            responses[str(cell_type)].append(mean_resp)

    return responses
