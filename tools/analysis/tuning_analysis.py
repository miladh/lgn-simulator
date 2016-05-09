"""
A collection of functions for analyzing tuning properties of the cells.
"""

import numpy as np

def average_response_vs_attr(sims, attr, rc=[0.5, 0.5]):
    """
    returns the mean response of every cell in
    each simulation with respect to a given simulation
    attribute (attr).

    Parameters
    ----------
    sims : list
        list of Simulation objects
    attr: str
        name of the attribute

    rc: array_like, optional
        normalized receptive field center position.

    Returns
    -------
    dict
        mean response of all cells with with
        respect to attr.

    """
    from collections import defaultdict
    responses = defaultdict(list)
    for sim in sims:
        responses[attr].append(sim.get_attribute(attr))
        for cell_type in sim.cell_types:
            resp = sim.single_cell_temporal_response(cell_type, rc)
            mean_resp = np.mean(resp)
            responses[str(cell_type)].append(mean_resp)

    return responses
