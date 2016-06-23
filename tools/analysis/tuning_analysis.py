"""
A collection of functions for analyzing tuning properties of the cells.
"""

import numpy as np

def average_response_vs_attr(sims, attr,
                            indices=None,
                            rc=[0.5, 0.5],):
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
            if(indices!=None):
                resp = sim.single_cell_temporal_response(cell_type, rc)[[indices]]
            else:
                resp = sim.single_cell_temporal_response(cell_type, rc)
                
            resp =  np.absolute(resp).max()
            responses[str(cell_type)].append(resp)

    return responses


def average_response_vs_attrA_vs_attrB(sims, cell_type, attrA, attrB, rc=[0.5, 0.5]):
    """
    returns the mean response with respect to
    attributes attrA and attrB.

    Parameters
    ----------
    sims : list
        list of Simulation objects
    cell_type: str
        cell type
    attrA: str
        name of the attribute A
    attrB: str
        name of the attribute B
    rc: array_like, optional
        normalized receptive field center position.

    Returns
    -------
    dict
        array with attrA values,
        array with attrB values,
        2d response array.

    """
    import data_extractor as de
    attrA_vec = de.extract_unique_simulation_attrs(sims, attrA)
    attrB_vec = de.extract_unique_simulation_attrs(sims, attrB)

    response = np.zeros([len(attrA_vec), len(attrB_vec)])

    print attrA, attrB
    for i, a in enumerate(attrA_vec):
        sims_ext = de.simulation_extractor(sims, attrA, a)
        data = average_response_vs_attr(sims_ext, attrB)

        sorted_indices = np.argsort(data[attrB])
        response[i,:] = np.array(data[cell_type])[sorted_indices]

    return attrA_vec, attrB_vec, response
