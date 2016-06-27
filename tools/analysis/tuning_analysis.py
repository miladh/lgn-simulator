"""
A collection of functions for analyzing tuning properties of the cells.
"""

import numpy as np

def resp_vs_attr(sims, attr, resp_type = "t_resp", rc=[0.0, 0.0], indices=None):
    """
    returns the max response of every cell in
    each simulation with respect to a given simulation
    attribute (attr).

    Parameters
    ----------
    sims : list
        list of Simulation objects
    attr: str
        name of the attribute
    resp_type: str
        response type: t_resp, w_resp
    rc: array_like
        neuron position indices
    indices:
        response indices

    Returns
    -------
    dict
        max response of all cells with with
        respect to attr.

    """
    from collections import defaultdict
    resps = defaultdict(list)
    for sim in sims:
        resps[attr].append(sim.get_attribute(attr))
        for cell_type in sim.cell_types:
            neuron =  getattr(sim, cell_type)
            if(indices is not None):
                resp = getattr(neuron, resp_type)(rc)[[indices]]
            else:
                resp = getattr(neuron, resp_type)(rc)

            #resp =  np.absolute(resp).max()
            resp = max(resp, key=np.absolute)
            resps[str(cell_type)].append(resp)

    return resps



def resp_vs_attrA_vs_attrB(sims, attrA, attrB, resp_type = "t_resp", rc=[0.0, 0.0], indices=None):
    """
    returns the max response with respect to
    attributes attrA and attrB.

    Parameters
    ----------
    sims : list
        list of Simulation objects
    attrA: str
        name of the attribute A
    attrB: str
        name of the attribute B
    resp_type: str
        response type: t_resp, w_resp
    rc: array_like
        neuron position indices
    indices:
        response indices

    Returns
    -------
    resp: dict
        2d response array (attrA x attrB).
    attrA_vec: array
        array with attrA values
    attrB_vec:
        array with attrB values
    """
    import data_extractor as de
    attrA_vec = de.extract_unique_simulation_attrs(sims, attrA)
    attrB_vec = de.extract_unique_simulation_attrs(sims, attrB)

    resp={}
    for cell_type in sims[0].cell_types:
        resp[str(cell_type)] = np.zeros([len(attrA_vec), len(attrB_vec)])

    for i, a in enumerate(attrA_vec):
        sims_ext = de.simulation_extractor(sims, attrA, a)
        data = resp_vs_attr(sims_ext, attrB, resp_type, rc, indices)

        sorted_indices = np.argsort(data[attrB])
        for cell_type in sims[0].cell_types:
            resp[str(cell_type)][i,:] = np.array(data[cell_type])[sorted_indices]

    return attrA_vec, attrB_vec, resp
