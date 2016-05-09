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


def spatial_summation_curve(sims, param, value):
        """
        returns the area summation vector for
        a list of simulations, all with the same
        setting parameter param=value.

        Parameters
        ----------
        sims : list
            list of response dictionaries.
            Dictionary keys are average responses
            for every cell, patch diameter, and
            any additional parameter.
        param: str
            parameter
        key: float
            value of the parameter key

        Returns
        -------
        dict
            mean response of all cells and
            patch stimulus diameter for a list of
            simulations., all with the same setting
            parameter param=value.

        See Also
        --------
        spatial_summation : returns mean response
        of all cells in one simulation and patch
        stimulus diameter.

        """
        data = {param: value}

        for key in sims[0]:
            if key != param:
                data[key] = []

        for (i, sim) in enumerate(sims):
            if sim[param] == value:
                for key in sim:
                    if key != param:
                        data[key].append(sim[key])

        return data
