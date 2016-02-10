import numpy as np
import h5py

class Integrator:
    """
    Class for integrator properties

    """

    def __init__(self, integrator_group):
        for attr in integrator_group.attrs.keys():
            setattr(self, attr, integrator_group.attrs[attr])

        grid_vectors = integrator_group.keys()
        for grid in grid_vectors:
            values = integrator_group[str(grid)]
            data = np.array(values)
            setattr(self, grid, data)
