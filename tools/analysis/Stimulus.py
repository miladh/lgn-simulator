import numpy as np
import h5py

class Stimulus:
    """
    # Class for stimulus used in a single simulation

    """

    def __init__(self, stim_group):
        for attr in stim_group.attrs.keys():
            setattr(self, attr, stim_group.attrs[attr])

        spaces = stim_group.keys()
        for space in spaces:
            values = stim_group[str(space)]
            data = np.array(values)
            setattr(self, space, data)

    def __getitem__(self, key):
        return self[key]
