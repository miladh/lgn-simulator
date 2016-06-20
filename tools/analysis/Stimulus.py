import numpy as np
import h5py

class Stimulus:
    """
    # Class for stimulus used in a single simulation

    """

    def __init__(self, h5_group):
        self.group = h5_group
        for attr in h5_group.attrs.keys():
            setattr(self, attr, h5_group.attrs[attr])

    def __getitem__(self, key):
        return self[key]


    def read_property(self,  space=None):
        if hasattr(self, str(space)):
            return 1

        if(space==None):
            spaces = self.group.keys()
            for sp in spaces:
                values = self.group[str(sp)]
                data = np.array(values)
                setattr(self, sp, data)
        else:
            values = self.group[space]
            data = np.array(values)
            setattr(self, space, data)
