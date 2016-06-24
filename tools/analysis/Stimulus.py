import numpy as np
import h5py

class Stimulus:
    """
    # Class for stimulus

    """
    def __init__(self, h5_group):
        self.group = h5_group
        self.__st = None
        self.__ft = None
        for attr in h5_group.attrs.keys():
            setattr(self, attr, h5_group.attrs[attr])

    def spatio_temporal(self):
        if(self.__st is None):
            self.__st = np.array(self.group["spatio_temporal"])

        return self.__st

    def fourier_transform(self):
        if(self.__ft is None):
            self.__ft = np.array(
                    np.array(self.group["fourier_transform"]["real"]), dtype=complex)
            self.__ft.imag =np.array(np.array(self.group["fourier_transform"]["complex"]))
        return self.__ft
