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
            if(isinstance(self.group["fourier_transform"], h5py.Dataset)):
                self.__ft = np.array(self.group["fourier_transform"])
            else:
                self.__ft = np.array(
                        np.array(self.group["fourier_transform"]["real"]), dtype=complex)
                self.__ft.imag =np.array(np.array(self.group["fourier_transform"]["complex"]))
        return self.__ft


    def t_domain(self, rc=[0., 0.]):
        return self.spatio_temporal()[:, int(rc[1]), int(rc[0])]

    def w_domain(self, k=[0., 0.]):
        return self.fourier_transform()[:, int(k[1]), int(k[0])]

    def k_domain(self, w=0, k=0, axis=1):
        if(axis==0):
            return self.fourier_transform()[int(w), :, int(k)]
        elif(axis==1):
            return self.fourier_transform()[int(w), int(k), :]
        else:
            raise IndexError("axis="+str(axis)+" > 1")
