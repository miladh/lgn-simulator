import numpy as np
import h5py
import argparse as ap

class Neuron:
    """
    Class for neurons
    """

    def __init__(self, h5_group):
        self.__resp = None
        self.__resp_ft = None
        self.__irf = None
        self.__irf_ft = None
        self.group = h5_group
        for attr in h5_group.attrs.keys():
            setattr(self, attr, h5_group.attrs[attr])

    def resp(self):
        if(self.__resp==None):
            self.__resp = np.array(self.group["response"]["spatio_temporal"])

        return self.__resp

    def resp_ft(self):
        if(self.__resp_ft==None):
            self.__resp_ft = np.array(
                    np.array(self.group["response"]["fourier_transform"]["real"]), dtype=complex)
            self.__resp_ft.imag =np.array(np.array(self.group["response"]["fourier_transform"]["complex"]))
        return self.__resp_ft


    def irf(self):
        if(self.__irf==None):
            self.__irf = np.array(self.group["impulse_response"]["spatio_temporal"])

        return self.__irf

    def irf_ft(self):
        if(self.__irf_ft==None):
            self.__irf_ft = np.array(
                   np.array(self.group["impulse_response"]["fourier_transform"]["real"]),dtype=complex)
            self.__irf_ft.imag = np.array(
                   np.array(self.group["impulse_response"]["fourier_transform"]["complex"]))
        return self.__irf_ft
