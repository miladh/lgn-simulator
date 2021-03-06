import numpy as np
import h5py
import argparse as ap

class Neuron:
    """
    Class for neuron populations
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
        if(self.__resp is None):
            self.__resp = np.array(self.group["response"]["spatio_temporal"])

        return self.__resp

    def resp_ft(self):
        if(self.__resp_ft is None):
            if(isinstance(self.group["response"]["fourier_transform"], h5py.Dataset)):
                self.__resp_ft = np.array(self.group["response"]["fourier_transform"])
            else:
                self.__resp_ft = np.array(
                np.array(self.group["response"]["fourier_transform"]["real"]), dtype=complex)
                self.__resp_ft.imag =np.array(np.array(self.group["response"]["fourier_transform"]["complex"]))

        return self.__resp_ft


    def irf(self):
        if(self.__irf is None):
            self.__irf = np.array(self.group["impulse_response"]["spatio_temporal"])

        return self.__irf

    def irf_ft(self):
        if(self.__irf_ft is None):
            if(isinstance(self.group["impulse_response"]["fourier_transform"], h5py.Dataset)):
                self.__irf_ft = np.array(self.group["impulse_response"]["fourier_transform"])
            else:
                self.__irf_ft = np.array(
                       np.array(self.group["impulse_response"]["fourier_transform"]["real"]),
                       dtype=complex)
                self.__irf_ft.imag = np.array(
                       np.array(self.group["impulse_response"]["fourier_transform"]["complex"]))
        return self.__irf_ft


    def t_resp(self, rc=[0., 0.]):
        return self.resp()[:,  int(rc[1]), int(rc[0])]

    def t_irf(self, rc=[0., 0.]):
        return self.irf()[:,  int(rc[1]), int(rc[0])]

    def w_resp(self, k=[0., 0.]):
        return self.resp_ft()[:,  int(k[1]), int(k[0])]

    def w_irf(self, k=[0., 0.]):
        return self.irf_ft()[:,  int(k[1]), int(k[0])]

    def k_resp(self, w=0, k=0, axis=1):
        if(axis==0):
            return self.resp_ft()[int(w), :, int(k)]
        elif(axis==1):
            return self.resp_ft()[int(w), int(k), :]
        else:
            raise IndexError("axis="+str(axis)+" > 1")

    def k_irf(self, w=0, k=0, axis=1):
        if(axis==0):
            return self.irf_ft()[int(w), :, int(k)]
        elif(axis==1):
            return self.irf_ft()[int(w), int(k), :]
        else:
            raise IndexError("axis="+str(axis)+" > 1")

    def st_resp(self, rx=0.0):
        return self.resp()[:, :, int(rx)]

    def st_irf(self, rx=0.0):
        return self.irf()[:, :, int(rx)]

    def kw_resp(self, k=0, axis=1):
        if(axis==0):
            return self.resp_ft()[:, :, int(k)]
        elif(axis==1):
            return self.resp_ft()[:, int(k), :]
        else:
            raise IndexError("axis="+str(axis)+" > 1")

    def kw_irf(self, k=0, axis=1):
        if(axis==0):
            return self.irf_ft()[:, :, int(k)]
        elif(axis==1):
            return self.irf_ft()[:, int(k), :]
        else:
            raise IndexError("axis="+str(axis)+" > 1")
