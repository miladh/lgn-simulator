import h5py
from glob import glob
import numpy as np

class Simulation:
    """
    Class for a single eDOG simulation

    """
    def __init__(self, h5_file):
        self.num_cell_types =  len(h5_file.keys())

        for item in h5_file.keys():
            if(item== "stimuli"):
                self.num_cell_types -=1
                spatioTemporal = np.array(h5_file.get("/" + str(item)+ "/real"))
                fourierTransform = np.array(h5_file.get("/" + str(item)+ "/complex"))
                data = {"spatioTemporal": spatioTemporal, "fourierTransform": fourierTransform}
                setattr(self, item, data)

            else:
                cellAttr = h5_file.get("/" + str(item)).keys()
                data = {}
                for attr in cellAttr:
                    spatioTemporal = np.array(h5_file.get("/" + str(item)+ "/" +str(attr) + "/real"))
                    fourierTransform =  np.array(h5_file.get("/" + str(item)+ "/" +str(attr) + "/complex"))
                    data[attr] = {"spatioTemporal": spatioTemporal, "fourierTransform": fourierTransform}

                setattr(self, item, data)


if __name__ == "__main__":
    outputFilePath = "/home/milad/Dropbox/projects/edog/extendedDOG/eDOG/DATA/*.h5"
    outputFile = glob(outputFilePath)[0]
    f = h5py.File(outputFile, "r")
    sim = Simulation(f)
    print sim.ganglion.keys()
    print sim.stimuli.keys()
