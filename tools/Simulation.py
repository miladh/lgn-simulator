import numpy as np

class Simulation:
    """
    Class for a single eDOG simulation

    """
    def __init__(self, h5_file):
        self.numCellTypes =  len(h5_file.keys())
        self.nSpatialPoints = 256 #READ FROM FILE!!!!!
        self.nTemporalSteps = 16  #READ FROM FILE!!!!!

        for item in h5_file.keys():
            if(item== "stimuli"):
                self.numCellTypes -=1
                spatioTemporal = np.array(h5_file.get("/" + str(item)+ "/real"))
                fourierTransform = np.array(h5_file.get("/" + str(item)+ "/complex"))
                data = {"spatioTemporal": spatioTemporal, "fourierTransform": fourierTransform}
                setattr(self, item, data)

            else:
                cellAttr = h5_file.get("/" + str(item)).keys()
                data = {}
                for attr in cellAttr:
                    spatioTemporal = np.array(h5_file.get("/" + str(item) + "/"  +str(attr) + "/real"))
                    fourierTransform =  np.array(h5_file.get("/" + str(item) + "/" + str(attr) + "/complex"))
                    data[attr] = {"spatioTemporal": spatioTemporal, "fourierTransform": fourierTransform}

                setattr(self, item, data)


    def singleCellTemporalResponse(self, cellType, idx=0 , idy=0):
        response = getattr(self, cellType)["response"]["spatioTemporal"][:,idx, idy]
        return response

    def singleCellFreqResponse(self, cellType, idx=0 , idy=0):
        FreqResponse = getattr(self, cellType)["response"]["fourierTransform"][:,idx, idy]
        return FreqResponse


if __name__ == "__main__":
    import h5py
    from glob import glob
    from pylab import*
    outputFilePath = "/home/milad/Dropbox/projects/edog/extendedDOG/eDOG/DATA/*.h5"
    outputFile = glob(outputFilePath)[0]
    f = h5py.File(outputFile, "r")
    sim = Simulation(f)
    print sim.ganglion["response"]["spatioTemporal"].shape
    print sim.stimuli.keys()
    plot(sim.singleCellFreqResponse("ganglion"))
    show()
