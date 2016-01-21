import numpy as np

class Simulation:
    """
    Class for a single eDOG simulation

    """
    def __init__(self, h5_file):
        self.num_steps = h5_file.attrs["nSteps"]
        self.num_points = h5_file.attrs["nPoints"]
        self.dt = h5_file.attrs["dt"]
        self.ds = h5_file.attrs["ds"]
        self.time_vec = np.arange(0, self.num_steps*self.dt, self.dt )

        self.cellTypes = []

        ########################## Read file ###################################
        for item in h5_file.keys():
            if(item == "stimulus"):
                stim_group = h5_file.get("/" + str(item))
                for attr in stim_group.attrs.keys():
                    setattr(self, attr, stim_group.attrs[attr])

                spaces = stim_group.keys()
                data = {}
                for space in spaces:
                    values = h5_file.get("/" + str(item) + "/" + str(space))
                    data.update({str(space):  np.array(values)})
                    setattr(self, item, data)

            else:
                self.cellTypes.append(item)
                cellAttr = h5_file.get("/" + str(item)).keys()
                data = {}
                for attr in cellAttr:
                    spaces = h5_file.get("/"+str(item)+"/" +str(attr)).keys()
                    data[attr] = {}
                    for space in spaces[0:]:
                        values = np.array(h5_file.get("/"+str(item)+"/"+str(attr)+"/"+str(space)))
                        data[attr].update({str(space): np.array(values)})
                setattr(self, item, data)
        ########################################################################
        self.numCellTypes  = len(self.cellTypes)
        self.normalize()

    def normalize(self, cellType=None):
        if cellType == None:
            for cell in self.cellTypes:
                temp =  getattr(self, cell)
                for attr in temp.keys():
                    for i in range(temp[attr]["spatioTemporal"].shape[0]):
                        if not abs(temp[attr]["spatioTemporal"][i,:,:]).max()==0:
                            temp[attr]["spatioTemporal"][i,:,:]/= \
                            abs(temp[attr]["spatioTemporal"][i,:,:]).max()
                        # temp[attr]["spatioTemporal"][i,:,:]+=1
                    setattr(self, cell,temp)
        else:
            temp =  getattr(self, cellType)
            for attr in temp.keys():
                for i in range(temp[attr]["spatioTemporal"].shape[0]):
                    if not abs(temp[attr]["spatioTemporal"][i,:,:]).max()==0:
                        temp[attr]["spatioTemporal"][i,:,:]/= \
                        abs(temp[attr]["spatioTemporal"][i,:,:]).max()
                        # temp[attr]["spatioTemporal"][i,:,:] +=1
                setattr(self, cellType,temp)


    def singleCellTemporalResponse(self, cellType, idx=0 , idy=0):
        response = getattr(self, cellType)["response"]["spatioTemporal"][:,idx, idy]
        return response

    def singleCellFreqResponse(self, cellType, idx=0 , idy=0):
        FreqResponse = getattr(self, cellType)["response"]["fourierTransform"][:,idx, idy]
        return FreqResponse

    def spikeTrain(self, cellType, idx=0 , idy=0, num_trails=2):
        response = self.singleCellTemporalResponse(cellType, idx , idy)
        spike_times = [[] for i in range(num_trails)]

        for k in range(num_trails):
            for i in range(self.num_steps):
                r = np.random.uniform(0,1)
                if(response[i]+1 * self.dt > r):
                    spike_times[k].append(i*self.dt)
        return spike_times

if __name__ == "__main__":
    import h5py
    from glob import glob
    from pylab import*

    outputFilePath = "/home/milad/Dropbox/projects/edog/extendedDOG/DATA/spatialSummation/tmp/*.h5"
    outputFile = glob(outputFilePath)[0]
    f = h5py.File(outputFile, "r")
    sim = Simulation(f)
    print sim.ganglion["response"].keys()
    print sim.maskSize
