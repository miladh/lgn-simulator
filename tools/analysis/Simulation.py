import numpy as np
import h5py

import Stimulus
import Integrator
import Cell

class Simulation:
    """
    Class for a single simulation

    """
    def __init__(self, h5_file):
        self.cell_types = []
        self.simulation_file = h5_file
        ########################## Read file ###################################
        for item in h5_file.keys():
            if(item == "stimulus"):
                stim_group = h5_file.get("/" + str(item))
                self.stimulus = Stimulus.Stimulus(stim_group)
            elif(item == "integrator"):
                integrator_group = h5_file.get("/" + str(item))
                self.integrator = Integrator.Integrator(integrator_group)
            else:
                cell_group = h5_file.get("/" + str(item))
                setattr(self, item, Cell.Cell(cell_group))
                self.cell_types.append(item)
        ########################################################################
        self.num_cell_types  = len(self.cell_types)
        # self.zeroOutsmallValues()
        # self.normalize()

    def zeroOutsmallValues(self):
        for cell in self.cell_types:
            cell_type = getattr(self, cell)
            cell_type.zeroOutsmallValues()

    def normalize(self):
        for cell in self.cell_types:
            cell_type = getattr(self, cell)
            cell_type.normalize()

    def singleCellTemporalResponse(self, cellType, idx=0 , idy=0):
        response = getattr(self, cellType).response["spatioTemporal"][:,idx, idy]
        return response

    def singleCellFreqResponse(self, cellType, idx=0 , idy=0):
        FreqResponse = getattr(self, cellType).response["fourierTransform"][:,idx, idy]
        return FreqResponse

    def temporalImpulseResponse(self, cellType, idx=0 , idy=0):
        impulseResponse = getattr(self, cellType).impulseResponse["spatioTemporal"][:,idx, idy]
        return impulseResponse

    def spikeTrain(self, cellType, idx=0 , idy=0, num_trails=2):
        response = self.singleCellTemporalResponse(cellType, idx , idy)
        spike_times = [[] for i in range(num_trails)]

        dt = self.integrator.temporalResolution
        for k in range(num_trails):
            for i in range(self.integrator.nPointsTemporal):
                r = np.random.uniform(0,1)
                if(response[i] * dt > r):
                    spike_times[k].append(i*dt)
            spike_times[k]=np.array(spike_times[k])
        return spike_times

if __name__ == "__main__":
    print "---Class for single lgn-simulator experiments---"
