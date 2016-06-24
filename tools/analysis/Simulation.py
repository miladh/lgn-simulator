import numpy as np
import h5py
import yaml

import Stimulus
import Integrator
import Neuron

class Simulation:
    """
    Class for a single simulation
    """
    def __init__(self, config_file, h5_file):
        self.cell_types = []
        self.config_file = config_file
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
                setattr(self, item, Neuron.Neuron(cell_group))
                self.cell_types.append(item)
        ########################################################################
        self.num_cell_types  = len(self.cell_types)


    def get_attribute(self, key):
        from operator import attrgetter
        try:
            return attrgetter(key)(self)
        except AttributeError:
            return self.__get_setting(key)

    def __get_setting(self, path):
        setting_path = path.split('.')
        try:
            with open(self.config_file, 'r') as f:
                config_data = yaml.load(f)
                key =  config_data[setting_path[0]]
                for node in setting_path[1:]:
                    key =  key[node]
                return key
        except KeyError:
            raise KeyError("key not found: " + path + ", options: " + str(key))


    def spike_train(self, cell_type, rc=[0.0, 0.0], num_trails=2):
        response = getattr(self, cell_type).t_resp(rc)
        spike_times = [[] for i in range(num_trails)]

        dt = self.integrator.dt
        for k in range(num_trails):
            for i in range(self.integrator.Nt):
                r = np.random.uniform(0,1)
                if(response[i] * dt > r):
                    spike_times[k].append(i*dt)
            spike_times[k]=np.array(spike_times[k])
        return spike_times


    def cross_correlogram(self, spike_times_A, spike_times_B):
        cross_correlogram = []
        num_trails = min(len(spike_times_A),len(spike_times_B))

        for k in range(num_trails):
            for ta in spike_times_A[k]:
                for tb in spike_times_B[k]:
                    cross_correlogram.append(ta - tb)
        return cross_correlogram


if __name__ == "__main__":
    print "---Class for single lgn-simulator experiments---"
