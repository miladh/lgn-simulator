import numpy as np
import yaml
import os, os.path

from analysis.Simulation import Simulation

def simulation_extractor(sims, attr, value):
    """
    extracts simulations where attribute
    attr=value

    Parameters
    ----------
    sims : list
        list of Simulation objects
    attr: str
        name of the attribute

    Returns
    -------
    list
        list of Simulation objects where attribute
        attr=value

    """
    extracted_sims = []
    for sim in sims:
        p = sim.get_attribute(attr)
        if(p==value):
            extracted_sims.append(sim)

    return extracted_sims




def extract_unique_simulation_attrs(sims, attr):
    """
    extracts unique values for attribute attr
    from a list of simulations

    Parameters
    ----------
    sims : list
        list of Simulation objects
    attr: str
        name of the attribute

    Returns
    -------
    array
        array of unique attr values in the sims

    """

    attr_values = set()
    for sim in sims:
        attr_values.add(sim.get_attribute(attr))

    return np.array(sorted(list(attr_values)))



def get_simulations(config_file):
    """
    returns list of simulations

    Parameters
    ----------
    config_file : str
        path to config file with simulation ids

    Returns
    -------
    list
        list of Simulation objects.

    """
    import h5py

    with open(config_file, 'r') as stream:
        config_data = yaml.load(stream)
        data_path = config_data["data_path"]
        data_ids  = config_data["simulation_ids"]

    print "Data path: ", data_path
    print "Reading simulation_ids:\n", data_ids

    sims=[]
    for data_id in data_ids:
        setting_file = os.path.join(data_path, data_id,"_"+data_id +".yaml")
        data_file = os.path.join(data_path, data_id, data_id +".h5")
        f = h5py.File(data_file, "r")
        sims.append(Simulation(setting_file, f))

    return sims
