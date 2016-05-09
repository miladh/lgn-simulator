import numpy as np
import yaml
import os, os.path

def simulation_extractor(sims, key, value):
    """
    extracts simulations with a specific run setting parameter

    Parameters
    ----------
    sims : list
        list of Simulation objects
    key: str
        attribute/setting which the extraction is based on

    Returns
    -------
    list
        list of Simulation objects with specific run attribute/setting parameter

    """
    from analysis.Simulation import Simulation
    extracted_sims = []
    for sim in sims:
        p = sim.get_attribute(key)
        if(p==value):
            extracted_sims.append(sim)

    return extracted_sims




def extract_unique_simulation_param(sims, key):
    """
    extracts simulations with a specific run setting parameter

    Parameters
    ----------
    sims : list
        list of Simulation objects
    key: str
        name of the attribute/setting parameter

    Returns
    -------
    array
        array of unique instances of the given paramter

    """
    from analysis.Simulation import Simulation
    key_values = set()
    for sim in sims:
        key_values.add(sim.get_attribute(key))

    return np.array(sorted(list(key_values)))



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
    from analysis.Simulation import Simulation

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
